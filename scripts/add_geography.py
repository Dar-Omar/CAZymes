#!/usr/bin/env python3
"""
Add ENA sample geography to the CAZyme table.

Inputs
  protein_erz.tsv  protein_accession <tab> erz, with header
                   (one row per CAZyme x assembly; from make_protein_erz.py)
  erz_geo.tsv      ENA results, no header (from fetch_ena_geo.sh):
                   analysis_accession, sample_accession, study_accession,
                   country, location[, checklist]

Outputs
  --long     one row per (CAZyme, assembly) with that assembly's geography
  --out      the CAZyme table with these columns appended:
               assemblies        every ERZ containing the CAZyme, ";"-separated
               sample_accession  distinct samples
               study_accession   distinct studies
               country           distinct values with assembly counts, most
                                 frequent first: "Canada:Alberta (3); Japan (1)"
               location          same format, for coordinates
               checklist         distinct ENA checklists
               n_assemblies, n_with_country, n_with_location
  --missing  ERZs that ENA did not return

ENA missing-value placeholders ("missing", "not collected", ...) are treated
as empty. Empty text columns are written as "-", as in the rest of the table.
"""

import argparse
import atexit
import os
import shutil
import sys
import tempfile
from pathlib import Path

import duckdb

PROJECT = Path(__file__).resolve().parents[1]

PLACEHOLDER = (r"^(missing.*|not (applicable|collected|provided).*"
               r"|restricted access.*|unknown|n/?a|null|none|)$")
TEXT_COLS = ["assemblies", "sample_accession", "study_accession",
             "country", "location", "checklist"]
COUNT_COLS = ["n_assemblies", "n_with_country", "n_with_location"]


def project_path(p):
    """Relative paths are taken from the project directory."""
    p = Path(p).expanduser()
    return p if p.is_absolute() else PROJECT / p


def distinct(col):
    """Sorted distinct non-null values joined by ';', or '-' if none."""
    return (f"COALESCE(NULLIF(array_to_string("
            f"list_sort(list_distinct(list({col}))), ';'), ''), '-') AS {col}")


def counted(col):
    """Per CAZyme: 'value (n assemblies); ...', most frequent first."""
    return f"""
        SELECT "Gene ID",
               string_agg(v || ' (' || n || ')', '; ' ORDER BY n DESC, v) AS {col}
        FROM (SELECT "Gene ID", {col} AS v, COUNT(*) AS n
              FROM long WHERE {col} IS NOT NULL GROUP BY 1, 2)
        GROUP BY 1
    """


def pct(a, b):
    return f"{100 * a / b:.1f}%" if b else "n/a"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--overview", default="dbcan_1071/overview_2_substrates_group_biomes_category.tsv")
    ap.add_argument("--protein-erz", default="dbcan_1071/protein_erz.tsv")
    ap.add_argument("--erz-geo", default="dbcan_1071/erz_geo.tsv")
    ap.add_argument("--out", default="dbcan_1071/overview_2_substrates_group_biomes_category_geo.tsv")
    ap.add_argument("--long", default="dbcan_1071/protein_geo_long.tsv")
    ap.add_argument("--missing", default="dbcan_1071/erz_missing.txt")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--memory", default="48GB", help="DuckDB memory limit")
    args = ap.parse_args()

    overview, protein_erz, erz_geo, out, long_out, missing_out = (project_path(p) for p in (
        args.overview, args.protein_erz, args.erz_geo, args.out, args.long, args.missing))
    for path in (overview, protein_erz, erz_geo):
        if not path.exists():
            sys.exit(f"missing input: {path}")

    # private temp folder for DuckDB, removed on any exit
    tmpdir = tempfile.mkdtemp(prefix="duckdb_tmp_", dir=out.parent)
    atexit.register(shutil.rmtree, tmpdir, True)

    con = duckdb.connect()
    con.execute(f"SET threads={args.threads}")
    con.execute(f"SET memory_limit='{args.memory}'")
    con.execute(f"SET temp_directory='{tmpdir}'")
    con.execute("SET preserve_insertion_order=false")
    con.execute(f"""
        CREATE TEMP MACRO clean(x) AS
        CASE WHEN x IS NULL OR regexp_matches(lower(trim(x)), '{PLACEHOLDER}')
             THEN NULL ELSE trim(x) END
    """)

    # ENA results; null_padding accepts files with or without the checklist column.
    con.execute(f"""
        CREATE TEMP TABLE geo AS
        SELECT DISTINCT trim(assemblies)        AS assemblies,
                        clean(sample_accession) AS sample_accession,
                        clean(study_accession)  AS study_accession,
                        clean(country)          AS country,
                        clean(location)         AS location,
                        clean(checklist)        AS checklist
        FROM read_csv('{erz_geo}', delim='\t', header=false, null_padding=true,
                      columns={{'assemblies': 'VARCHAR', 'sample_accession': 'VARCHAR',
                               'study_accession': 'VARCHAR', 'country': 'VARCHAR',
                               'location': 'VARCHAR', 'checklist': 'VARCHAR'}})
        WHERE trim(assemblies) LIKE 'ERZ%'
    """)
    dup = con.execute("SELECT COUNT(*) - COUNT(DISTINCT assemblies) FROM geo").fetchone()[0]
    if dup:
        print(f"WARNING: {dup} ERZ(s) have conflicting rows in {erz_geo}; "
              f"assembly counts for their CAZymes will be inflated", file=sys.stderr)

    con.execute(f"""
        CREATE TEMP TABLE pe AS
        SELECT protein_accession, erz
        FROM read_csv('{protein_erz}', delim='\t', header=true, all_varchar=true)
    """)

    # ERZs ENA did not return
    con.execute("""
        CREATE TEMP TABLE missing AS
        SELECT DISTINCT erz FROM pe WHERE erz NOT IN (SELECT assemblies FROM geo)
    """)
    n_missing = con.execute("SELECT COUNT(*) FROM missing").fetchone()[0]

    # Long table: one row per (CAZyme, assembly)
    con.execute("""
        CREATE TEMP TABLE long AS
        SELECT p.protein_accession AS "Gene ID", p.erz AS assemblies,
               g.sample_accession, g.study_accession, g.country, g.location, g.checklist
        FROM pe p LEFT JOIN geo g ON p.erz = g.assemblies
    """)

    # One row per CAZyme
    con.execute(f"""
        CREATE TEMP TABLE per_protein AS
        SELECT b."Gene ID", b.assemblies, b.sample_accession, b.study_accession,
               COALESCE(c.country, '-')  AS country,
               COALESCE(l.location, '-') AS location,
               b.checklist, b.n_assemblies, b.n_with_country, b.n_with_location
        FROM (SELECT "Gene ID",
                     array_to_string(list_sort(list_distinct(list(assemblies))), ';')
                         AS assemblies,
                     {distinct('sample_accession')},
                     {distinct('study_accession')},
                     {distinct('checklist')},
                     COUNT(DISTINCT assemblies) AS n_assemblies,
                     COUNT(DISTINCT assemblies) FILTER (WHERE country IS NOT NULL)
                         AS n_with_country,
                     COUNT(DISTINCT assemblies) FILTER (WHERE location IS NOT NULL)
                         AS n_with_location
              FROM long GROUP BY 1) b
        LEFT JOIN ({counted('country')})  c USING ("Gene ID")
        LEFT JOIN ({counted('location')}) l USING ("Gene ID")
    """)

    # Write the three outputs, each to .tmp first
    text = ",\n".join(f"COALESCE(p.{c}, '-') AS {c}" for c in TEXT_COLS)
    nums = ",\n".join(f"COALESCE(p.{c}, 0) AS {c}" for c in COUNT_COLS)
    outputs = [
        (missing_out, "SELECT erz FROM missing ORDER BY 1", "HEADER false"),
        (long_out, "SELECT * FROM long ORDER BY 1, 2", "HEADER, DELIMITER '\t'"),
        (out, f"""
            SELECT o.*, {text}, {nums}
            FROM read_csv('{overview}', delim='\t', header=true, all_varchar=true) o
            LEFT JOIN per_protein p ON o."Gene ID" = p."Gene ID"
            ORDER BY o."Gene ID"
        """, "HEADER, DELIMITER '\t'"),
    ]
    for path, query, opts in outputs:
        tmp = f"{path}.tmp"
        con.execute(f"COPY ({query}) TO '{tmp}' ({opts})")
        os.replace(tmp, path)

    # Coverage summary
    n_erz, e_c, e_l = con.execute(
        "SELECT COUNT(*), COUNT(country), COUNT(location) FROM geo").fetchone()
    n_prot, p_c, p_l = con.execute("""
        SELECT COUNT(*),
               COUNT(*) FILTER (WHERE n_with_country > 0),
               COUNT(*) FILTER (WHERE n_with_location > 0)
        FROM per_protein
    """).fetchone()
    con.close()

    print(f"assemblies returned by ENA  : {n_erz:,}")
    print(f"  with country              : {e_c:,} ({pct(e_c, n_erz)})")
    print(f"  with location             : {e_l:,} ({pct(e_l, n_erz)})")
    print(f"CAZymes                     : {n_prot:,}")
    print(f"  >=1 assembly with country : {p_c:,} ({pct(p_c, n_prot)})")
    print(f"  >=1 assembly with location: {p_l:,} ({pct(p_l, n_prot)})")
    if n_missing:
        print(f"{n_missing:,} ERZ(s) not returned by ENA; listed in {missing_out}",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
