#!/usr/bin/env python3
"""
Map each CAZyme to the MGnify assemblies (ERZ) it was found in, using
mgy_assemblies.tsv.gz, and count the distinct assemblies per substrate group.

Outputs (in dbcan_1071/ by default)
  protein_erz.tsv              one row per protein x assembly
  erz_list.txt                 distinct ERZ accessions, one per line
  assemblies_per_group.tsv     Group, CAZymes, Assemblies[, Pct_of_MGnify_assemblies]
                               (a CAZyme counts in every group it belongs to; an
                               assembly counts once per group however many of the
                               group's CAZymes it contains; the row "All groups"
                               counts assemblies with any substrate-annotated CAZyme)

--total-assemblies is the number of assemblies in the MGnify Proteins release
(51,274 in 2026_07, from the release's STATS.md); it is used for the percentage
of all MGnify assemblies that contain each group. Set to 0 to omit that column.
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


def project_path(p):
    """Relative paths are taken from the project directory."""
    p = Path(p).expanduser()
    return p if p.is_absolute() else PROJECT / p


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--overview", default="dbcan_1071/overview_2_substrates_group.tsv",
                    help="table with 'Gene ID' and 'Substrate_group' columns")
    ap.add_argument("--assemblies", default="mgy_assemblies.tsv.gz")
    ap.add_argument("--out", default="dbcan_1071/protein_erz.tsv")
    ap.add_argument("--erz-list", default="dbcan_1071/erz_list.txt")
    ap.add_argument("--groups-out", default="dbcan_1071/assemblies_per_group.tsv")
    ap.add_argument("--total-assemblies", type=int, default=51274,
                    help="assemblies in the MGnify release (2026_07: 51,274), for percentages; 0 to omit")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--memory", default="48GB", help="DuckDB memory limit; keep below the job's --mem")
    args = ap.parse_args()

    overview, assemblies, out, erz_list, groups_out = (project_path(p) for p in (
        args.overview, args.assemblies, args.out, args.erz_list, args.groups_out))
    for path in (overview, assemblies):
        if not path.exists():
            sys.exit(f"missing input: {path}")
    with open(overview) as fh:
        header = next(fh).rstrip("\n").split("\t")
    for col in ("Gene ID", "Substrate_group"):
        if col not in header:
            sys.exit(f"{overview} has no '{col}' column")

    # private temp folder for DuckDB, removed on any exit
    tmpdir = tempfile.mkdtemp(prefix="duckdb_tmp_", dir=out.parent)
    atexit.register(shutil.rmtree, tmpdir, True)

    con = duckdb.connect()
    con.execute(f"SET threads={args.threads}")
    con.execute(f"SET memory_limit='{args.memory}'")
    con.execute(f"SET temp_directory='{tmpdir}'")
    con.execute("SET preserve_insertion_order=false")

    # 1) CAZyme -> group (one row per CAZyme and group)
    con.execute(f"""
        CREATE TEMP TABLE pg AS
        SELECT DISTINCT protein_accession, trim(grp) AS grp
        FROM (SELECT "Gene ID" AS protein_accession,
                     unnest(string_split("Substrate_group", ';')) AS grp
              FROM read_csv('{overview}', delim='\t', header=true, all_varchar=true))
        WHERE trim(grp) NOT IN ('', '-')
    """)

    # 2) CAZyme -> assembly (the slow scan of all MGnify proteins)
    con.execute(f"""
        CREATE TEMP TABLE pe AS
        SELECT protein_accession,
               unnest(string_split(assembly_accessions, ';')) AS erz
        FROM read_csv('{assemblies}', delim='\t', header=false,
                      columns={{'protein_accession': 'VARCHAR',
                               'assembly_accessions': 'VARCHAR'}})
        WHERE protein_accession IN (SELECT protein_accession FROM pg)
    """)

    for path, query in ((out, "SELECT * FROM pe ORDER BY 1, 2"),
                        (erz_list, "SELECT DISTINCT erz FROM pe ORDER BY 1")):
        tmp = f"{path}.tmp"
        header_opt = "HEADER" if path == out else "HEADER false"
        con.execute(f"COPY ({query}) TO '{tmp}' ({header_opt}, DELIMITER '\t')")
        os.replace(tmp, path)

    # 3) distinct assemblies per group
    pct = (f", round(100.0 * Assemblies / {args.total_assemblies}, 2) AS Pct_of_MGnify_assemblies"
           if args.total_assemblies else "")
    tmp = f"{groups_out}.tmp"
    con.execute(f"""
        COPY (
            SELECT *{pct} FROM (
                SELECT grp AS "Group",
                       COUNT(DISTINCT pg.protein_accession) AS CAZymes,
                       COUNT(DISTINCT erz) AS Assemblies
                FROM pg JOIN pe USING (protein_accession)
                GROUP BY grp
                UNION ALL
                SELECT 'All groups', COUNT(DISTINCT protein_accession), COUNT(DISTINCT erz)
                FROM pe
            )
            ORDER BY ("Group" = 'All groups'), Assemblies DESC
        ) TO '{tmp}' (HEADER, DELIMITER '\t')
    """)
    os.replace(tmp, groups_out)

    n_caz, n_with = con.execute(
        "SELECT (SELECT COUNT(DISTINCT protein_accession) FROM pg), "
        "(SELECT COUNT(DISTINCT protein_accession) FROM pe)").fetchone()
    con.close()
    if n_with < n_caz:
        print(f"WARNING: {n_caz - n_with:,} of {n_caz:,} CAZymes have no assembly in {assemblies}",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
