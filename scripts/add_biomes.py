#!/usr/bin/env python3
"""
Add MGnify biomes to dbCAN overview.

A new last column, Biomes, lists every biome the protein was observed in, with
the number of observations in that biome, most frequent first ("-" if none).
Biomes come from mgy_biomes.parquet (MGnify Proteins), where each protein's
biome is the biome of the assembly it was found in.

Defaults:
    --in      dbcan_1071/overview_2_substrates_group.tsv
    --out     dbcan_1071/overview_2_substrates_group_biomes.tsv
    --biomes  mgy_biomes.parquet
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
    ap.add_argument("--in", dest="infile", default="dbcan_1071/overview_2_substrates_group.tsv")
    ap.add_argument("--out", dest="outfile", default="dbcan_1071/overview_2_substrates_group_biomes.tsv")
    ap.add_argument("--biomes", default="mgy_biomes.parquet")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--memory", default="32GB")
    args = ap.parse_args()

    infile, outfile, biomes = (project_path(p) for p in (args.infile, args.outfile, args.biomes))
    for path in (infile, biomes):
        if not path.exists():
            sys.exit(f"missing input: {path}")
    with open(infile) as fh:
        if "Gene ID" not in next(fh).rstrip("\n").split("\t"):
            sys.exit(f"{infile} has no 'Gene ID' column")

    # private temp folder for DuckDB, removed on any exit
    tmpdir = tempfile.mkdtemp(prefix="duckdb_tmp_", dir=outfile.parent)
    atexit.register(shutil.rmtree, tmpdir, True)

    con = duckdb.connect()
    con.execute(f"SET threads={args.threads}")
    con.execute(f"SET memory_limit='{args.memory}'")
    con.execute(f"SET temp_directory='{tmpdir}'")
    con.execute("SET preserve_insertion_order=false")

    tmp_out = f"{outfile}.tmp"
    con.execute(f"""
        COPY (
            WITH ov AS (
                SELECT *,
                       CAST(substr("Gene ID", 5) AS BIGINT) AS protein_id
                FROM read_csv('{infile}', delim='\t', header=true, all_varchar=true)
            ),
            bio AS (
                SELECT protein_id,
                       string_agg(biome || ' (' || lineage_count || ')', '; '
                                  ORDER BY lineage_count DESC, biome) AS Biomes
                FROM read_parquet('{biomes}')
                WHERE protein_id IN (SELECT protein_id FROM ov)
                GROUP BY protein_id
            )
            SELECT ov.* EXCLUDE (protein_id), COALESCE(bio.Biomes, '-') AS Biomes
            FROM ov LEFT JOIN bio USING (protein_id)
            ORDER BY "Gene ID"
        ) TO '{tmp_out}' (HEADER, DELIMITER '\t')
    """)
    os.replace(tmp_out, outfile)

    n, missing = con.execute(f"""
        SELECT COUNT(*), COUNT(*) FILTER (WHERE Biomes = '-')
        FROM read_csv('{outfile}', delim='\t', header=true, all_varchar=true)
    """).fetchone()
    con.close()
    if missing:
        print(f"WARNING: {missing:,} of {n:,} proteins have no biome", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
