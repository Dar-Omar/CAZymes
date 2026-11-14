#!/usr/bin/env python3
"""
make_hits_and_ids.py

From a list of Pfam families, find every MGnify protein with a hit to any of
them and write:

  1. parquet file
       one row per distinct (protein_id, pfam_accession); reused if it already
       exists (--rebuild forces a new scan)
  2. ids
       each distinct protein_id once, sorted, as input for extract_shards_gz.sh
"""

import argparse
import atexit
import hashlib
import os
import shutil
import sys
import tempfile
import time

import duckdb


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def load_families(path):
    """Sorted Pfam integers from column 1 (PF00128.32 -> 128)."""
    fams = set()
    with open(path) as fh:
        for line in fh:
            acc = line.split("\t")[0].strip().split(".")[0]
            if acc.upper().startswith("PF"):
                try:
                    fams.add(int(acc[2:]))
                except ValueError:
                    pass
    return sorted(fams)


def main():
    slurm_mem = os.environ.get("SLURM_MEM_PER_NODE")
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--families-file", default="cazyme_pfam_min100.txt")
    ap.add_argument("--pfam-parquet", default="mgy_proteins_pfam.parquet")
    ap.add_argument("--cache-dir", default="cache")
    ap.add_argument("--ids", default="cazy1071_ids.txt")
    ap.add_argument("--rebuild", action="store_true",
                    help="redo the scan even if the hits file exists")
    ap.add_argument("--threads", type=int,
                    default=int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count())))
    ap.add_argument("--memory",
                    default=f"{int(int(slurm_mem) * 0.8)}MB" if slurm_mem else "128GB",
                    help="DuckDB memory limit (default: 80%% of the SLURM allocation)")
    ap.add_argument("--temp", default=".",
                    help="parent directory for DuckDB temporary files")
    args = ap.parse_args()

    if not os.path.exists(args.families_file):
        sys.exit(f"missing input: {args.families_file}")
    families = load_families(args.families_file)
    if not families:
        sys.exit(f"no Pfam accessions found in {args.families_file}")
    digest = hashlib.md5(",".join(map(str, families)).encode()).hexdigest()[:10]
    hits_path = os.path.join(args.cache_dir, f"hits_cazy{digest}_n{len(families)}.parquet")
    os.makedirs(args.cache_dir, exist_ok=True)

    tmpdir = tempfile.mkdtemp(prefix="duckdb_tmp_", dir=args.temp)
    atexit.register(shutil.rmtree, tmpdir, True)

    log(f"{len(families)} Pfam families from {args.families_file}")
    log(f"threads {args.threads}, memory limit {args.memory}")

    con = duckdb.connect()
    con.execute(f"SET threads={args.threads}")
    con.execute(f"SET memory_limit='{args.memory}'")
    con.execute(f"SET temp_directory='{tmpdir}'")
    con.execute("SET preserve_insertion_order=false")

    # hits parquet: reused if already built
    if os.path.exists(hits_path) and not args.rebuild:
        log(f"reusing {hits_path}")
    else:
        if not os.path.exists(args.pfam_parquet):
            sys.exit(f"missing input: {args.pfam_parquet}")
        log(f"scanning {args.pfam_parquet} (slow) ...")
        t0 = time.time()
        fam_list = ",".join(map(str, families))
        tmp = hits_path + ".tmp"
        con.execute(f"""
            COPY (SELECT DISTINCT protein_id, pfam_accession
                  FROM read_parquet('{args.pfam_parquet}')
                  WHERE pfam_accession IN ({fam_list}))
            TO '{tmp}' (FORMAT parquet)
        """)
        os.replace(tmp, hits_path)
        log(f"wrote {hits_path} in {time.time() - t0:,.0f}s")

    pairs, distinct = con.execute(f"""
        SELECT COUNT(*), COUNT(DISTINCT protein_id) FROM read_parquet('{hits_path}')
    """).fetchone()
    present = {r[0] for r in con.execute(
        f"SELECT DISTINCT pfam_accession FROM read_parquet('{hits_path}')").fetchall()}
    absent = [f"PF{f:05d}" for f in families if f not in present]

    print(f"(protein, family) pairs : {pairs:,}")
    print(f"distinct proteins       : {distinct:,}")
    print(f"families with hits      : {len(present):,} of {len(families):,}")
    if absent:
        print(f"families with no hits   : {', '.join(absent[:20])}"
              f"{' ...' if len(absent) > 20 else ''}")

    # distinct protein IDs, sorted
    tmp_ids = args.ids + ".tmp"
    con.execute(f"""
        COPY (SELECT DISTINCT protein_id FROM read_parquet('{hits_path}') ORDER BY protein_id)
        TO '{tmp_ids}' (FORMAT csv, HEADER false)
    """)
    os.replace(tmp_ids, args.ids)
    print(f"wrote {args.ids}: {distinct:,} IDs, {os.path.getsize(args.ids) / 1e9:.2f} GB")
    con.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
