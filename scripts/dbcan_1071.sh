#!/bin/bash
#SBATCH --job-name=dbcan_1071
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=05:00:00
#SBATCH --array=1-300%20
#SBATCH --output=logs/dbcan_1071_%A_%a.out

# dbCAN CAZyme annotation of the full-length candidate sequences, run as a job
# array with one 500,000-sequence chunk per task (at most 20 tasks at once).
#
# Before submitting run:
#   seqkit split2 -s 500000 -O chunks_1071 cazy1071_fl.faa.gz
#   mkdir -p logs
# The array size must equal the number of chunks, override at submission if needed

# Output: dbcan_1071/chunk_NNN/ .done marks a finished chunk.
# Resubmitting skips finished chunks and reruns failed or unfinished ones.
#
# Settings (environment variables defaults):
#   DBCAN_DB=$HOME/dbcan_db     dbCAN database directory
#   CHUNK_DIR=chunks_1071       input chunks
#   OUT_DIR=dbcan_1071          output directory

set -euo pipefail
shopt -s nullglob

command -v run_dbcan >/dev/null || { echo "Missing prerequisite: run_dbcan" >&2; exit 1; }

DB="${DBCAN_DB:-$HOME/dbcan_db}"
CHUNK_DIR="${CHUNK_DIR:-chunks_1071}"
OUT_DIR="${OUT_DIR:-dbcan_1071}"

FILES=("$CHUNK_DIR"/*.faa.gz)
N=${#FILES[@]}
ID=$SLURM_ARRAY_TASK_ID

[ "$N" -gt 0 ] || { echo "no chunk files in $CHUNK_DIR/" >&2; exit 1; }
if [ "$ID" -gt "$N" ]; then
    echo "task $ID: only $N chunk files exist, nothing to do"
    exit 0
fi
for f in CAZy.dmnd dbCAN.hmm dbCAN-sub.hmm fam-substrate-mapping.tsv; do
    [ -s "$DB/$f" ] || { echo "missing database file: $DB/$f" >&2; exit 1; }
done

IN=${FILES[$((ID - 1))]}
OUT=$OUT_DIR/chunk_$(printf '%03d' "$ID")

if [ -e "$OUT/.done" ]; then
    echo "$OUT already finished, skipping"
    exit 0
fi
rm -rf "$OUT"
mkdir -p "$OUT_DIR"

echo "task $ID of $N: $IN -> $OUT on $(hostname), started $(date '+%F %T')"

time run_dbcan CAZyme_annotation \
  --input_raw_data "$IN" \
  --mode protein \
  --db_dir "$DB" \
  --output_dir "$OUT" \
  --threads "$SLURM_CPUS_PER_TASK" \
  --large --large_dbsub \
  --batch_size 5000 --batch_size_dbsub 5000 \
  --log-level INFO

# save space: drop run_dbcan copy of the input, compress the large raw outputs
rm -f "$OUT/uniInput.faa"
for f in dbCANsub_hmm_raw.tsv diamond.out; do
    if [ -e "$OUT/$f" ]; then gzip -f "$OUT/$f"; fi
done

touch "$OUT/.done"
echo "finished $(date '+%F %T')"
