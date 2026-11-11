#!/bin/bash
#SBATCH --job-name=cazy_pfam
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=05:00:00
#SBATCH --output=logs/%x_%j.out

set -eo pipefail

for tool in hmmsearch seqkit; do
    command -v "$tool" >/dev/null || { echo "Missing prerequisite: $tool (see README)"; exit 1; }
done

THREADS=${SLURM_CPUS_PER_TASK:-32}
NCHUNKS=16
PER_CHUNK=$(( THREADS / NCHUNKS ))

log()     { echo "[$(date '+%F %T')] $*"; }
done_ok() { [ -s "$1" ] && tail -n 1 "$1" | grep -q '\[ok\]'; }


for f in cazy_nr.fa Pfam-A.hmm; do
    [ -e "$f" ] || { echo "Missing input: $f"; exit 1; }
done
log "Unique sequences: $(grep -c '>' cazy_nr.fa)"

[ -d chunks ] || seqkit split2 -p "$NCHUNKS" cazy_nr.fa -O chunks
todo=$(for c in chunks/*.fa; do done_ok "$c.pfam.domtbl" || echo "$c"; done)
if [ -n "$todo" ]; then
    log "Running Pfam hmmsearch on $(echo "$todo" | wc -l) chunk(s)"
    printf '%s\n' $todo | xargs -P "$NCHUNKS" -I{} \
        hmmsearch --cut_ga --cpu "$PER_CHUNK" --noali -o /dev/null \
        --domtblout {}.pfam.domtbl Pfam-A.hmm {}
fi
for c in chunks/*.fa; do
    done_ok "$c.pfam.domtbl" || { echo "Incomplete: $c"; exit 1; }
done
cat chunks/*.pfam.domtbl > pfam.domtbl.tmp
mv pfam.domtbl.tmp pfam.domtbl
rm -rf chunks
log "All done"
