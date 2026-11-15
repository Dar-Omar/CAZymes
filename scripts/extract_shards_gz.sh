#!/bin/bash
# Extract MGnify proteins whose IDs are listed in the .npy file
# from the shards of MGnify Proteins, keeping full-length sequences.

# For each shard: download (resumable, with retries) -> filter -> append to the
# gzipped output -> delete the shard.
#
# Resumable: finished shards are listed in <out>.done and skipped on re-run.
# An append interrupted before its shard was marked done is undone on re-run,
# so no shard is duplicated in the output.
#
# Requires internet access to ftp.ebi.ac.uk. Run it inside a screen session.
# prefetch_shards.sh can download shards ahead in parallel.
# Shards are locked atomically while downloading, so the two scripts don't
# download the same shard at the same time.
#
# Usage (from the project directory):
#   scripts/extract_shards_gz.sh [ids.npy] [out.faa.gz] [any]
#     ids.npy      default: cazy1071_ids.npy
#     out.faa.gz   default: cazy1071_fl.faa.gz
#     any          also keep partial sequences (default: FL=1 only)

set -uo pipefail

IDS="${1:-cazy1071_ids.npy}"
OUT="${2:-cazy1071_fl.faa.gz}"
KEEP="${3:-full}"

MGNIFY_RELEASE="${MGNIFY_RELEASE:-2026_07}"
BASE="${BASE:-https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/$MGNIFY_RELEASE}"
WORK="${WORK:-${TMPDIR:-$HOME/tmp}/mgy_shards}"
N_SHARDS="${N_SHARDS:-48}"
MAX_ATTEMPTS="${MAX_ATTEMPTS:-200}"
RETRY_WAIT="${RETRY_WAIT:-30}"
STALE_MIN="${STALE_MIN:-120}"

DONE="${OUT}.done"
COUNTS="${OUT}.counts"
PENDING="${OUT}.pending"

log()  { echo "$(date '+%F %T')  $*"; }
warn() { echo "$(date '+%F %T')  $*" >&2; }
sum_counts() { awk -F'\t' '{c[$1]=$2} END{s=0; for(k in c) s+=c[k]; print s}' "$COUNTS"; }

# checks
test -s "$IDS" || { warn "ERROR: $IDS is missing or empty"; exit 1; }
for tool in wget curl python3; do
    command -v "$tool" >/dev/null || { warn "ERROR: missing prerequisite: $tool"; exit 1; }
done
python3 -c "import numpy" 2>/dev/null || { warn "ERROR: python3 cannot import numpy"; exit 1; }
n_ids=$(python3 -c "import numpy as np, sys; print(len(np.load(sys.argv[1], mmap_mode='r')))" "$IDS") \
    || { warn "ERROR: cannot read $IDS as a .npy file"; exit 1; }

mkdir -p "$WORK"
touch "$DONE" "$COUNTS" "$OUT"

# recover from an interrupted previous run
# Locks this script left behind (it writes "extract" into its own locks).
for l in "$WORK"/*.lock; do
    [ -e "$l" ] && grep -q '^extract' "$l" && rm -f "$l" && log "removed own stale lock $(basename "$l")"
done
# An append that was not followed by marking its shard done: cut it off again.
if [ -s "$PENDING" ]; then
    read -r pname psize < "$PENDING"
    if ! grep -qxF "$pname" "$DONE"; then
        truncate -s "$psize" "$OUT"
        log "undid interrupted append of $pname"
    fi
    rm -f "$PENDING"
fi

CUR_LOCK=""
trap '[ -n "$CUR_LOCK" ] && rm -f "$CUR_LOCK"' EXIT

# is EBI reachable?
log "checking EBI availability ..."
probe_ok=0
for _ in 1 2 3; do
    code=$(curl -s -o /dev/null -m 20 -w '%{http_code}' "$BASE/md5sum.txt" || true)
    log "  probe: $code"
    [ "$code" = "200" ] && probe_ok=$((probe_ok + 1))
    sleep 3
done
if [ "$probe_ok" -eq 0 ]; then
    warn "ERROR: EBI returned no successful response in 3 probes; try again later."
    exit 1
fi
[ "$probe_ok" -lt 3 ] && warn "WARNING: only $probe_ok/3 probes succeeded; expect retries."

log "release: $MGNIFY_RELEASE   ids: $n_ids   keep: $KEEP"
log "scratch: $WORK"
log "retry policy: up to $MAX_ATTEMPTS attempts, ${RETRY_WAIT}s apart; stale lock after ${STALE_MIN}m"

# main loop
for i in $(seq 1 "$N_SHARDS"); do
    NAME="mgy_proteins_${i}.fasta.gz"
    SHARD="$WORK/$NAME"
    LOCK="$SHARD.lock"
    tag="[$i/$N_SHARDS]"

    if grep -qxF "$NAME" "$DONE"; then
        log "$tag $NAME already done, skipping"
        continue
    fi

    # Claim the shard. The lock is created atomically, so it fails while
    # another process (prefetch_shards.sh) is downloading the shard: then wait.
    # If the file stops growing for STALE_MIN minutes, that lock is stale.
    last=-1; idle=0; announced=0
    until ( set -o noclobber; echo "extract $$" > "$LOCK" ) 2>/dev/null; do
        if [ "$announced" -eq 0 ]; then
            log "$tag $NAME is being downloaded by another process, waiting"
            announced=1
        fi
        size=$(stat -c%s "$SHARD" 2>/dev/null || echo 0)
        if [ "$size" = "$last" ]; then idle=$((idle + 1)); else idle=0; last=$size; fi
        if [ "$idle" -ge "$STALE_MIN" ]; then
            warn "$tag no download progress for ${STALE_MIN}m; removing stale lock"
            rm -f "$LOCK"
            idle=0
            continue
        fi
        sleep 60
    done
    CUR_LOCK="$LOCK"

    # Download (resumes a partial file; returns at once if already complete).
    log "$tag downloading $NAME"
    attempt=0
    until wget -c --tries=20 --waitretry=30 --retry-connrefused \
               --read-timeout=120 -O "$SHARD" "$BASE/$NAME"; do
        attempt=$((attempt + 1))
        if [ "$attempt" -ge "$MAX_ATTEMPTS" ]; then
            warn "$tag giving up after $attempt attempts"
            rm -f "$SHARD" "$LOCK"; CUR_LOCK=""
            continue 2
        fi
        warn "$tag attempt $attempt failed, waiting ${RETRY_WAIT}s"
        sleep "$RETRY_WAIT"
    done
    rm -f "$LOCK"; CUR_LOCK=""

    # Filter into a separate part file.
    log "$tag filtering ($(du -h "$SHARD" | cut -f1))"
    PART="${OUT}.part.${i}"
    rm -f "$PART"
    n=$(python3 - "$IDS" "$SHARD" "$PART" "$KEEP" <<'PYEOF'
import gzip, sys
import numpy as np

ids_path, shard_path, out_path, keep = sys.argv[1:5]
full_only = (keep != "any")
wanted = np.load(ids_path, mmap_mode="r")
n_wanted = len(wanted)
BATCH = 200_000

written = 0
batch_ids, batch_recs = [], []

def flush(out):
    global written
    if not batch_ids:
        return
    ids = np.asarray(batch_ids, dtype=np.int64)
    idx = np.searchsorted(wanted, ids)
    idx[idx >= n_wanted] = n_wanted - 1
    hit = wanted[idx] == ids
    for h, rec in zip(hit, batch_recs):
        if h:
            out.write(rec)
            written += 1
    batch_ids.clear()
    batch_recs.clear()

def take(header, seq, out):
    acc = header[1:].split(None, 1)[0]
    try:
        pid = int(acc[4:]) if acc.upper().startswith("MGYP") else int(acc)
    except ValueError:
        return
    batch_ids.append(pid)
    batch_recs.append(header + "".join(seq))
    if len(batch_ids) >= BATCH:
        flush(out)

with gzip.open(shard_path, "rt") as stream, \
     gzip.open(out_path, "wt", compresslevel=3) as out:
    header, seq, skip = None, [], True
    for line in stream:
        if line.startswith(">"):
            if not skip:
                take(header, seq, out)
            header, seq = line, []
            skip = full_only and "FL=1" not in line
        elif not skip:
            seq.append(line)
    if header is not None and not skip:
        take(header, seq, out)
    flush(out)
print(written)
PYEOF
)
    rc=$?
    if [ $rc -ne 0 ]; then
        warn "$tag FILTER FAILED (exit $rc), discarding shard and partial output"
        rm -f "$PART" "$SHARD"
        continue
    fi

    size_before=$(stat -c%s "$OUT")
    printf '%s\t%s\n' "$NAME" "$size_before" > "$PENDING"
    if ! cat "$PART" >> "$OUT"; then
        truncate -s "$size_before" "$OUT"
        rm -f "$PART" "$PENDING"
        warn "$tag APPEND FAILED (disk full?); output restored, shard kept for re-run"
        continue
    fi
    printf '%s\t%s\n' "$NAME" "$n" >> "$COUNTS"
    echo "$NAME" >> "$DONE"
    rm -f "$PENDING" "$PART" "$SHARD"
    log "$tag matched $n   running total $(sum_counts)"
done

remaining=$(( N_SHARDS - $(grep -c . "$DONE") ))
expected=$(sum_counts)
log "counting sequences in $OUT (one full pass) ..."
actual=$(zcat "$OUT" | grep -c '^>')
log "complete: $actual sequences in $OUT (expected $expected)"
if [ "$actual" -ne "$expected" ]; then
    warn "WARNING: $((actual - expected)) sequence(s) differ from the per-shard counts; check $OUT"
fi
if [ "$remaining" -gt 0 ]; then
    warn "$remaining shard(s) still unfinished - re-run to retry"
    exit 2
fi
rmdir "$WORK" 2>/dev/null
exit 0
