#!/bin/bash
# Download MGnify Proteins shards in parallel ahead of extract_shards_gz.sh,
# so the extractor rarely has to wait for a download.
# N parallel downloads are about N times faster than one.
# Shards listed in the extract <out>.done are skipped.
# Run inside a screen session.

set -uo pipefail

OUT="${1:-cazy1071_fl.faa.gz}"
DONE="${OUT}.done"

MGNIFY_RELEASE="${MGNIFY_RELEASE:-2026_07}"
BASE="${BASE:-https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/$MGNIFY_RELEASE}"
WORK="${WORK:-${TMPDIR:-$HOME/tmp}/mgy_shards}"
N_SHARDS="${N_SHARDS:-48}"
FIRST="${FIRST:-1}"
LAST="${LAST:-$N_SHARDS}"
JOBS="${JOBS:-5}"
MAX_STAGED="${MAX_STAGED:-3}"
MIN_FREE_GB="${MIN_FREE_GB:-60}"
MAX_ATTEMPTS="${MAX_ATTEMPTS:-100}"

log() { echo "$(date '+%F %T')  $*"; }

for tool in wget curl timeout; do
    command -v "$tool" >/dev/null || { log "ERROR: missing prerequisite: $tool" >&2; exit 1; }
done
mkdir -p "$WORK"

free_gb() { df -BG --output=avail "$WORK" | tail -1 | tr -dc '0-9'; }

staged() {
    # complete shards waiting for the extractor (shard files without a lock)
    local n=0 f
    for f in "$WORK"/mgy_proteins_*.fasta.gz; do
        [ -e "$f" ] || continue
        [ -e "$f.lock" ] && continue
        n=$((n + 1))
    done
    echo "$n"
}

is_done() { [ -f "$DONE" ] && grep -qxF "$1" "$DONE"; }

# Create a lock only if none exists (atomic), writing a token into it.
claim() { ( set -o noclobber; echo "$2" > "$1" ) 2>/dev/null; }
owns()  { [ "$(cat "$1" 2>/dev/null)" = "$2" ]; }

fetch_one() {
    local i="$1" token="$2"
    local name="mgy_proteins_${i}.fasta.gz"
    local file="$WORK/$name" lock="$WORK/$name.lock" attempt=0
    while :; do
        # stop if the extractor declared this lock stale and took the shard over
        if ! owns "$lock" "$token"; then
            log "  [$i] lock taken over by another process, stopping"
            return 1
        fi
        timeout 7200 wget -c -q --read-timeout=120 -O "$file" "$BASE/$name" && break
        attempt=$((attempt + 1))
        if [ "$attempt" -ge "$MAX_ATTEMPTS" ]; then
            log "  [$i] gave up after $attempt attempts (the extractor will retry it)"
            if owns "$lock" "$token"; then rm -f "$file" "$lock"; fi
            return 1
        fi
        sleep 30
    done
    owns "$lock" "$token" && rm -f "$lock"
    log "  [$i] complete ($(du -h "$file" | cut -f1))"
}

# stop running downloads and remove this run's locks.
cleanup() {
    trap - EXIT INT TERM
    local p l
    for p in $(jobs -p); do
        pkill -TERM -P "$p" 2>/dev/null
        kill "$p" 2>/dev/null
    done
    wait 2>/dev/null
    for l in "$WORK"/*.lock; do
        [ -e "$l" ] && grep -q "^prefetch $$ " "$l" && rm -f "$l"
    done
}
trap cleanup EXIT
trap 'exit 130' INT TERM

# Locks left by an earlier prefetcher that was killed (only one runs at a time).
for l in "$WORK"/*.lock; do
    [ -e "$l" ] && grep -q '^prefetch ' "$l" && rm -f "$l" && log "removed stale lock $(basename "$l")"
done

log "prefetching shards $FIRST-$LAST of release $MGNIFY_RELEASE with $JOBS streams"
log "scratch: $WORK   done list: $DONE"
log "limits: <= $MAX_STAGED staged shards, >= ${MIN_FREE_GB}G free"

running=0
for i in $(seq "$FIRST" "$LAST"); do
    name="mgy_proteins_${i}.fasta.gz"

    # Back-pressure: wait for the extractor to catch up, or for disk space.
    while [ "$(staged)" -ge "$MAX_STAGED" ] || [ "$(free_gb)" -lt "$MIN_FREE_GB" ]; do
        log "  waiting: $(staged) staged, $(free_gb)G free"
        sleep 300
    done

    if is_done "$name"; then
        log "  [$i] already extracted, skipping"
        continue
    fi
    remote_size=$(curl -sI -m 30 "$BASE/$name" | tr -d '\r' \
        | awk 'BEGIN{IGNORECASE=1} /^content-length:/ {print $2}')
    local_size=$(stat -c%s "$WORK/$name" 2>/dev/null || echo 0)
    if [ -n "$remote_size" ] && [ "$local_size" = "$remote_size" ]; then
        log "  [$i] already on disk, skipping"
        continue
    fi

    token="prefetch $$ $i"
    if ! claim "$WORK/$name.lock" "$token"; then
        log "  [$i] locked by another process, skipping"
        continue
    fi

    log "  [$i] starting"
    fetch_one "$i" "$token" &
    running=$((running + 1))
    if [ "$running" -ge "$JOBS" ]; then
        wait -n 2>/dev/null
        running=$((running - 1))
    fi
done

wait
log "prefetch finished"
