#!/bin/bash

# Query the ENA Portal API for the sample, study, country, location and
# checklist of each assembly (ERZ), 500 assemblies per request.

# Needs internet access.

# Resumable: each batch's result is saved separately and skipped on re-runs.
# A failed request (HTTP error after curl's own retries) leaves no file, so
# re-running retries only the failed batches.

# Defaults:
#   erz_list.txt  dbcan_1071/erz_list.txt   (from make_protein_erz.py)
#   out.tsv       dbcan_1071/erz_geo.tsv
#
# Output: no header; columns are
#   analysis_accession, sample_accession, study_accession, country, location, checklist
# Batch files are kept in <out>_batches/ for resuming; can be deleted once the run is complete.

set -uo pipefail

PROJECT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT"

LIST="${1:-dbcan_1071/erz_list.txt}"
OUT="${2:-dbcan_1071/erz_geo.tsv}"
DIR="${OUT%.tsv}_batches"
URL="https://www.ebi.ac.uk/ena/portal/api/search"
FIELDS="sample_accession,study_accession,country,location,checklist"

command -v curl >/dev/null || { echo "ERROR: curl not found" >&2; exit 1; }
test -s "$LIST" || { echo "ERROR: $LIST is missing or empty" >&2; exit 1; }
mkdir -p "$DIR"
# Same list -> same batch files, so re-splitting on a re-run is harmless.
split -d -a 5 -l 500 "$LIST" "$DIR/batch_"

n=$(ls "$DIR"/batch_????? | wc -l)
done_n=0; failed=0; i=0
for f in "$DIR"/batch_?????; do
    i=$((i + 1))
    if [ -s "$f.tsv" ]; then done_n=$((done_n + 1)); continue; fi
    code=$(curl -sS --retry 5 --retry-delay 30 -X POST "$URL" \
            --data-urlencode result=analysis --data-urlencode format=tsv \
            --data-urlencode "fields=$FIELDS" \
            --data-urlencode "includeAccessions=$(paste -sd, "$f")" \
            -o "$f.tmp" -w '%{http_code}')
    rc=$?
    if [ "$rc" -eq 0 ] && [ "$code" = "200" ]; then
        mv "$f.tmp" "$f.tsv"
        done_n=$((done_n + 1))
    else
        failed=$((failed + 1))
        echo "$(date '+%F %T')  $(basename "$f") failed: curl exit $rc, HTTP $code" >&2
        [ -s "$f.tmp" ] && head -c 300 "$f.tmp" >&2 && echo >&2
        rm -f "$f.tmp"
    fi
    [ $((i % 50)) -eq 0 ] && echo "$(date '+%F %T')  $i/$n batches processed"
    sleep 1
done

# Each batch result starts with a header line; drop it when combining.
tail -q -n +2 "$DIR"/batch_*.tsv > "$OUT"
echo "$(date '+%F %T')  $done_n/$n batches done; $failed failed"
echo "$(date '+%F %T')  returned $(cut -f1 "$OUT" | sort -u | wc -l) of $(wc -l < "$LIST") assemblies -> $OUT"
if [ "$failed" -gt 0 ]; then
    echo "re-run the same command to retry the failed batches" >&2
    exit 2
fi
