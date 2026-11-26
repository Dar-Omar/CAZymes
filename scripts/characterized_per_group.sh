#!/bin/bash

# Characterisation per substrate group: how many experimentally characterised
# CAZy proteins (and distinct EC numbers) support each group.

# For each characterised protein (CAZy header with an EC number), every
# family x EC combination in its header is looked up in dbCAN's
# fam-substrate-mapping.tsv (subfamilies reduced to family level). The substrates
# of the matching pairs are converted to groups (substrates.tsv, plus EXTRA_GROUPS
# below). A protein counts once per group. Groups cover synthesis as well as
# degradation, as in dbCAN's substrate mapping.
#
# Inputs (environment variables; defaults shown):
#   CHAR=cazy_characterized.faa                          characterised CAZy proteins
#   FAM_MAP=$HOME/dbcan_db/fam-substrate-mapping.tsv     dbCAN family/EC -> substrate
#   SUBS=dbcan_1071/substrates.tsv                       substrate -> group (Count, Substrate, Group)
#   OVERVIEW=dbcan_1071/overview_2_substrates_group.tsv  predicted CAZymes (Substrate_group column)
#   OUT=dbcan_1071/characterized_per_group.tsv
# Output columns (sorted by Characterized_sequences):
#   Group                    substrate group
#   CAZymes                  predicted CAZymes in the group (for reference)
#   Characterized_proteins   characterised CAZy proteins linked to the group (accessions)
#   Characterized_sequences  the same, counting identical sequences once (main measure)
#   Distinct_ECs             distinct EC numbers behind those links (breadth of known activities)

set -euo pipefail

PROJECT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT"

CHAR="${CHAR:-cazy_characterized.faa}"
FAM_MAP="${FAM_MAP:-$HOME/dbcan_db/fam-substrate-mapping.tsv}"
SUBS="${SUBS:-dbcan_1071/substrates.tsv}"
OVERVIEW="${OVERVIEW:-dbcan_1071/overview_2_substrates_group.tsv}"
OUT="${OUT:-dbcan_1071/characterized_per_group.tsv}"

for f in "$CHAR" "$FAM_MAP" "$SUBS" "$OVERVIEW"; do
    [ -s "$f" ] || { echo "missing input: $f" >&2; exit 1; }
done

# Groups for substrates that occur in fam-substrate-mapping.tsv but not among the
# 46 predicted substrates in substrates.tsv (substrate|group):
#   host glycans   spelling variant of "host glycan"
#   funoran        red-algal galactan
#   gellan         bacterial exopolysaccharide
#   lactose, gentiobiose, solabiose   disaccharides
#   galactan       grouped with beta-galactan
#   acholetin      bacterial poly-beta-1,3-GlcNAc (Macdonald et al., 2022)
# Any other substrate without a group is reported as UNMAPPED and not counted;
# add it here.
EXTRA_GROUPS='
host glycans|Glycoconjugate-degrading enzymes
funoran|Algal polysaccharidases
gellan|Microbial EPS & OM polysaccharidases
lactose|Di-/trisaccharidases
gentiobiose|Di-/trisaccharidases
solabiose|Di-/trisaccharidases
galactan|Arabinogalactanases
acholetin|Microbial EPS & OM polysaccharidases
'

awk -v extra="$EXTRA_GROUPS" '
function trim(x) { gsub(/^[ \t\r]+|[ \t\r]+$/, "", x); return x }
function family(x) { x = trim(x); sub(/_.*/, "", x); return x }

# one characterised protein: families x ECs -> substrates -> groups
function flush(   i, j, k, n, ns, s, g, key, pg, got) {
    if (hdr == "") return
    total++; got = 0
    for (i = 1; i <= nf; i++) for (j = 1; j <= ne; j++) {
        key = fams[i] SUBSEP ecs[j]
        if (!(key in pairsub)) continue
        ns = split(pairsub[key], s, /[;,]/)
        for (k = 1; k <= ns; k++) {
            s[k] = trim(s[k]); if (s[k] == "" || s[k] == "-") continue
            got = 1
            g = (s[k] in grp) ? grp[s[k]] : (s[k] in xgrp) ? xgrp[s[k]] : ""
            if (g == "") { unmapped[s[k]]++; continue }
            pg[g]
            if (!((g SUBSEP ecs[j]) in gec)) { gec[g SUBSEP ecs[j]]; nec[g]++ }
        }
    }
    if (!got) nopair++
    n = 0
    for (g in pg) {
        n++; cprot[g]++
        if (!((g SUBSEP seq) in cseq)) { cseq[g SUBSEP seq]; useq[g]++ }
    }
    if (n) withgroup++
    hdr = ""; seq = ""
}

BEGIN {
    n = split(extra, l, "\n")
    for (i = 1; i <= n; i++) if (split(l[i], p, "|") == 2) xgrp[p[1]] = p[2]
}
{ sub(/\r$/, "") }                       # tolerate Windows (CRLF) line endings
FNR == 1 { if (fi == 3) flush(); fi++ }

# 1) substrates.tsv: Count, Substrate, Group
fi == 1 { if (FNR > 1) { split($0, c, "\t"); grp[c[2]] = c[3] }; next }

# 2) fam-substrate-mapping.tsv: Substrate_high_level, Substrate_curated, Family, Name, EC_Number
fi == 2 {
    if (FNR == 1 || /^#/) next
    split($0, c, "\t")
    f = family(c[3]); e = trim(c[5]); sb = trim(c[1])
    if (f == "" || e == "" || sb == "") next
    key = f SUBSEP e
    pairsub[key] = (key in pairsub) ? pairsub[key] ";" sb : sb
    next
}

# 3) characterised proteins (FASTA)
fi == 3 {
    if (/^>/) {
        flush()
        hdr = $0; nf = 0; ne = 0; delete fseen
        n = split(substr($0, 2), h, "|")
        for (i = 2; i <= n; i++) {
            if (h[i] ~ /^[0-9]+\.[0-9]+\.[0-9]+\.([0-9]+|-)$/) ecs[++ne] = h[i]
            else if (h[i] ~ /^(GH|GT|PL|CE|AA|CBM)[0-9]/) {
                f = family(h[i]); if (!(f in fseen)) { fseen[f]; fams[++nf] = f }
            }
        }
    } else seq = seq $0
    next
}

# 4) predicted CAZymes per group
fi == 4 {
    split($0, c, "\t")
    if (FNR == 1) { for (i = 1; i in c; i++) if (c[i] == "Substrate_group") col = i; next }
    if (c[col] == "" || c[col] == "-") next
    n = split(c[col], parts, ";")
    for (i = 1; i <= n; i++) { g = trim(parts[i]); if (g != "") caz[g]++ }
}

END {
    if (fi == 3) flush()
    for (g in caz) all[g]; for (g in cprot) all[g]
    for (g in all)
        printf "%s\t%d\t%d\t%d\t%d\n", g, caz[g], cprot[g], useq[g], nec[g]
    printf "characterised proteins: %d; with >= 1 group: %d; no family/EC pair in the mapping: %d\n", \
        total, withgroup, nopair > "/dev/stderr"
    for (u in unmapped) printf "UNMAPPED substrate (add to EXTRA_GROUPS): %s (%d proteins)\n", u, unmapped[u] > "/dev/stderr"
}' "$SUBS" "$FAM_MAP" "$CHAR" "$OVERVIEW" \
  | sort -t$'\t' -k4,4nr -k1,1 \
  | { printf 'Group\tCAZymes\tCharacterized_proteins\tCharacterized_sequences\tDistinct_ECs\n'; cat; } > "$OUT.tmp"

mv "$OUT.tmp" "$OUT"
