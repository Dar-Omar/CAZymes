#!/bin/bash
# Count, for each predicted substrate, the CAZymes (rows) whose Substrate field
# lists it, and add the substrate's group (mapping below).
# Rows listing several substrates (separated by ";") count once under each.
#
# Usage (from the project directory):
#   scripts/substrate_counts.sh [overview.tsv] [out.tsv]
# Defaults:
#   overview.tsv  dbcan_1071/overview_2_substrates.tsv
#   out.tsv       dbcan_1071/substrates.tsv
#
# Substrates missing from the mapping get the group "UNMAPPED" and are reported

set -euo pipefail

OVERVIEW="${1:-dbcan_1071/overview_2_substrates.tsv}"
OUT="${2:-dbcan_1071/substrates.tsv}"

[ -s "$OVERVIEW" ] || { echo "missing input: $OVERVIEW" >&2; exit 1; }

# Substrate -> group mapping (substrate|group). Edit here to change the grouping.
SUBSTRATE_GROUPS='
cellulose|Cellulases
xylan|Xylanases/xyloglucanases
xyloglucan|Xylanases/xyloglucanases
beta-mannan|Mannanases
arabinan|Arabinogalactanases
arabinogalactan|Arabinogalactanases
acetylmannan|Mannanases
alpha-galactan|Alpha-galactosidase
beta-galactan|Arabinogalactanases
pectin|Pectinases
rhamnose|Pectinases
alpha-rhamnoside|Pectinases
beta-glucan|Beta-glucanases
chitin|Chitinases
chitosan|Chitinases
agarose|Algal polysaccharidases
carrageenan|Algal polysaccharidases
porphyran|Algal polysaccharidases
fucoidan|Algal polysaccharidases
alginate|Algal polysaccharidases
ulvan|Algal polysaccharidases
peptidoglycan|Peptidoglycanases
host glycan|Glycoconjugate-degrading enzymes
human milk oligosaccharides|Glycoconjugate-degrading enzymes
human milk polysaccharide|Glycoconjugate-degrading enzymes
glycolipid|Glycoconjugate-degrading enzymes
alpha-mannan|Glycoconjugate-degrading enzymes
galactosaminogalactan|Microbial EPS & OM polysaccharidases
beta-fucosides|Glycoconjugate-degrading enzymes
beta-glucuronan|Glycoconjugate-degrading enzymes
alpha-glucan|Alpha-glucanases
sucrose|Di-/trisaccharidases
fructan|Fructanases
trehalose|Di-/trisaccharidases
raffinose|Di-/trisaccharidases
exo-polysaccharide|Microbial EPS & OM polysaccharidases
bacterial outer membrane (OM) polysaccharides|Microbial EPS & OM polysaccharidases
xanthan|Microbial EPS & OM polysaccharidases
lignin|Ligninases
glucosylglycerate|Other
mannosylglycerate|Other
polyphenol|Other
plant glycoside|Other
glycoalkaloid|Other
cephalosporin C|Other
triglyceride|Other
'

awk -F'\t' -v OFS='\t' -v map="$SUBSTRATE_GROUPS" '
    BEGIN {
        n = split(map, lines, "\n")
        for (i = 1; i <= n; i++)
            if (split(lines[i], p, "|") == 2) group[p[1]] = p[2]
    }
    # find the Substrate column, then count rows per substrate
    FNR == 1 { for (i = 1; i <= NF; i++) if ($i == "Substrate") c = i
               if (!c) { print "no Substrate column in " FILENAME > "/dev/stderr"; exit 1 }
               next }
    { delete seen
      n = split($c, s, ";")
      for (j = 1; j <= n; j++) {
          gsub(/^ +| +$/, "", s[j])
          if (s[j] != "" && s[j] != "-" && !(s[j] in seen)) { seen[s[j]]; count[s[j]]++ }
      } }
    END {
        for (k in count) {
            g = (k in group) ? group[k] : "UNMAPPED"
            if (g == "UNMAPPED") print "WARNING: no group for substrate: " k > "/dev/stderr"
            print count[k], k, g
        }
        for (k in group) if (!(k in count)) print "note: mapped substrate not in data: " k > "/dev/stderr"
    }' "$OVERVIEW" \
  | sort -t$'\t' -k1,1nr -k2,2 | { printf 'Count\tSubstrate\tGroup\n'; cat; } > "$OUT.tmp"

mv "$OUT.tmp" "$OUT"
echo "wrote $OUT: $(( $(wc -l < "$OUT") - 1 )) substrates"
