#!/bin/bash

# Count CAZymes per MGnify biome lineage and per biome category.

# Reads the Biomes column of the CAZyme table (made by add_biomes.py), e.g.
#   root:Host-associated:Human:Digestive system (12); root:Environmental:Aquatic:Marine (3)
# and writes one row per biome lineage:
#   CAZymes            number of CAZymes observed in this lineage
#   Biome              the MGnify biome lineage
#   Category           its category (mapping below)
#   Category_CAZymes   distinct CAZymes in the category: a CAZyme observed in several
#                      lineages of the same category counts once
# A CAZyme observed in lineages of different categories counts in each of them.
#
# Lineages not in the exact mapping are assigned by the longest matching prefix rule
# (reported, so they can be checked); lineages matching neither get "UNMAPPED".

# Defaults:
#   in.tsv   dbcan_1071/overview_2_substrates_group_biomes.tsv
#   out.tsv  dbcan_1071/unique_biomes_aggregated.tsv

set -euo pipefail

PROJECT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT"

IN="${1:-dbcan_1071/overview_2_substrates_group_biomes.tsv}"
OUT="${2:-dbcan_1071/unique_biomes_aggregated.tsv}"
[ -s "$IN" ] || { echo "missing input: $IN" >&2; exit 1; }

# Biome lineage -> category (lineage|category)
EXACT='
root|Excluded
root:Engineered|Engineered ecosystems
root:Engineered:Biogas plant|Engineered ecosystems
root:Engineered:Biogas plant:Wet fermentation|Engineered ecosystems
root:Engineered:Bioreactor|Engineered ecosystems
root:Engineered:Bioreactor:Continuous culture:Marine intertidal flat sediment inoculum:Wadden Sea-Germany|Engineered ecosystems
root:Engineered:Bioreactor:Continuous culture:Marine sediment inoculum:Wadden Sea-Germany|Engineered ecosystems
root:Engineered:Bioremediation:Terephthalate:Wastewater|Engineered ecosystems
root:Engineered:Built environment|Engineered ecosystems
root:Engineered:Food production|Engineered ecosystems
root:Engineered:Food production:Dairy products|Engineered ecosystems
root:Engineered:Food production:Fermented beverages|Engineered ecosystems
root:Engineered:Food production:Fermented vegetables|Engineered ecosystems
root:Engineered:Lab enrichment:Defined media|Engineered ecosystems
root:Engineered:Lab enrichment:Defined media:Marine media|Engineered ecosystems
root:Engineered:Lab Synthesis:Genetic cross|Engineered ecosystems
root:Engineered:Modeled:Simulated communities (microbial mixture)|Engineered ecosystems
root:Engineered:Solid waste|Engineered ecosystems
root:Engineered:Solid waste:Composting|Engineered ecosystems
root:Engineered:Solid waste:Composting:Bioreactor|Engineered ecosystems
root:Engineered:Solid waste:Landfill|Engineered ecosystems
root:Engineered:Wastewater|Engineered ecosystems
root:Engineered:Wastewater:Activated Sludge|Engineered ecosystems
root:Engineered:Wastewater:Industrial wastewater|Engineered ecosystems
root:Engineered:Wastewater:Nutrient removal:Biological phosphorus removal:Activated sludge|Engineered ecosystems
root:Engineered:Wastewater:Water and sludge|Engineered ecosystems
root:Environmental:Air|Excluded
root:Environmental:Aquatic|Excluded
root:Environmental:Aquatic:Aquaculture|Aquaculture
root:Environmental:Aquatic:Estuary|Marine
root:Environmental:Aquatic:Estuary:Sediment|Marine
root:Environmental:Aquatic:Freshwater|Freshwater
root:Environmental:Aquatic:Freshwater:Drinking water|Freshwater
root:Environmental:Aquatic:Freshwater:Drinking water:Delivery networks|Freshwater
root:Environmental:Aquatic:Freshwater:Groundwater|Freshwater
root:Environmental:Aquatic:Freshwater:Groundwater:Biofilm|Freshwater
root:Environmental:Aquatic:Freshwater:Ice:Glacier|Freshwater
root:Environmental:Aquatic:Freshwater:Lake|Freshwater
root:Environmental:Aquatic:Freshwater:Lentic|Freshwater
root:Environmental:Aquatic:Freshwater:Lentic:Sediment|Freshwater
root:Environmental:Aquatic:Freshwater:Lotic|Freshwater
root:Environmental:Aquatic:Freshwater:Lotic:Acidic|Freshwater
root:Environmental:Aquatic:Freshwater:Lotic:Low land river systems|Freshwater
root:Environmental:Aquatic:Freshwater:Lotic:Mid stream|Freshwater
root:Environmental:Aquatic:Freshwater:Lotic:Sediment|Freshwater
root:Environmental:Aquatic:Freshwater:Pond|Freshwater
root:Environmental:Aquatic:Freshwater:Pond:Sediment|Freshwater
root:Environmental:Aquatic:Freshwater:Sediment|Freshwater
root:Environmental:Aquatic:Freshwater:Storm water:Drainage pipe biofilm|Freshwater
root:Environmental:Aquatic:Lentic:Brackish|Non-marine saline & alkaline
root:Environmental:Aquatic:Marine|Marine
root:Environmental:Aquatic:Marine:Brackish|Marine
root:Environmental:Aquatic:Marine:Coastal|Marine
root:Environmental:Aquatic:Marine:Coastal:Sediment|Marine
root:Environmental:Aquatic:Marine:Hydrothermal vents:Microbial mats|Marine
root:Environmental:Aquatic:Marine:Hydrothermal vents:Sediment|Marine
root:Environmental:Aquatic:Marine:Intertidal zone|Marine
root:Environmental:Aquatic:Marine:Intertidal zone:Coral reef|Marine
root:Environmental:Aquatic:Marine:Intertidal zone:Estuary|Marine
root:Environmental:Aquatic:Marine:Intertidal zone:Oil-contaminated|Marine
root:Environmental:Aquatic:Marine:Intertidal zone:Salt marsh|Marine
root:Environmental:Aquatic:Marine:Marginal Sea|Marine
root:Environmental:Aquatic:Marine:Oceanic|Marine
root:Environmental:Aquatic:Marine:Oceanic:Aphotic zone|Marine
root:Environmental:Aquatic:Marine:Oceanic:Oil-contaminated|Marine
root:Environmental:Aquatic:Marine:Oceanic:Photic zone|Marine
root:Environmental:Aquatic:Marine:Oceanic:Sediment|Marine
root:Environmental:Aquatic:Marine:Oil-contaminated sediment|Marine
root:Environmental:Aquatic:Marine:Pelagic|Marine
root:Environmental:Aquatic:Marine:Sediment|Marine
root:Environmental:Aquatic:Non-marine Saline and Alkaline|Non-marine saline & alkaline
root:Environmental:Aquatic:Non-marine Saline and Alkaline:Hypersaline:Microbial mats|Non-marine saline & alkaline
root:Environmental:Aquatic:Non-marine Saline and Alkaline:Saline|Non-marine saline & alkaline
root:Environmental:Aquatic:Non-marine Saline and Alkaline:Saline:Microbial mats|Non-marine saline & alkaline
root:Environmental:Aquatic:Non-marine Saline and Alkaline:Salt crystallizer pond|Non-marine saline & alkaline
root:Environmental:Aquatic:Sediment|Excluded
root:Environmental:Aquatic:Thermal springs:Hot (42-90C)|Thermal springs
root:Environmental:Aquatic:Thermal springs:Near-boiling (>90C):Alkaline|Thermal springs
root:Environmental:Aquatic:Thermal springs:Sediment|Thermal springs
root:Environmental:Terrestrial:Soil|Soil
root:Environmental:Terrestrial:Soil:Agricultural|Soil
root:Environmental:Terrestrial:Soil:Contaminated|Soil
root:Environmental:Terrestrial:Soil:Forest soil|Soil
root:Environmental:Terrestrial:Soil:Oil-contaminated|Soil
root:Environmental:Terrestrial:Soil:Permafrost|Soil
root:Environmental:Terrestrial:Soil:Sand:Oil-contaminated|Soil
root:Environmental:Terrestrial:Soil:Tropical rainforest|Soil
root:Host-associated:Algae|Other host-associated
root:Host-associated:Algae:Red algae|Other host-associated
root:Host-associated:Amphibia|Other host-associated
root:Host-associated:Annelida|Other host-associated
root:Host-associated:Birds|Other host-associated
root:Host-associated:Cnidaria|Other host-associated
root:Host-associated:Fungi|Other host-associated
root:Host-associated:Human|Other host-associated
root:Host-associated:Insecta|Other host-associated
root:Host-associated:Invertebrates:Cnidaria:Coral|Other host-associated
root:Host-associated:Mammals|Other host-associated
root:Host-associated:Mammals:Respiratory system|Other host-associated
root:Host-associated:Microbial:Bacteria|Other host-associated
root:Host-associated:Mollusca|Other host-associated
root:Host-associated:Porifera|Other host-associated
root:Host-associated:Tunicates:Ascidians|Other host-associated
root:Host-associated:Animal:Digestive system:Fecal|Animal digestive system
root:Host-associated:Arthropoda:Digestive system:Gut|Animal digestive system
root:Host-associated:Arthropoda:Digestive system:Hindgut|Animal digestive system
root:Host-associated:Arthropoda:Digestive system:Oral|Animal digestive system
root:Host-associated:Birds:Digestive system|Bird digestive system
root:Host-associated:Birds:Digestive system:Ceca|Bird digestive system
root:Host-associated:Birds:Digestive system:Digestive tube|Bird digestive system
root:Host-associated:Birds:Digestive system:Digestive tube:Cecum|Bird digestive system
root:Host-associated:Birds:Digestive system:Fecal|Bird digestive system
root:Host-associated:Fish:Digestive system|Fish digestive system
root:Host-associated:Fish:Digestive system:Intestine|Fish digestive system
root:Host-associated:Human:Digestive system|Human digestive system
root:Host-associated:Human:Digestive system:Hindgut|Human digestive system
root:Host-associated:Human:Digestive system:Hindgut:Rectum|Human digestive system
root:Host-associated:Human:Digestive system:Intestine|Human digestive system
root:Host-associated:Human:Digestive system:Large intestine|Human digestive system
root:Host-associated:Human:Digestive system:Large intestine:Fecal|Human digestive system
root:Host-associated:Human:Digestive system:Large intestine:Sigmoid colon|Human digestive system
root:Host-associated:Human:Digestive system:Oral|Human digestive system
root:Host-associated:Human:Digestive system:Oral:Saliva|Human digestive system
root:Host-associated:Human:Digestive system:Oral:Subgingival plaque|Human digestive system
root:Host-associated:Human:Excretory system:Urethra:Urine|Human other body sites
root:Host-associated:Human:Reproductive system:Vagina|Human other body sites
root:Host-associated:Human:Respiratory system:Nasopharyngeal|Human other body sites
root:Host-associated:Human:Skin|Human other body sites
root:Host-associated:Human:Skin:Naris|Human other body sites
root:Host-associated:Insecta:Digestive system|Insect digestive system
root:Host-associated:Mammals:Digestive system|Mammal digestive system
root:Host-associated:Mammals:Digestive system:Fecal|Mammal digestive system
root:Host-associated:Mammals:Digestive system:Foregut:Rumen|Mammal digestive system
root:Host-associated:Mammals:Digestive system:Large intestine|Mammal digestive system
root:Host-associated:Mammals:Digestive system:Large intestine:Cecum|Mammal digestive system
root:Host-associated:Mammals:Digestive system:Large intestine:Fecal|Mammal digestive system
root:Host-associated:Mammals:Digestive system:Midgut|Mammal digestive system
root:Host-associated:Mammals:Digestive system:Stomach:Rumen|Mammal digestive system
root:Host-associated:Mammals:Gastrointestinal tract|Mammal digestive system
root:Host-associated:Mammals:Gastrointestinal tract:Intestine|Mammal digestive system
root:Host-associated:Mammals:Gastrointestinal tract:Intestine:Fecal|Mammal digestive system
root:Host-associated:Plants|Plants
root:Host-associated:Plants:Phylloplane|Plants
root:Host-associated:Plants:Phylloplane:Epiphytes|Plants
root:Host-associated:Plants:Rhizoplane:Endophytes|Plants
root:Host-associated:Plants:Rhizosphere|Plants
root:Host-associated:Plants:Rhizosphere:Soil|Plants
root:Host-associated:Plants:Root|Plants
root:Mixed|Excluded
'

# Fallback for lineages not listed above: every lineage below these prefixes
# gets the category (prefix|category); the longest matching prefix wins.
PREFIX='
root:Engineered|Engineered ecosystems
root:Environmental:Aquatic:Aquaculture|Aquaculture
root:Environmental:Aquatic:Estuary|Marine
root:Environmental:Aquatic:Freshwater|Freshwater
root:Environmental:Aquatic:Marine|Marine
root:Environmental:Aquatic:Non-marine Saline and Alkaline|Non-marine saline & alkaline
root:Environmental:Aquatic:Thermal springs|Thermal springs
root:Environmental:Terrestrial:Soil|Soil
root:Host-associated:Plants|Plants
root:Host-associated:Animal:Digestive system|Animal digestive system
root:Host-associated:Arthropoda:Digestive system|Animal digestive system
root:Host-associated:Birds:Digestive system|Bird digestive system
root:Host-associated:Fish:Digestive system|Fish digestive system
root:Host-associated:Human:Digestive system|Human digestive system
root:Host-associated:Insecta:Digestive system|Insect digestive system
root:Host-associated:Mammals:Digestive system|Mammal digestive system
root:Host-associated:Mammals:Gastrointestinal tract|Mammal digestive system
'

awk -F'\t' -v OFS='\t' -v exact="$EXACT" -v prefix="$PREFIX" '
function load(text, arr,   n, l, i, p) {
    n = split(text, l, "\n")
    for (i = 1; i <= n; i++) if (split(l[i], p, "|") == 2) arr[p[1]] = p[2]
}
function category(lin,   p, best) {
    if (lin in cat) return cat[lin]
    if (lin in ex) { cat[lin] = ex[lin]; return cat[lin] }
    best = ""
    for (p in px) if (index(lin, p ":") == 1 && length(p) > length(best)) best = p
    if (best != "") { cat[lin] = px[best]; byprefix[lin] }
    else cat[lin] = "UNMAPPED"
    return cat[lin]
}
BEGIN { load(exact, ex); load(prefix, px) }

NR == 1 {
    for (i = 1; i <= NF; i++) if ($i == "Biomes") b = i
    if (!b) { print "no Biomes column in " FILENAME > "/dev/stderr"; exit 1 }
    next
}
{
    total++
    if ($b == "" || $b == "-") { nobiome++; next }
    delete cs
    n = split($b, e, "; ")
    for (i = 1; i <= n; i++) {
        lin = e[i]
        sub(/ \([0-9]+\)$/, "", lin)                  # drop the observation count
        gsub(/^[ \t]+|[ \t]+$/, "", lin)
        if (lin == "") continue
        cnt[lin]++
        cs[category(lin)]
    }
    for (c in cs) ccnt[c]++
}
END {
    for (lin in cnt) print cnt[lin], lin, cat[lin], ccnt[cat[lin]]
    printf "CAZymes: %d; without a biome: %d\n", total, nobiome > "/dev/stderr"
    for (lin in byprefix)
        printf "note: assigned by prefix rule (check): %s -> %s (%d CAZymes)\n", lin, cat[lin], cnt[lin] > "/dev/stderr"
    for (lin in cnt) if (cat[lin] == "UNMAPPED")
        printf "WARNING: UNMAPPED lineage (add to EXACT): %s (%d CAZymes)\n", lin, cnt[lin] > "/dev/stderr"
}' "$IN" \
  | sort -t$'\t' -k2,2 \
  | { printf 'CAZymes\tBiome\tCategory\tCategory_CAZymes\n'; cat; } > "$OUT.tmp"

mv "$OUT.tmp" "$OUT"
