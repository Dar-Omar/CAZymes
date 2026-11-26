# CAZymes Tutorial

## Requirements

- [dbCAN3](https://dbcan.readthedocs.io/en/latest/installation.html)
- [HMMER 3.4](http://hmmer.org/) (`hmmsearch`)
- [Pfam-A.hmm 38.1](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.1/)
- [SeqKit](https://bioinf.shenwei.me/seqkit/)
- [DuckDB](https://duckdb.org) 1.5.5
- [`mgy_proteins_pfam.parquet`](https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2026_07/mgy_proteins_pfam.parquet)
- [pandas](https://pandas.pydata.org/pandas-docs/version/3.0/) 3.0.6
- [matplotlib](https://matplotlib.org/stable/index.html) 3.11.2
- [`mgy_biomes.parquet`](https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2026_07/)
- [`mgy_assemblies.tsv.gz`](https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2026_07/)


## Preparation

The analysis started by obtaining the reference sequences from dbCAN3 and extracting them in FASTA format using DIAMOND:

```bash
diamond getseq --db ~/dbcan_db/CAZy.dmnd --out CAZyDB.fa
```

This produces `CAZyDB.fa`, a FASTA file containing all CAZy reference sequences from the dbCAN3 database. This fasta file contains 4,098,879 that are not all unique sequences. To obtain the non-redundant set the [`make_nr.py`](scripts/make_nr.py) script is used:

```bash
python3 scripts/make_nr.py
```

These non-redundant sequences were searched against [Pfam-A.hmm 38.1](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.1/) using `hmmsearch` from [HMMER v3.4](http://hmmer.org/) using the shell script [`cazy_pfam.sh`](scripts/cazy_pfam.sh).


```bash
sbatch scripts/cazy_pfam.sh
```

### Select CAZyme-associated Pfam families

Count how many non-redundant CAZy sequences contain each Pfam family, then keep families found in at least 100 sequences.

```bash
grep -v '^#' pfam.domtbl \
  | awk -v OFS='\t' '{split($5,a,"."); print a[1], $4, $1}' \
  | sort -u \
  | cut -f1,2 | sort | uniq -c | sort -rn \
  | awk -v OFS='\t' '{print $2, $3, $1}' > cazyme_pfam_counts.tsv
  
awk -F'\t' '$3 >= 100 {print $1}' cazyme_pfam_counts.tsv > cazyme_pfam_min100.txt
```

### Retrieve candidate protein IDs from MGnify

Find MGnify proteins with at least hit for the selected Pfam families, using MGnify's precomputed Pfam annotations.

`mgy_proteins_pfam.parquet` from [MGnify Proteins 2026_07](https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2026_07/), should be placed in the working directory.

```bash
srun --partition=cpu --cpus-per-task=16 --mem=100G --time=00:10:00 \
    python3 scripts/make_hits_and_ids.py
```

Outputs:
- `cazy1071_ids.npy`: sorted, unique MGnify protein IDs as a NumPy int64 array. 668,145,215 protein identifier the input for sequence extraction.
- `cache/hits_cazy<hash>_n1071.parquet`: protein–family pairs, all 1,071 families have hits.

Users are advised to adjust the `srun` options to their clusters.

### Extract full-length candidate sequences

[`extract_shards_gz.sh`](scripts/extract_shards_gz.sh) downloads the 48 MGnify Proteins FASTA shards (~14 GB each) one at a time and keeps full-length sequences (`FL=1`) whose IDs are in `cazy1071_ids.npy`. [`prefetch_shards.sh`](scripts/prefetch_shards.sh) downloads next shards in parallel so the extractor. Each shard is deleted after extracting to save disk space.

Run each in its own screen session:

```bash
scripts/extract_shards_gz.sh        # session 1
scripts/prefetch_shards.sh          # session 2
```

Outputs:
- `cazy1071_fl.faa.gz`: full-length candidate sequences
- `cazy1071_fl.faa.gz.done` and `.counts`: finished shards and sequences kept per shard

Both scripts are resumable: rerun the same command after an interruption. The run is complete when the extractor exits with status 0 and `.done` lists all 48 shards. To check use the following two commands which shall return 48:

```bash
grep -c . cazy1071_fl.faa.gz.done

cut -f1 cazy1071_fl.faa.gz.counts | sort -u | wc -l
```

149,792,781 sequences to be retained after this step .


### Annotate candidates with dbCAN

Split the full-length candidates into chunks of 500,000 sequences and annotate each chunk with run_dbcan as a SLURM job array.

```bash
seqkit split2 -s 500000 -O chunks_1071 cazy1071_fl.faa.gz
mkdir -p logs
sbatch scripts/dbcan_1071.sh
```

The array size in the script must equal the number of chunks (`ls chunks_1071/*.faa.gz | wc -l`), or override it with `sbatch --array=1-<n_chunks>%20 scripts/dbcan_1071.sh`. The dbCAN database directory defaults to `$HOME/dbcan_db` (users can change it to what suits their systems). Resubmitting reruns only unfinished chunks.

Merge the per-chunk results with:

```bash
awk 'FNR == 1 && NR != 1 {next} 1' dbcan_1071/chunk_*/overview.tsv > dbcan_1071/overview.tsv.tmp \
  && mv dbcan_1071/overview.tsv.tmp dbcan_1071/overview.tsv
```
Then keep proteins supported by at least two of dbCAN's three methods:

```bash
awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "#ofTools") c=i; print; next} $c >= 2' overview.tsv > overview_2.tsv
```
Keep CAZymes with predicted substrates:

```bash
awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "Substrate") c=i; print; next} $c ~ /[^-;|, ]/' \
    overview_2.tsv > overview_2_substrates.tsv
```

- `overview.tsv`: proteins annotated by at least one method (expected 30,757,223)
- `overview_2.tsv`: CAZymes supported by at least two methods (expected 24,848,976)
- `overview_2_substrates.tsv`: CAZymes with known predicted substrates (expected 9,762,736)

```bash
awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "Substrate") c=i; next}
            {n = split($c, s, ";"); for (j=1; j<=n; j++) {gsub(/^ +| +$/, "", s[j]); if (s[j] != "" && s[j] != "-") u[s[j]]}}
            END {print length(u)}' overview_2_substrates.tsv
```
The previous command counts unique substrates (expected 46)

### Count CAZymes per substrate and group

Count, for each predicted substrate, the CAZymes whose `Substrate` field lists it, and add group the substrates. A CAZyme with several substrates (separated by ;) counts once under each.

```bash
scripts/substrate_counts.sh
```

- Output: `dbcan_1071/substrates.tsv` with columns `Count`, `Substrate`, `Group`, sorted by count

The substrate-to-group mapping (46 substrates, 16 groups plus `Other`) is defined in the script, you can edit it there to change the grouping. Substrates missing from the mapping are labelled `UNMAPPED` and reported. 

Expected: 46 substrates, no `UNMAPPED` rows.



```bash
awk -F'\t' -v OFS='\t' '
    NR == FNR { if (FNR > 1) grp[$2] = $3; next }
    FNR == 1  { for (i = 1; i <= NF; i++) if ($i == "Substrate") c = i
                print $0, "Substrate_group"; next }
    { delete seen; out = ""
      n = split($c, s, ";")
      for (j = 1; j <= n; j++) {
          gsub(/^ +| +$/, "", s[j])
          if (s[j] == "" || s[j] == "-") continue
          g = (s[j] in grp) ? grp[s[j]] : "UNMAPPED"
          if (g == "UNMAPPED") miss++
          if (!(g in seen)) { seen[g]; out = (out == "" ? g : out ";" g) }
      }
      print $0, (out == "" ? "-" : out) }
    END { if (miss) print "WARNING: " miss " substrate(s) without a group" > "/dev/stderr" }
' dbcan_1071/substrates.tsv dbcan_1071/overview_2_substrates.tsv > dbcan_1071/overview_2_substrates_group.tsv.tmp \
  && mv dbcan_1071/overview_2_substrates_group.tsv.tmp dbcan_1071/overview_2_substrates_group.tsv
```
The previous command is used to add the grouping column to dbCAN3 output

You can generate a bar plot to see the count of each group since `substrates.tsv` shows the count for individual substrates:

```bash
python3 plots/plot_groups_total.py
```

### Add MGnify biomes

Add a `Biomes` column listing every MGnify biome each CAZyme was observed in, with the number of observations, most frequent first, e.g.
`root:Host-associated:Human:Digestive system (12); root:Environmental:Aquatic:Marine (3)`.

```bash
python3 add_biomes.py
```

- Output: `dbcan_1071/overview_2_substrates_group_biomes.tsv`

A warning is printed if any CAZyme has no biome (`-`), no warnings expected.

Count CAZymes per MGnify biome and assign the biomes to one of 16 habitat categories (plus `Excluded` in this analysis). The biomes to categories mapping can be done using the script:

```bash
./biome_categories.sh
```
- Output: `dbcan_1071/unique_biomes_aggregated.tsv` with columns:
  - `CAZymes`: CAZymes observed in that biome
  - `Biome`: the biome
  - `Category`: its assigned category
  - `Category_CAZymes`: distinct CAZymes in the category (a CAZyme seen in several biomes of the same category counts once)

Then add the category (or categories) of each CAZyme to the table, using the biomes to categories assignments from the output above with:

```bash
awk -F'\t' -v OFS='\t' '
    NR == FNR { if (FNR > 1) cat[$2] = $3; next }
    FNR == 1  { for (i = 1; i <= NF; i++) if ($i == "Biomes") b = i
                if (!b) { print "no Biomes column" > "/dev/stderr"; exit 1 }
                print $0, "Category"; next }
    { if ($b == "" || $b == "-") { print $0, "-"; next }
      delete seen; out = ""
      n = split($b, e, "; ")
      for (j = 1; j <= n; j++) {
          lin = e[j]; sub(/ \([0-9]+\)$/, "", lin); gsub(/^[ \t]+|[ \t]+$/, "", lin)
          if (lin == "") continue
          c = (lin in cat) ? cat[lin] : "UNMAPPED"
          if (c == "UNMAPPED") miss++
          if (!(c in seen)) { seen[c]; out = (out == "" ? c : out ";" c) }
      }
      print $0, (out == "" ? "-" : out) }
    END { if (miss) print "WARNING: " miss " lineage(s) without a category" > "/dev/stderr" }
' dbcan_1071/unique_biomes_aggregated.tsv dbcan_1071/overview_2_substrates_group_biomes.tsv \
  > dbcan_1071/overview_2_substrates_group_biomes_category.tsv.tmp \
  && mv dbcan_1071/overview_2_substrates_group_biomes_category.tsv.tmp dbcan_1071/overview_2_substrates_group_biomes_category.tsv
```

Categories are listed in the order of the their biomes and separated by `;`.
Expected: 9,557,523 of 9,762,736 CAZymes (97.9%) in at least one habitat category; 205,213 only in excluded lineages.

### Extract characterised CAZymes

Keep the CAZy entries annotated with an EC number, used as a proxy for experimentally characterised enzymes:

```bash
awk '/^>/{keep = ($0 ~ /\|[0-9]+\.[0-9]+\.[0-9]+\.([0-9]+|-)/)} keep' CAZyDB.fa > cazy_characterized.faa
grep -c '^>' cazy_characterized.faa
```

- Output: `cazy_characterized.faa`

Headers like `>CAC83725.1|GT2|2.4.1.16` are kept; complete and partial EC numbers (e.g. `3.2.1.-`) both count.

Expected: 26,776 entries (16,166 unique sequences).


Count how many characterised CAZymes support each substrate group, as a measure of how well studied each group is.

```bash
./characterized_per_group.sh
```

It uses `$HOME/dbcan_db/fam-substrate-mapping.tsv` for dbCAN family, EC, substrate mapping.

- Output: `characterized_per_group.tsv` with columns:
  - `Group`: substrate group
  - `CAZymes`: predicted CAZymes in the group
  - `Characterized_proteins`: characterised entries linked to the group
  - `Characterized_sequences`: the same, counting identical sequences once (main measure)
  - `Distinct_ECs`: distinct EC numbers in the group (breadth of known activities)

Substrates found in the mapping but not among the predicted substrates are assigned in the `EXTRA_GROUPS` block of the script. The script reports how many entries could be linked to a group, and warns about any substrate without a group (can be added it to `EXTRA_GROUPS`).

Expected: 17,553 of 26,776 entries (65.6%) linked to at least one group; no `UNMAPPED` substrates.


### Map CAZymes to assemblies

Link each substrate-annotated CAZyme to the MGnify assemblies (ERZ) it was found in, and count how many assemblies contain each substrate group.

```bash
python3 make_protein_erz.py
```

- Outputs (in `dbcan_1071/`):
  - `protein_erz.tsv`: one row per CAZyme × assembly
  - `erz_list.txt`: distinct assemblies, one per line
  - `assemblies_per_group.tsv`: per group, the CAZymes, the number of distinct assemblies containing the group, and their percentage of all MGnify assemblies. The last row, `All groups`, counts assemblies with any substrate-annotated CAZyme.

 Percentages use the 51,274 assemblies of release 2026_07 (from the release's `STATS.md`).

Expected: 34,615,481 CAZyme–assembly pairs; 48,522 assemblies (94.6%) containing at least one substrate-annotated CAZyme.


### Fetch sampling locations from ENA

Retrieve the sample, study, country and geographic coordinates of each assembly from the [ENA Portal API](https://www.ebi.ac.uk/ena/portal/api/), 500 assemblies per request, needs internet access.

```bash
./fetch_ena_geo.sh
```

- Output: `dbcan_1071/erz_geo.tsv`, no header; columns: assembly, sample, study, country, location, checklist

The run is resumable: each batch is saved in `dbcan_1071/erz_geo_batches/`, and failed batches are reported. If any fail (the script exits with status 2), rerun the same command and only the failed batches are retried. The batch folder can be deleted once the run is complete.

Expected: 48,478 of 48,522 assemblies returned. 35,551 (73.3%) with a country and 31,861 (65.7%) with coordinates. ENA metadata change over time.
