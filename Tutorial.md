# CAZymes Tutorial

## Requirements

- [dbCAN3](https://dbcan.readthedocs.io/en/latest/installation.html)
- [HMMER 3.4](http://hmmer.org/) (`hmmsearch`)
- [Pfam-A.hmm 38.1](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.1/)
- [SeqKit](https://bioinf.shenwei.me/seqkit/)
- [DuckDB](https://duckdb.org) 1.5.5 (Python package)
- [`mgy_proteins_pfam.parquet`](https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2026_07/mgy_proteins_pfam.parquet) Placed in the working directory


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
