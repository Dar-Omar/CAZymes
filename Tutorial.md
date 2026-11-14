# CAZymes Tutorial

## Requirements

- [dbCAN3](https://dbcan.readthedocs.io/en/latest/installation.html)
- [HMMER 3.4](http://hmmer.org/) (`hmmsearch`)
- [Pfam-A.hmm 38.1](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.1/)
- [SeqKit](https://bioinf.shenwei.me/seqkit/)
- [DuckDB](https://duckdb.org) 1.5.5 (Python package)


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
srun --partition=cpu --cpus-per-task=32 --mem=180G --time=00:30:00 \
    python3 scripts/make_hits_and_ids.py
```

Outputs:
- `cazy1071_ids.txt`: sorted, unique MGnify protein IDs, the input for sequence extraction. 668,145,215 protein identifier. 

- `cache/hits_cazy<hash>_n1071.parquet`: protein–family pairs, all 1,071 families have hits.

Users are advised to adjust the `srun` options to their clusters.



