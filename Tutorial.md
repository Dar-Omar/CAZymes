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




