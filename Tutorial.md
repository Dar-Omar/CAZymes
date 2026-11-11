# CAZymes Tutorial

## Requirements

- [dbCAN3](https://dbcan.readthedocs.io/en/latest/installation.html)
- [HMMER 3.4](http://hmmer.org/) (`hmmsearch`)
- [Pfam-A.hmm 38.1](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.1/)
- [SeqKit](https://bioinf.shenwei.me/seqkit/)

## Preparation

The analysis started by obtaining the reference sequences from dbCAN3 and extracting them in FASTA format using DIAMOND:

```bash
diamond getseq --db ~/dbcan_db/CAZy.dmnd --out CAZyDB.fa
```

This produces `CAZyDB.fa`, a FASTA file containing all CAZy reference sequences from the dbCAN3 database. This fasta file contains 4,098,879 that are not all unique sequences. To obtain the non-redundant set the [`make_nr.py`](scripts/make_nr.py) script is used:

```bash
python3 scripts/make_nr.py
```
