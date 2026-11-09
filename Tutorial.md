# CAZymes Tutorial

## Requirements

- [dbCAN3](https://dbcan.readthedocs.io/en/latest/installation.html)

## Preparation

The analysis started by obtaining the reference sequences from dbCAN3 and extracting them in FASTA format using DIAMOND:

```bash
diamond getseq --db ~/dbcan_db/CAZy.dmnd --out CAZyDB.fa
```

This produces `CAZyDB.fa`, a FASTA file containing all CAZy reference sequences from the dbCAN3 database.
