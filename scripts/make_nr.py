#!/usr/bin/env python3
"""Deduplicate CAZyDB sequences.

Outputs:
  cazy_nr.fa         unique sequences, each named by its md5
  cazy_md5_map.tsv   md5 hash, original accession, CAZy label (one row per original sequence)
"""
import hashlib

infile = "CAZyDB.fa"
nr_path = "cazy_nr.fa"
map_path = "cazy_md5_map.tsv"

seen = set()

def process(header, seq_lines, nr_out, map_out):
    if header is None:
        return
    s = "".join(seq_lines).upper().rstrip("*")
    if not s:
        return
    h = hashlib.md5(s.encode()).hexdigest()

    # map: every original sequence gets a row
    parts = header[1:].split("|")
    acc, fams = parts[0], ",".join(p for p in parts[1:] if p)
    map_out.write(f"{h}\t{acc}\t{fams}\n")

    # nr fasta: each distinct sequence written once
    if h not in seen:
        seen.add(h)
        nr_out.write(f">{h}\n{s}\n")

with open(infile) as f, open(nr_path, "w") as nr_out, open(map_path, "w") as map_out:
    map_out.write("md5\taccession\tcazy_label\n")
    header, seq = None, []
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">"):
            process(header, seq, nr_out, map_out)
            header, seq = line, []
        else:
            seq.append(line.strip())
    process(header, seq, nr_out, map_out)

print(f"Input: {infile}")
print(f"{len(seen):,} unique sequences written to {nr_path}")
print(f"md5-to-accession map written to {map_path}")
