# sars-cov-2-lineage-protein-sequences
A tool to determine protein sequences for SARS-CoV-2 lineages

Prerequisites
=============
Python version 3 must be installed, and having Git installed is strongly recommended.

This tool is best installed by cloning from GitHub with git clone.

The dominant mutations in each lineage are synced by git cloning the repository https://github.com/alurqu/sars-cov-2-lineage-dominant-mutations with
```bash
git clone https://github.com/alurqu/sars-cov-2-lineage-dominant-mutations.git
```

Installation
============

Install this tool by running
```bash
git clone https://github.com/alurqu/sars-cov-2-lineage-protein-sequences.git
```

from the same directory from which the git clone for sars-cov-2-lineage-dominant-mutations was run


Usage
=====
To get the protein amino acid sequence for a protein from a SARS-CoV-2 lineage
```bash
python3 find_protein_sequence.py --protein [protein] --lineage [lineage] --additional-muts [additional nucleotide mutations]
```

For example, to get the Spike amino acid sequence for NB.1.8.1 
```bash
python3 find_protein_sequence.py --protein S --lineage NB.1.8.1
```

or to get the Envelope amino acid sequence for BA.2
```bash
python3 find_protein_sequence.py --protein E --lineage BA.2
```

To get the Spike amino acid sequence for XFG with an additional nucleotide mutation to C at 21567:
```bash
python3 find_protein_sequence.py --protein S --lineage XFG --additional-muts Nuc:21566C
```

To get the Spike amino acid sequence for XBB.1.5 with an additional nucleotide deletion at 21569-21571:
```bash
python3 find_protein_sequence.py --protein S --lineage XBB.1.5 --additional-muts Del:21569-21571
```

To get the Spike amino acid sequence for JN.1 with an additional nucleotide insertion of AAA just after 21565:
```bash
python3 find_protein_sequence.py --protein S --lineage JN.1 --additional-muts Ins:21565:AAA
```

To get the Spike amino acid sequence for KP.2 with an additional nucleotide mutation C23039G (S:Q493E) and a nucleotide deletion at 21653-21655 (S:S31del):
```bash
python3 find_protein_sequence.py --protein S --lineage KP.2 --additional-muts Nuc:23039G,Del:21653-21655
```

The protein can also be specified as a nucleotide range or even just a starting nucleotide. If just a starting
nucleotide is specified, translation will procede until the first stop codon.

To get the translated ORF8 amino acid sequence for EG.5
```bash
python3 find_protein_sequence.py --start 27894 --lineage EG.5
```

To get a the amino acid sequence for custom frame from 27756 to 27887 (ORF7b) for P.1:
```bash
python3 find_protein_sequence.py --start 27756 --end 27788 --lineage P.1
```

To get the list of all proteins supported, including some whose expression is controversial,
```bash
python3 find_protein_sequence.py
```

To help on additional options including reading FASTA files and scanning for start codons or Transcription Regulatory Sequences,
```bash
python3 find_protein_sequence.py --help
```

Lineage Data Updates
====================
To update the lineage data, on a command line change into the sars-cov-2-lineage-dominant-mutations directory and run
```bash
git pull
```

Usage Restriction
=================
This tool may not be used for any purpose in support of offensive biowarfare.

Warranty
========
This tool comes with no warranty or other guarantee of correct operation. Use at your own risk.
