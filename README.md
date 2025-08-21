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
python3 find_protein_sequence.py [protein] [lineage]
```

For example, to get the Spike amino acid sequence for NB.1.8.1 
```bash
python3 find_protein_sequence.py S NB.1.8.1
```

or to get the Envelope amino acid sequence for BA.2
```bash
python3 find_protein_sequence.py E BA.2
```

To get the list of all proteins supported, including some whose expression is controversial,
```bash
python3 find_protein_sequence.py
```

Lineage Data Updates
====================
To update the lineage data, on a command line change into the sars-cov-2-lineage-dominant-mutations and run
```bash
git pull
```

Warranty
========
This tool comes with no warranty or other guarantee of correct operation. Use at your own risk.
