# segfaultin
Segfaultin is a fictional protein with several exotic features that break in different tools and across formats.

The protein features:

* isopeptide bond between E34 and K11
* disulfide bond between D24S and E52S
* M1X[NLE] (norleucine is mostly identical)
* phosphorylation of S20

To do:

* Ions
* covalent ligand
* other ligand with bonds not specified

Once complete, I will do a roundtrip in Pyrosetta, PyMOL, RDKit, NGL etc. as PDB, CIF, mmTF
and see what gets lost.

For example I know PyMOL removes LINK records.
Pyrosetta makes ions disappear from CIF.
NGL wants CONECT records as opposed to LINK records in PDBs.

## Model
I keep using ubiquitin (PDB:1UBQ) for experiments.
It is smaller than GFP (PDB:1GFL), but isn't an NMR like trp-cage (PDB:1L2Y),
which takes a few lines to import in Pyrosetta.

However, the PDB oldie myoglobin (PDB:1MBN) has a haem group and is small.
I started giving it a go, but had all sorts of issues with it.

See [myoglobin](myoglobin.md).

## Making

> See [making notes](making.md)

It's actually a cool example of different movers in Pyrosetta, which may make a nice tutorial.





