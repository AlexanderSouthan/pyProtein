[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![build workflow](https://github.com/AlexanderSouthan/pyProtein/actions/workflows/python-package.yml/badge.svg)](https://github.com/AlexanderSouthan/pyProtein/actions/workflows/python-package.yml)
[![codecov](https://codecov.io/gh/AlexanderSouthan/pyProtein/branch/master/graph/badge.svg?token=ZZ74G67EFQ)](https://codecov.io/gh/AlexanderSouthan/pyProtein)

# pyProtein
pyProtein does basic calculations on polypeptides/
proteins/enzymes, *e.g.* the isoelectric point similar to *e.g.* the "ExPASy
Compute pI/Mw tool". Methods are contained in class protein. Is intended
to be included into scripts for automated calculations.

Currently only calculations for compounds with amino acid residues (proteins/
enzymes/polypeptides) are supported (so no synthetic polymers). As an input, 
the protein composition as amino acid abundance or amino acid sequence must be 
given.

The package contains functions for charge calculations of proteins:
* The calculations are done with the Henderson-Hasselbalch equation, assuming
that the individual amino acid residues do not show interactions which alter
the p*K*<sub>a</sub> values to a relevant extent.
* Calculates the charge states of proteins at a given pH value. 
* Calculates the charge curve of proteins, *i.e.* the protein charge in a given
pH value interval.
* Calculates the isoelectric point (IEP) of proteins.
* The charges can be given as total charge, separated into postive/negative
charges, or as charges carried by the individual amino acids.
* Currently contains four different p*K*<sub>a</sub> data tables for IEP and
charge calculations (Bjellqvist, IPC, IPC2, EMBOSS). It is relatively
straight-forward to implement other p*K*<sub>a</sub> tables is necessary.

The package can calculate some other protein properties:
* Nitrogen content
* Molar mass
* Mean residue molar mass
* Elemental composition

The package can handle chemical protein modifications:
* Can take modifications into account for calculated protein properties, such
as IEP, molar mass, or nitrogen content.
* Currently implemented modifications are methacryl modifications of
amino/hydroxy groups and aminoethyl modifications of carboxylic acid groups.
* Other modifications must be listed with their characteristics in the
DataFrame chain_modifications (see amino_acid_properties.py). Implemenation is
rather easy.

For imformation on how to use it:
* See docstrings.
* Or look at how the tests are done in the tests folder.

# Requirements
Requirements are listed in the requirements.txt file.

# Installation
Download repository and run:
```
pip install -e .
```

pyProtein is also part of PyPI, so install with:

```
pip install pyProtein
```
However, the newest version can always be found here.
