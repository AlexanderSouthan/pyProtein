[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://app.travis-ci.com/AlexanderSouthan/pyProtein.svg?branch=master)](https://app.travis-ci.com/AlexanderSouthan/pyProtein)
[![codecov](https://codecov.io/gh/AlexanderSouthan/pyProtein/branch/master/graph/badge.svg?token=ZZ74G67EFQ)](https://codecov.io/gh/AlexanderSouthan/pyProtein)

# pyProtein
pyProtein does basic calculations on polypeptides/
proteins/enzymes, *e.g.* the isoelectric point similar to *e.g.* the "ExPASy
Compute pI/Mw tool". Methods are contained in class polyampholyte. Is intended
to be included into scripts for automated calculations.

Currently only calculations for compounds with amino acid residues (proteins/
enzymes/polypeptides) are supported: 
* Input as amino acid abundance or amino acid sequence.
* Currently contains *p*K<sub>a</sub> tables from three sources for IEP and
charge calculation.
* Calculates charge states of polyamphlytes depending on pH values.

For imformation on how to use it:
* See docstrings.
* Or look at how the tests are done in the tests folder.

# Requirements
Requirements are listed in the requirements.txt file. Please also click on the
Travis CI badge to find further details on tested conditions.

# Installation
pyProtein is part of PyPI, so install with:

```
pip install pyProtein
```
