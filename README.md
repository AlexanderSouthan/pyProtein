# pyPolyampholyte
Does basic calculations on polyampholytes such as polypeptides/proteins/
enzymes, *e.g.* the isoelectric point similar to e.g. the "ExPASy Compute pI/Mw
tool". Methods are contained in class polyampholyte. Is intended to be included
into scripts for automated calculations.

Currently only calculations for compounds with amino acid residues (proteins/
enzymes/polypeptides) are supported: 
* Input as amino acid abundance or amino acid sequence.
* Currently contains *p*K<sub>a</sub> tables from three sources for IEP and
charge calculation.
* Calculates charge states of polyamphlytes depending on pH values.

For imformation on how to use it:
* See examples files in examples folder. 
* See docstrings.
