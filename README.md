# Coalescent Tree

Tools for handling coalescent trees obtained with SLiM or msprime.

## Usage

All the functions take a .trees file as an input (generated with SLiM or msprime for instance).

The module was created to handle haploid coalescent tree optained with SLiM. In the folder 'SLiM_templates' there are two examples of SLiM code that generate such coalescent trees.

## Installation

`python setup.py install`

The following python modules have to be installed:
- numpy
- random
- matplotlib
- msprime
- tskit
- cairosvg

## Potential issues

### The function handle_haploid always return the Exception: Tree not coalesced !

The simulations were run using SLiM 3.6 build directly from the source code on git (https://github.com/MesserLab/SLiM/).
One can install SLiM following the installation guide from the web site (https://messerlab.org/slim/) but there may be some issues on the behaviour of the SLiM function 'addRecombinant' (it may not create NULL genomes as it is supposed to).
- Check that the version of SLiM is SLiM 3.6.
- If there is still the issue then build SLiM directly from the source code (https://github.com/MesserLab/SLiM/)
