# EPG-isochromat
### EPG vs isochromat-summation simulations

This repository contains Matlab code for implementing EPG and isochromat-summation simulations of FSE and SPGR pulse sequences. The methods are compared and shown to be equivalent under certain conditions in an abstract now submitted for the ISMRM conference in 2016. The EPG and isochromat-summation code can be used for any purpose. Please cite this code:[![DOI](https://zenodo.org/badge/15049/mriphysics/EPG-isochromat.svg)](https://zenodo.org/badge/latestdoi/15049/mriphysics/EPG-isochromat) if you find it useful

Author: Shaihan Malik (twitter: @shaihanmalik), November 2015

### Usage
Results for ISMRM abstract can be replicated by running:
* `SPGR_runme_1.m`: Compares transient SPGR simulations
* `FSE_runme.m`: Compares transient FSE simulations
* `SPGR_runme_2.m`: Compares steady-state SPGR simulations with and without inclusion of diffusion
