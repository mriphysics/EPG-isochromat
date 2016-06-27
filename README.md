<script type="text/javascript" async
  src="//cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-MML-AM_CHTML">
</script>

# EPG-isochromat
### EPG vs isochromat-summation simulations

This repository contains Matlab code for implementing EPG and isochromat-summation simulations of FSE and SPGR pulse sequences. The methods are compared and shown to be equivalent under certain conditions. The EPG and isochromat-summation code can be used for any purpose.

This project was first presented as an ISMRM abstract at the 2016 conference and has since been written up as a full paper. For details of the ISMRM abstract version, [please see this page](docs/abstract.md)

The following document describes the full version of the code, which has been written up and submitted to *journal* and *doi link here*

Author: Shaihan Malik (twitter: [@shaihanmalik](https://twitter.com/shaihanmalik)), June 2016



## Code Description

The results in the paper are divided into 3 experiments, which are generated from three different scripts:

* **Experiment 1** concerns transient SPGR simulations with differing numbers of isochromats and is run by `SPGR_runme_1.m`
* **Experiment 2** concerns FSE simulations with different numbers of isochromats, and is run by `FSE_runme.m`
* **Experiment 3** concerns steady-state SPGR simulations with and without considering diffusion effects, and is run by `SPGR_runme_2.m`

The core functions to simulate SPGR and FSE sequences respectively are `SPGR_isochromat_sim.m` and `FSE_isochromat_sim.m` for isochromat summation and `SPGR_EPG_sim.m` and `FSE_EPG_sim.m` for Extended Phase Graphs. These have the same underlying syntax, and angles are always in radians, times in ms.


## [Detailed description of implementations](docs/detailed_description.html)

For a detailed description of the implementations please click [here](docs/detailed_description.html). Please note that the linked page uses [strapdown](http://strapdownjs.com/) for formatting and [mathjax](https://www.mathjax.org/) for rendering the maths. Please contact us if it doesn't display properly.
