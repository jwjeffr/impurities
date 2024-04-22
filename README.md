# Vacancy Concentration in FeNiCrMnCo and FeAl

![image](https://img.shields.io/badge/python-3.10_or_greater-blue?logo=python&logoColor=white&labelColor=blue&color=grey)

[![image](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT) ![pylint](https://img.shields.io/badge/Pylint-9.61-blue?logo=python&logoColor=white)
=======

## About

This repository holds all of the code and LAMMPS inputs to compute vacancy concentration, formation energy, and formation volume in equiatomic FeNiCrMnCo and FeAl as outlined in our [paper on arXiv](https://arxiv.org/abs/2402.07324).

## Tutorial

Here, we assume the LAMMPS executable is named `lmp`.

The crucial parts of the code are the files in the `inputs/` directory and some variables that need to be specified when calling `lmp`.

Here we need a `tag` variable, which specifies the name of the subdirectory in ``inputs/``. Here, I use the name of the system of interest.

In each subdirectory in `inputs/*`, we have specified a composition, lattice, masses, and a potential. Note that these are just LAMMPS commands that are directly inputted into LAMMPS input files, so they can be any commands as long as those commands are placed in the proper order corresponding to when `inputs/${tag}/*.in` is called. The intended purpose of each file is:

-   `composition.in` defines the composition of the solution
-   `lattice.in` defines the lattice
-   `masses.in` defines the mass of each type
-   `potential.in` defines the interatomic potential

Any reference to other files in these `*.in` files must take into account that the `*.in` files are referenced in the parent directory of `inputs/`, not `inputs/` itself.

First, to run the MC-NPT equilibration, run the command `lmp -in mc.in -var tag ${tag}`, where `tag` is defined as above. Then, the MC-NPT analysis scripts need to be run:

-   `mc_md_visualization.py`
-   `order_parameter_plots.py`

This will generate a file `time.txt` which contains all recorded timesteps, sampled at a logarithmic frequency, as well as the MC-NPT visualization and order parameter plots.

Then, LAMMPS must be invoked again to perform the vacancy insertions. This is performed by the input file `insertions.in`. To run them for all timesteps, loop through all integers in `time.txt`:

``` bash
for t in $(cat time.txt)
    do
    lmp -in insertions.in -var tag ${tag} -var t ${t} -log logs/${tag}/insertions${t}.log
done
```

These runs are quite slow, especially for a pair-style without a GPU/Kokkos accelerator variant. For calculations on an HPC system, consider writing a short script to perform parallel job submissions for each timestep.

Once these runs are finished, the insertions analysis scripts need to be run:

-   `histograms.py`
-   `formation_plots.py`
-   `order_thermo.py`
-   `fluctuation.py`

This will generate text files with chemical potentials, the histograms with distributions of local formation enthalpies and volumes, plots of global formation enthalpies and volumes, plots of vacancy thermodynamics and order, and a plot of occupation fluctuations using the calculated chemical potentials.

For your own usage, modify `config.json` accordingly. Run `python validate_config.py config.json` to validate any created config file.