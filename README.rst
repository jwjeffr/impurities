.. _paper: https://google.com

Vacancy Concentration in FeNiCrMnCo
-----------------------------------

About
#####

This repository holds all of the code and LAMMPS inputs to compute vacancy concentration, formation energy, and formation volume in equiatomic FeNiCrMnCo as outlined in `paper`_.

Tutorial
########

Here, we assume the LAMMPS executable is named ``lmp``.

The crucial parts of the code are the files in the ``inputs/`` directory and some variables that need to be specified when calling ``lmp``.

Here we need a ``tag`` variable, which specifies the name of the subdirectory in ``inputs/``. Here, I use the name of the system of interest.

In each subdirectory in ``inputs/*``, we have specified a composition, lattice, masses, and a potential. Note that these are just LAMMPS commands that are directly inputted into LAMMPS input files, so they can be any commands as long as those commands are placed in the proper order corresponding to when ``inputs/${tag}/*.in`` is called. The intended purpose of each file is:

- ``composition.in`` defines the composition of the solution
- ``lattice.in`` defines the lattice
- ``masses.in`` defines the mass of each type
- ``potential.in`` defines the interatomic potential

Any reference to other files in these ``*.in`` files must take into account that the ``*.in`` files are referenced in the parent directory of ``inputs/``, not ``inputs/`` itself.

First, to run the MC-NPT equilibration, run the command ``lmp -in mc.in -var tag ${tag}``, where ``tag`` is defined as above. Then, the MC-NPT analysis scripts need to be run:

- ``