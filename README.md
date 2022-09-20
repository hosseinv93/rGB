# rGB
soft repulsive Gay-Berne potential for LAMMPS

# Summary

The rgayberne package allows users to perform LAMMPS MD simulations with soft repuslive Gay-Berne potential, which is obtained by shifting and truncating the original Gay-Berne potential. 

# Dependencies

This package is only guaranteed to work for the 3Mar2020 stable version of LAMMPS. Although it may function with other versions, we do not support them and the user may need to do additional compilation steps.

# Installation instructions

First, set up ASPHERE package for inclusion using `make yes-ASPHERE` command. This package consists the potential file for original Gay-berne potential. Then, copy-paste `pair_rgayberne.cpp` and `pair_rgayberne.h` into the LAMMPS `src` directory. Filnaly, compile LAMMPS using the legacy Make procedure `make machine-name`, e.g., `make mpi`. 

# Usage instructions
To use this potential in your simulation, use `pair_style rgayberne gamma upsilon mu cutoff`. `gamma`, `upsilon`, `mu` are re the same as the explained one in LAMMPS manual. `cutoff` is automatically set to the position of the original Gay-berne potential and, here, is just definded to satisfy LAMMPS by assigining a large number to it. For example, two times of the largest diameter of the ellipsoids. 
