# rGB plugin
soft repulsive Gay-Berne potential for LAMMPS

# Summary

The rgayberne package allows users to perform LAMMPS molecular dynamics (MD) simulations with soft repuslive Gay-Berne potential, which is obtained by shifting and truncating the original Gay-Berne potential.

# Citing this plugin

If you have used this plugin, please use the following citations:

Interactions between polyelectrolytes mediated by ordering and orientation of multivalent non-spherical ions in salt solutions [ 	
https://doi.org/10.48550/arXiv.2210.03492]( 	
https://doi.org/10.48550/arXiv.2210.03492).

# Dependencies

This package is only guaranteed to work for the 3Mar2020 stable version of LAMMPS. Although it may function with other versions, we do not support them and the user may need to do additional compilation steps.

# Installation instructions

First, set up ASPHERE package for inclusion using `make yes-ASPHERE` command. This package consists the potential file for original Gay-Berne potential. Then, copy-paste `pair_rgayberne.cpp` and `pair_rgayberne.h` into the LAMMPS `src` directory. Filnaly, compile LAMMPS using the legacy Make procedure `make machine-name`, e.g., `make mpi`. 

# Usage instructions
To use this potential in your simulation, use `pair_style rgayberne gamma upsilon mu cutoff`. `gamma`, `upsilon`, and `mu` are the same as the explained ones in the LAMMPS manual. `cutoff` is automatically set to the position of the original Gay-Berne potential minimum. To satisfy LAMMPS, assign a large number to `cutoff` since it doesn't affect the calculations. For example, use two times of the major diameter of the ellipsoids. 


