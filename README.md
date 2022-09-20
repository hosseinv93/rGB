# rGB
soft repulsive Gay-Berne potential for LAMMPS

# Summary

The rgayberne package allows users to perform LAMMPS MD simulations with soft repuslive Gay-Berne potential.

# Dependencies

This package is only guaranteed to work for the 3Mar2020 stable version of LAMMPS. Although it may function with other versions, we do not support them and the user may need to do additional compilation steps.

# Installation instructions

compile LAMMPS using the legacy Make procedure, including `make yes-ASPHERE` to set up ASPHERE package for inclusion. This package consists the potential file for original AGy-berne potential. Then, copy 
