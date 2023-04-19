# rGB plugin
soft repulsive Gay-Berne potential for LAMMPS

# Summary

The rgayberne package allows users to perform LAMMPS molecular dynamics (MD) simulations with soft repuslive Gay-Berne potential, which is obtained by shifting and truncating the original Gay-Berne potential.

# Citing this plugin

When using this plugin, please cite:

Hossein Vahid, Alberto Scacchi, Maria Sammalkorpi, and Tapio Ala-Nissila "Interactions between rigid polyelectrolytes mediated by ordering and orientation of multivalent non-spherical ions in salt solutions", 
Physical Review Letters (2023). https://doi.org/10.1103/PhysRevLett.130.158202 . (https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.130.158202)

# Dependencies

This package is only guaranteed to work for the 3Mar2020 stable version of LAMMPS. Although it may function with other versions, we do not support them and the user may need to do additional compilation steps.

# Installation instructions

- When building LAMMPS (https://docs.lammps.org), set up the ASPHERE package for inclusion using `make yes-ASPHERE`. This package contains the script of the original Gay-Berne potential. 
- Remove `pair_gayberne.cpp` and `pair_gayberne.h` from the LAMMPS `src` directory. 
- Add `pair_rgayberne.cpp` and `pair_rgayberne.h` into the LAMMPS `src` directory. 
- Build LAMMPS. 

# Usage instructions

To employ this potential in your simulation, use `pair_style rgayberne gamma upsilon mu cutoff`, where`gamma`, `upsilon`, and `mu` are the same variables as in the original Gay-Berne potential (see LAMMPS manual). The cutoff value of the potential is automatically set to the position of the original Gay-Berne potential minimum by the plugin, however, LAMMPS requires the assignment of a value to `cutoff`. We suggest setting this value larger than the the major diameter of the particles. Note that this action does not affect the computational efficiency, computational accuracy nor the equations of motion within the simulations.

