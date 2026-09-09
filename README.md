# Gay-Berne/WCA pair style for LAMMPS

Soft, purely repulsive Gay-Berne interactions for anisotropic ellipsoids.

## Summary

The `gayberne/wca` pair style implements a Weeks-Chandler-Andersen-like
version of the Gay-Berne potential. It truncates the radial Gay-Berne
interaction at its orientation-dependent minimum and shifts the energy so
that both energy and radial force go continuously to zero there.

For an ellipsoidal pair,

```text
rho     = sigma_ij / (r - sigma_12 + gamma*sigma_ij)
r_min   = sigma_12 + (2^(1/6) - gamma)*sigma_ij

U_WCA = eta_12*chi_12 * [4*epsilon_ij*(rho^12-rho^6) + epsilon_ij]
        for r < r_min,
        0 otherwise.
```

Here `sigma_12`, `eta_12`, and `chi_12` depend on particle orientations in
the same way as in LAMMPS `pair_style gayberne`. For spherical particles
(equal shape axes and isotropic well-depth parameters), the implementation
uses the standard Lennard-Jones WCA form with cutoff
`2^(1/6)*sigma_ij`.

## How to cite

If you use this pair style, please cite:

> Hossein Vahid, Alberto Scacchi, Maria Sammalkorpi, and Tapio Ala-Nissila,
> **“Interactions between rigid polyelectrolytes mediated by ordering and
> orientation of multivalent non-spherical ions in salt solutions”**,
> *Physical Review Letters* **130**, 158202 (2023).
> [https://doi.org/10.1103/PhysRevLett.130.158202](https://doi.org/10.1103/PhysRevLett.130.158202)

The implementation is based on the standard LAMMPS Gay-Berne pair style and
also emits its Brown et al. citation through LAMMPS.

## Requirements

- LAMMPS built with the `ASPHERE` package.
- `atom_style ellipsoid`.
- All atoms of a given type must have identical shape dimensions.
- Tested with LAMMPS 10 Dec 2025. Older LAMMPS releases may require minor
  changes to neighbor-list API calls.

## Installation

Copy these files into `lammps/src/ASPHERE`:

```text
pair_gayberne_wca.cpp
pair_gayberne_wca.h
```

This adds `gayberne/wca`; it does not replace the standard `gayberne` style.

For a CMake build configured from the LAMMPS source directory:

```bash
cmake -S cmake -B build -D PKG_ASPHERE=yes
cmake --build build -j
```

For an existing configured build, the second command is sufficient. With
the traditional make build:

```bash
make yes-ASPHERE
make mpi
```

Confirm that the style is present:

```bash
build/lmp -h | grep 'gayberne/wca'
```

## Usage

### Pair style

```lammps
pair_style gayberne/wca gamma upsilon mu cutoff
```

- `gamma`: radial shift/softness parameter, as in `pair_style gayberne`.
- `upsilon`: exponent of the orientation-dependent `eta` factor.
- `mu`: exponent of the orientation-dependent `chi` factor.
- `cutoff`: global upper bound for the neighbor-list cutoff.

The attractive branch is never evaluated. A supplied cutoff larger than the
maximum possible WCA range is automatically reduced internally, avoiding an
unnecessarily large neighbor list. A cutoff that is too small still truncates
the repulsive potential early, so it must cover the largest WCA range needed
by the model.

For an anisotropic type pair, a conservative required cutoff is

```text
sqrt(2*(a_i,max^2 + a_j,max^2))
  + (2^(1/6) - gamma)*sigma_ij,
```

where `a_i,max` and `a_j,max` are the largest semi-axes. LAMMPS `set ...
shape` takes full diameters, so divide the largest specified shape value by
two when evaluating this expression. For a sphere-sphere pair, use at least
`2^(1/6)*sigma_ij`.

### Pair coefficients

Coefficients follow the standard Gay-Berne syntax:

```lammps
pair_coeff i j epsilon sigma eia eib eic eja ejb ejc [cut_ij]
```

- `epsilon`: base pair energy scale.
- `sigma`: base pair length scale; it is not necessarily a particle diameter.
- `eia eib eic`: relative well depths along the three body axes of type `i`.
- `eja ejb ejc`: corresponding values for type `j`.
- `cut_ij`: optional pair-specific upper cutoff.

For example:

```lammps
atom_style ellipsoid
set type 1 mass 1.0
set type 1 shape 3.0 1.4 1.0

pair_style gayberne/wca 1.0 3.0 1.0 4.0
pair_coeff 1 1 1.0 1.0 1.0 0.7 0.4 1.0 0.7 0.4
```

The three `shape` values above are full diameters. The internal Gay-Berne
calculation uses the corresponding semi-axes.

## Tested example and OVITO output

The accompanying input file `in.gayberne_wca_ovito` creates 125 prolate
ellipsoids, integrates 1,000 timesteps, and writes
`dump.gayberne_wca.lammpstrj` in the directory from which LAMMPS is run.

Run it with, for example:

```bash
mkdir -p /tmp/gayberne_wca_example
cd /tmp/gayberne_wca_example
/path/to/lmp -in /path/to/in.gayberne_wca_ovito
```

Load `dump.gayberne_wca.lammpstrj` in OVITO. Its atom records contain

```text
id type x y z ix iy iz quatw quati quatj quatk shapex shapey shapez
```

The `quat*` fields give the orientation quaternion and `shapex`, `shapey`,
and `shapez` are full ellipsoid diameters. Current OVITO versions normally
recognize these column names automatically. If particles initially appear as
spheres, enable ellipsoidal/aspherical particle rendering in the particle
visual element.

The relevant dump commands are:

```lammps
compute orient all property/atom quatw quati quatj quatk shapex shapey shapez
dump trajectory all custom 20 dump.gayberne_wca.lammpstrj \
     id type x y z ix iy iz c_orient[*]
dump_modify trajectory sort id \
     colname c_orient[1] quatw  colname c_orient[2] quati \
     colname c_orient[3] quatj  colname c_orient[4] quatk \
     colname c_orient[5] shapex colname c_orient[6] shapey \
     colname c_orient[7] shapez
```

## Notes

- Do not add an attractive interaction accidentally through a
  `pair_style hybrid` combination if a purely repulsive model is intended.
- `pair_modify shift` is unnecessary: the WCA energy shift is intrinsic to
  this pair style.
- The optimized CPU implementation caches orientation matrices once per atom
  and exits early outside the orientation-dependent WCA range.
- This implementation has no `single()` evaluator and currently has no
  accelerator-specific `/omp`, `/kk`, or GPU variant.
