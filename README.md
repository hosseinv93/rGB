# Gay-Berne/WCA plugin for LAMMPS

Soft repulsive (WCA-like) Gay–Berne potential for anisotropic ellipsoids.

---

## Summary

The `gayberne/wca` pair style provides a **purely repulsive**, Weeks–Chandler–Andersen–like version of the Gay–Berne potential for anisotropic ellipsoidal particles in LAMMPS.

Conceptually, it is obtained by

* starting from the standard anisotropic Gay–Berne potential,
* **truncating it at the orientation–dependent minimum**, and
* **shifting it upward so that the minimum is at zero energy**.

The result is a **soft excluded-volume interaction** for ellipsoids, analogous to the usual LJ–WCA potential for spheres.
For spherical particles, the style reduces to a WCA-truncated Lennard–Jones.

---

## How to cite

If you use this plugin, please cite:

> Hossein Vahid, Alberto Scacchi, Maria Sammalkorpi, and Tapio Ala-Nissila,
> **“Interactions between rigid polyelectrolytes mediated by ordering and orientation of multivalent non-spherical ions in salt solutions”**,
> *Physical Review Letters* **130**, 158202 (2023).
> [https://doi.org/10.1103/PhysRevLett.130.158202](https://doi.org/10.1103/PhysRevLett.130.158202)

---

## Dependencies

* Tested with **LAMMPS 2 Aug 2023** (and nearby 2023 versions).
* Requires LAMMPS **ASPHERE** package (for ellipsoids and original Gay–Berne).

---

## Installation

1. Enable `ASPHERE` when building LAMMPS, e.g.:

   ```bash
   make yes-ASPHERE
   ```

2. Copy the plugin files into the LAMMPS `src` directory:

   ```text
   pair_gayberne_wca.cpp
   pair_gayberne_wca.h
   ```

   > Unlike older `rgayberne` plugins, this **does not** replace `pair_gayberne.*`.
   > It adds a new style called `gayberne/wca`.

3. Rebuild LAMMPS.

   Classic build:

   ```bash
   make mpi
   ```

---

## Usage

### Pair style

```lammps
pair_style gayberne/wca gamma upsilon mu cutoff
```

* `gamma`   – shift parameter controlling the “softness” of the core (same role as in `pair_style gayberne`)
* `upsilon` – exponent in the orientation–dependent energy factor η(…)
* `mu`      – exponent in the orientation–dependent energy factor χ(…)
* `cutoff`  – global neighbor-list cutoff (distance units)

**Important:**
The WCA–like truncation is handled internally at the orientation-dependent minimum of the Gay–Berne potential.
The `cutoff` here only controls how far neighbor lists are built; it does **not** reintroduce attractions.

### Pair coefficients

Per-type coefficients are set in the same way as for the standard Gay–Berne potential:

```lammps
pair_coeff i j epsilon sigma eia eib eic eja ejb ejc [cut_ij]
```

* `epsilon` – base energy scale ε_ij
* `sigma`   – base length scale σ_ij (not necessarily equal to the physical diameter)
* `eia,eib,eic` – relative well depths for side/face/end orientations of type **i**
* `eja,ejb,ejc` – relative well depths for side/face/end orientations of type **j**
* `cut_ij`  – optional pair-specific cutoff (used only for neighbor lists; the WCA truncation is still at the GB minimum)

For spherical particles, choose equal shapes and equal well-depth parameters to recover standard LJ–WCA.

---


## Notes & best practices

* Always use `atom_style ellipsoid` and set `shape` per type consistently with your intended aspect ratios, e.g.:

  ```lammps
  set type 1 shape 0.5 0.5 1.5
  ```

* Choose `cutoff` at least as large as the **largest** relevant interaction length (major axis), e.g.:

  ```lammps
  pair_style gayberne/wca 1.0 1.0 2.0 4.0
  ```

* For **purely repulsive** models, make sure you are not inadvertently adding attractive terms via a `pair_style hybrid` combination.

* For spherical particles (equal semi-axes, isotropic wells), `gayberne/wca` reduces to an **LJ–WCA** potential with a minimum at
  ( r = 2^{1/6} \sigma ) and zero beyond that.

---
