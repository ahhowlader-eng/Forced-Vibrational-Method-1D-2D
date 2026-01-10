## Overview

This codebase implements a **percolation-based honeycomb lattice generator** and a **lattice-dynamics / phonon force simulation**.

The workflow is designed for **defective 2D materials (e.g. graphene-like lattices)** and supports:

* Random bond percolation
* Chiral vector–based lattice cutting
* Defect and neighbor tagging
* Periodic replication of unit cells
* Harmonic force calculations up to **4th nearest neighbors**
* Visualization of lattice geometry and displacement fields

---

## Code Structure & Purpose

### 1. Percolation Lattice Generator

**Purpose**
Generates a large honeycomb lattice with random bond percolation, applies a chiral cut, marks neighboring defect sites, and prepares the lattice for phonon simulations.

**Key Features**

* Bond occupation probability: `p = 0.95`
* Honeycomb geometry using parity-based deletions
* Chiral vector `(n,m)` definition
* Rectangular cut using `inpolygon`
* Boundary cleanup and defect expansion
* Periodic replication of lattice columns

**Main Outputs**

* `M`  : Replicated lattice mask (0 = empty, 1 = atom, 2 = neighbor)
* `mi` : Base lattice mask
* `xx` : Atomic x-coordinates
* `yy` : Atomic y-coordinates
* `LN` : Number of atoms inside chiral cut

**Saved as**
`n100_per_10000_5_exact.npz`

---

### 2. Lattice Subdomain Extraction & Visualization

**Purpose**
Extracts a sub-region of the lattice and visualizes atomic displacements and scalar fields using scatter plots.

**Key Operations**

* Windowed extraction of atoms from larger lattice
* Mapping of displacements (`UX`, `UY`, `U`) onto atomic positions
* Absolute-value filtering of directional displacements
* Flattening of 2D lattice blocks into 1D arrays for plotting

**Visualization**

* Colored scatter plot (`matplotlib.scatter`)
* Color mapped to scalar displacement magnitude
* Publication-style axes and ticks

---

### 3. Harmonic Force & Phonon Simulation

**Purpose**
Computes lattice forces and velocities using harmonic force constants up to **4th nearest neighbors**, supporting phonon DOS and transport calculations.

**Physical Model**

* 3 displacement components: `X, Y, Z`
* Neighbor shells:

  * 1st nearest
  * 2nd nearest
  * 3rd nearest
  * 4th nearest
* Direction-dependent force constants:

  * In-plane longitudinal (`Kin`)
  * In-plane transverse (`Kra`)
  * Out-of-plane (`Kout`)

**Core Arrays**

* Displacements: `U0X, U0Y, U0Z`
* Velocities: `VDX, VDY, VDZ`
* Forces: `FLX, FLY, FLZ`
* Force constant tensors per neighbor shell

**Key Characteristics**

* Exact parity rules preserved (`mod(X,2)`, `mod(Y,3)`)
* Explicit neighbor indexing (no implicit convolution)
* Time-harmonic external driving force
* Periodic boundary handling

---

## Numerical & Implementation Notes

* Loop-based implementation retained for correctness
* No vectorization applied unless explicitly requested
* Random number generation affects percolation realizations

---

## Dependencies

Python environment requires:

```text
numpy
scipy
matplotlib
```

files are loaded using `scipy.io.loadmat`.

---

## Typical Workflow

1. **Generate percolation lattice**
   → produces `M, mi, xx, yy`

2. **Extract simulation window**
   → selects active atoms for dynamics

3. **Run phonon / force simulation**
   → computes velocities and displacements

4. **Visualize or post-process**
   → DOS, displacement maps, transport analysis

---

## Intended Use Cases

* Phonon transport in defective 2D materials
* Percolation effects on vibrational modes
* Disorder-induced localization studies
* Validation against molecular dynamics or density functional phonons

---

## Reproducibility

For reproducible percolation networks, set:

```python
np.random.seed(<integer>)
```

---

## Notes for Future Extensions

* Vectorization or GPU acceleration (CuPy / JAX)
* Cluster connectivity analysis
* Phonon DOS & spectral function extraction
* Export to LAMMPS / VASP formats


## Publications

* Vibrational Properties of Disordered Carbon Nanotube (DOI: 10.13140/RG.2.2.23999.56481/2)
* Vacancy and Curvature Effects on The Phonon Properties of Single Wall Carbon Nanotube (DOI: 10.7567/JJAP.57.02CB08)
* Phonon localization in single wall carbon nanotube: Combined effect of 13C isotope and vacancies (DOI: 10.1063/5.0011810)
* Localization of the Optical Phonon Modes in Boron Nitride Nanotubes: Mixing Effect of 10 B Isotopes and Vacancies (DOI: 10.1021/acsomega.2c02792)

