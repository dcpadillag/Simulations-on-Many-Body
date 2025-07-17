# Simulating Lattices with Fermionic and Bosonic Degrees of Freedom

This project focuses on the simulation of quantum many-body systems on periodic lattices with both fermionic and bosonic degrees of freedom. It employs tensor network methods, particularly **Matrix Product States (MPS)** and **Matrix Product Operators (MPO)**, implemented in both **iTensor (C++)** and **TeNPy (Python)**.

## Overview

The objective is to explore strongly correlated systems defined on hexagonal lattices, with the capability to simulate:
- Pure fermionic or bosonic systems.
- Mixed systems with boson-fermion interactions.
- Time evolution using SVD-based truncation.

Simulations are built using custom MPO constructions for Hamiltonians and time-evolution operators, and utilize linear-scaling algorithms enabled by MPS representations.

## Core Features

- Construction of custom hexagonal lattice topologies.
- Fermionic and bosonic operators implemented within tensor network formalism.
- MPO representations for Hubbard-like and Bose-Hubbard models.
- Real- or imaginary-time evolution via **Trotter-Suzuki decomposition** and **SVD truncation**.
- Implementations in:
  - `iTensor` (C++): high-performance routines for low-level tensor manipulations.
  - `TeNPy` (Python): fast prototyping and advanced tensor algorithms.



