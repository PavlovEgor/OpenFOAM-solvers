# OpenFOAM-solvers

Custom OpenFOAM v2412 Solvers and Libraries

## Description

OpenFOAM-solvers is a collection of specialized solvers and libraries for OpenFOAM v2412. This repository contains custom implementations and extensions of CFD solvers.


## Repository Structure

```
OpenFOAM-solvers/
├── solvers/          # Custom OpenFOAM solvers
├── src/              # Additional libraries and classes
├── tutorials/        # Test cases and tutorials
└── README.md         # This file
```

### Solvers
Pre-built executable solvers for various applications.

### Src
Custom C++ libraries and classes that extend OpenFOAM functionality:

### Tutorials
Complete tutorial cases with:
- Pre-configured meshes
- Input parameter files
- Expected results for validation

## Available Solvers

| Solver Name | Description | Features | Status |
|-------------|-------------|----------|--------|
| **buoyantBoussinesqKinematicPimpleFoam** | Buoyant flow with kinematic particle tracking | Heat transfer, turbulence modeling, Lagrangian particles, dynamic mesh, radiation | ✓ Stable |

### Solver Details

#### buoyantBoussinesqKinematicPimpleFoam
Transient solver for buoyant, turbulent flow of incompressible fluids with kinematic Lagrangian particle tracking. Combines buoyant heat transfer (Boussinesq approximation) with kinematic particle cloud evolution for multiphase flow simulations.

**Key Features:**
- Boussinesq approximation for buoyancy-driven flows
- Kinematic particle cloud tracking
- PIMPLE algorithm for pressure-velocity coupling
- Optional dynamic mesh capabilities
- Mesh topology changes support
- Radiation heat transfer modeling
- Turbulence modeling (incompressible)


## Installation

### Prerequisites
- OpenFOAM v2412

### Building

Default OpenFOAM solver building.

## License

This work is based on OpenFOAM, which is distributed under the GNU General Public License v3.

Copyright (C) 2026 Egor Pavlov

---

**Author:** Egor Pavlov  
**OpenFOAM Version:** 2412  
**Last Updated:** July 2026
