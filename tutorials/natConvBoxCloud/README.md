# Natural Convection with Dust Particles (natConvBoxCloud)

## Overview

This tutorial case demonstrates **natural convection with Lagrangian particle tracking** in an open-top domain with periodic side boundaries. It simulates dust/particle lift-off driven by thermal convection on a hot day.

## Physical Setup

The simulation domain is a square box (100 m × 100 m in 2D) with:

- **Bottom boundary**: Heated wall at temperature T_bottom = 313 K (simulating hot ground)
- **Top boundary**: Open to atmosphere at T_top = 310 K (zero gradient conditions)
- **Side boundaries (left/right)**: Periodic boundary conditions (infinite domain assumption)
- **Front/back**: Empty (2D simulation)

This configuration represents a vertical cross-section through the atmosphere on a hot day, where convective updrafts created by ground heating lift dust particles upward.

## Key Features

### Thermal Convection
- Buoyancy-driven flow due to temperature difference (ΔT = 3 K)
- Boussinesq approximation for density variations
- Turbulence modeling with k-ε (RANS)

### Particle Tracking
- Lagrangian particle cloud system
- Particles follow local flow field with drag forces
- Useful for tracking dust, pollutants, or other tracers
- See `constant/kinematicCloudProperties` for particle parameters

### Domain Periodicity
- Periodic side boundaries allow for representing a slice of infinite atmosphere
- Reduces boundary effects and computational cost
- Particles exiting one side re-enter on opposite side

## Numerical Parameters

| Parameter | Value | Unit |
|-----------|-------|------|
| Domain size | 100 × 100 | m |
| Grid resolution | 100 × 100 | cells |
| Temperature (bottom) | 313 | K |
| Temperature (top) | 310 | K |
| Temperature difference | 3 | K |
| Kinematic viscosity | 1.5×10⁻⁵ | m²/s |

## File Structure

```
natConvBoxCloud/
├── 0.orig/                    # Initial field definitions
│   ├── U                      # Velocity field
│   ├── T                      # Temperature field
│   ├── p_rgh                  # Pressure field
│   ├── p                      # Total pressure
│   └── alphat                 # Turbulent thermal diffusivity
├── constant/
│   ├── transportProperties    # Fluid properties
│   ├── turbulenceProperties   # Turbulence model settings
│   ├── kinematicCloudProperties  # Particle properties
│   └── g                      # Gravitational acceleration
├── system/
│   ├── blockMeshDict          # Mesh generation
│   ├── controlDict            # Solver settings
│   ├── fvSchemes              # Discretization schemes
│   └── fvSolution             # Linear solver settings
├── Allrun                     # Execution script
├── Allclean                   # Cleanup script
└── README.md                  # This file
```

## Running the Case

### Quick Start
```bash
./Allrun
```

This will:
1. Clean any previous results
2. Generate the mesh using blockMeshDict
3. Copy initial fields from 0.orig to 0
4. Run the CFD solver with particle tracking

### Manual Steps
```bash
# Clean previous run
./Allclean

# Generate mesh
blockMesh

# Initialize
cp -r 0.orig 0

# Run simulation
buoyantBoussinesqPimpleFoam
```

## Boundary Conditions

### Velocity (U)
- **Bottom**: No-slip (fixed value = 0)
- **Top**: Zero gradient (free slip, open)
- **Sides**: Cyclic (periodic)

### Temperature (T)
- **Bottom**: Fixed value = 313 K (hot surface)
- **Top**: Zero gradient (ambient temperature)
- **Sides**: Cyclic (periodic)

### Pressure (p_rgh)
- **Bottom**: Fixed flux pressure (wall constraint)
- **Top**: Fixed value = 0 Pa (atmospheric reference)
- **Sides**: Cyclic (periodic)

## Physical Interpretation

### Expected Behavior
1. **Heating phase**: Hot air near the bottom becomes less dense
2. **Convection onset**: Buoyancy forces drive warm air upward in plumes
3. **Updrafts**: Strong vertical velocities (typically 1-3 m/s)
4. **Particle lift-off**: Dust particles are entrained and transported upward
5. **Periodic recirculation**: Due to side periodicity, flow patterns repeat horizontally

### Typical Results
- Maximum vertical velocity: ~2-3 m/s
- Convection cell dimensions: 20-30 m
- Particle residence time: 100-200 s (depending on particle size)
- Dust plume height: Extends throughout domain

## Customization

### Modify Temperature Gradient
Edit `0.orig/T`:
```c
T_top  310;     // Change top temperature
dT 3;           // Change temperature difference
```

### Change Domain Size
Edit `system/blockMeshDict`:
```c
L 100;          // Domain extent (m)
N 100;          // Number of cells
```

### Adjust Particle Properties
Edit `constant/kinematicCloudProperties`:
- Particle diameter
- Density
- Number of particles
- Initial positions

### Modify Turbulence Model
Edit `constant/turbulenceProperties`:
- Switch between laminar/RANS/LES
- Change k-ε coefficients

## Solver Details

This case uses **buoyantBoussinesqPimpleFoam**:
- Incompressible flow solver with Boussinesq approximation
- Coupled energy equation for temperature
- PIMPLE algorithm (merged PISO-SIMPLE)
- Lagrangian particle tracking via kinematic cloud

## References

- OpenFOAM user guide on thermal convection
- Lagrangian particle tracking documentation
- Boussinesq approximation theory

## Notes

- This is a **2D simulation** (very thin z-dimension) suitable for teaching and initial studies
- For production work, consider extending to full 3D
- Mesh resolution (100×100) is coarse; refine for higher accuracy
- Simulation is **computational fluid dynamics (CFD)**, not actual dust transport; dust properties must be specified in `kinematicCloudProperties`
- Periodic boundaries assume infinite extent—use with caution for near-domain effects

## Troubleshooting

**Simulation diverges**: 
- Reduce time step in `system/controlDict`
- Check particle parameters in `kinematicCloudProperties`

**No particle motion**:
- Verify particles are initialized inside domain
- Check if drag models are enabled in cloud properties

**Weak convection**:
- Increase temperature difference (ΔT)
- Increase grid resolution

---

**Author**: OpenFOAM-solvers  
**Version**: v2412  
**Last Updated**: 2026
