# SpaceAGORA.jl
[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://Space-FALCON-Lab.github.io/SpaceAGORA.jl/stable/)
[![Dev Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://Space-FALCON-Lab.github.io/SpaceAGORA.jl/dev/)

**Space Aerobraking Global Orbital Research and Analysis**

A high-fidelity Julia-based trajectory simulator for aerobraking mission analysis and design.

## Overview

SpaceAGORA.jl provides a comprehensive framework for simulating aerobraking maneuvers across planetary bodies. While aerobraking has been extensively validated at Mars and Venus through successful missions, this tool extends rigorous analysis capabilities to underexplored targets including Earth-return trajectories and Titan. This work addresses this gap with a systematic perturbation analysis of the aero-braking environment performed with the Julia Aerobraking Trajectory Simulator. The Julia Aerobraking Trajectory Simulator is a high-fidelity tool that couples the GRAM atmospheric suite with high-degree gravitational harmonics, third-body gravitational effects, and solar-radiation pressure. 

### Key Capabilities

- **High-Fidelity Atmosphere Modeling**: Integration with NASA's Global Reference Atmospheric Model (GRAM) suite
- **Advanced Dynamics**: High-degree gravitational harmonics, third-body perturbations, and solar radiation pressure
- **Multi-Target Support**: Extensible framework for various planetary bodies
- **Performance**: Built in Julia for computational efficiency and ease of use
using SpaceAGORA

## Installation

For detailed installation instructions, including Docker environment setup and GRAM configuration, see the [documentation](https://Space-FALCON-Lab.github.io/SpaceAGORA.jl/stable/installation/).

## Documentation

Comprehensive documentation is available at [Space-FALCON-Lab.github.io/SpaceAGORA.jl](https://Space-FALCON-Lab.github.io/SpaceAGORA.jl/stable/):

- [Installation Guide](https://Space-FALCON-Lab.github.io/SpaceAGORA.jl/stable/installation/)
- [Quick Start Tutorial](https://Space-FALCON-Lab.github.io/SpaceAGORA.jl/stable/quickstart/)
- [Docker Environment Setup](https://Space-FALCON-Lab.github.io/SpaceAGORA.jl/stable/docker/)
- [GRAM Configuration](https://Space-FALCON-Lab.github.io/SpaceAGORA.jl/stable/gram/)
- [API Reference](https://Space-FALCON-Lab.github.io/SpaceAGORA.jl/stable/api/)

## Docker Environment
See the [Docker setup guide](https://Space-FALCON-Lab.github.io/SpaceAGORA.jl/stable/docker/) for details.

## GRAM Setup

SpaceAGORA.jl requires NASA's GRAM atmospheric models:

**Space-FALCON Lab members**: Access via lab's shared drive at `ABTS/GRAM`

**External users**: Request GRAM from [NASA Software Catalog](https://software.nasa.gov/software/MFS-33888-1)

See the [GRAM setup guide](https://Space-FALCON-Lab.github.io/SpaceAGORA.jl/stable/gram/) for complete instructions.

## Features

| Category | Capabilities | Application |
|----------|--------------|-------------|
| **Atmospheric Modeling** | NASA GRAM suite integration for Mars, Venus, Earth, and Titan with time-varying conditions | High-fidelity density profiles essential for accurate aerobraking predictions |
| **Gravitational Dynamics** | Spherical harmonics up to degree 100, third-body perturbations, relativistic corrections | Captures cumulative gravitational effects over multi-orbit campaigns |
| **Solar Radiation** | Solar pressure and thermal radiation modeling | Accounts for non-conservative forces during long-duration missions |
| **Aerodynamic Forces** | Variable drag coefficients with orientation-dependent profiles | Enables spacecraft-specific configuration analysis |
| **Trajectory Analysis** | Orbital element evolution, ΔV budgets, and heat load integration | Quantifies mission performance and thermal constraints |
| **Campaign Simulation** | Multi-pass aerobraking with corridor control strategies | Complete mission design from capture to operational orbit |

## Example Applications

| Application | Target | Objective |
|-------------|--------|-----------|
| **Mars Orbit Insertion** | Mars | Propellant savings analysis for orbital capture and circularization |
| **Venus Aerocapture** | Venus | Entry corridor assessment in dense atmospheric environment |
| **Titan Exploration** | Titan | Aerobraking feasibility for outer solar system missions |
| **Earth Return Missions** | Earth | Reentry trajectory optimization for sample return architectures |
| **Mission Campaign Design** | Multi-body | End-to-end aerobraking strategy from arrival to science orbit |


## Citation

If you use SpaceAGORA.jl in your research, please cite:

```bibtex
@software{spaceagora2024,
  title = {SpaceAGORA.jl: Julia Aerobraking Trajectory Simulator},
  author = {Space-FALCON Lab},
  year = {2024},
  url = {https://github.com/Space-FALCON-Lab/SpaceAGORA.jl}
}
```

## Contributing

We welcome contributions! Please see our [contributing guidelines](https://Space-FALCON-Lab.github.io/SpaceAGORA.jl/stable/contributing/) for details.

## License

MIT License

## Acknowledgments

This work builds upon NASA's Global Reference Atmospheric Model (GRAM) suite and benefits from the Julia ecosystem's powerful numerical computing capabilities.

## Contact

- **Issues**: [GitHub Issues](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/issues)
- **Discussions**: [GitHub Discussions](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/discussions)
- **Lab**: [Space-FALCON Lab](https://www.spacefalconlab.com/)
