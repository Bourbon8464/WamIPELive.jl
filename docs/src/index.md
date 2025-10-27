# WamIPELive.jl

**WamIPELive.jl** is a Julia package for programmatic access to the **NOAA NOMADS WAM–IPE v1.2** data product.  
It provides a unified, high-level interface for retrieving, interpolating, and analysing space-weather and ionospheric parameters from the **Whole Atmosphere Model / Ionosphere Plasmasphere Electrodynamics (WAM–IPE)** coupled system.

The package is part of the *SpaceAgora* toolchain and is intended for use in upper-atmosphere and space-environment research, including:
- **Thermospheric and ionospheric modelling**
- **Satellite drag and density prediction**
- **Data-driven trajectory analysis**
- **Assimilation and model validation**

## Features

- **Direct access to NOMADS**: Automatically downloads and caches WAM–IPE NetCDF files from the NOAA public data archive.  
- **Spatio-temporal interpolation**: Evaluate any WAM–IPE variable at arbitrary latitude, longitude, altitude, and time using SciML-based interpolation.  
- **Batch queries**: Efficiently compute parameter values along trajectories or grids with a single function call.  
- **Transparent caching**: Uses a local cache directory to avoid redundant downloads and speed up repeated evaluations.  
- **Self-contained API**: Only depends on `NCDatasets.jl`, `DataInterpolations.jl`, and core Julia libraries.  

## Example

```julia
using WamIPELive
itp = WFSInterpolator(
    product="wfs", stream="ipe10",
    varname="ion_temperature", cache_dir="cache"
)

val = get_value(itp, DateTime(2025,9,30,18), -153.2, -33.4, 400.0)
println(val)
```

WamIPELive/
├── src/              → core module
├── cache/            → local NOMADS data cache
└── docs/             → documentation source
