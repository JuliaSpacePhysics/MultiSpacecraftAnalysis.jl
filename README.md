# MultiSpacecraftAnalysis

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSpacePhysics.github.io/MultiSpacecraftAnalysis.jl/dev/)
[![DOI](https://zenodo.org/badge/1095241504.svg)](https://doi.org/10.5281/zenodo.17686011)
[![version](https://juliahub.com/docs/General/MultiSpacecraftAnalysis/stable/version.svg)](https://juliahub.com/ui/Packages/General/MultiSpacecraftAnalysis)

[![Build Status](https://github.com/JuliaSpacePhysics/MultiSpacecraftAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSpacePhysics/MultiSpacecraftAnalysis.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSpacePhysics/MultiSpacecraftAnalysis.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSpacePhysics/MultiSpacecraftAnalysis.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Multi-spacecraft analysis methods in Julia.

## Quickstart

Linear gradient estimation (`lingradest`) from four-spacecraft magnetic-field samples — single point:

```julia
using Pkg; Pkg.add("MultiSpacecraftAnalysis")
using MultiSpacecraftAnalysis

R1 = [64278.6, 17683.5, -2512.4]   # km
R2 = [64276.4, 17704.7, -2514.06]
R3 = [64291.1, 17697.7, -2514.05]
R4 = [64282.6, 17697.4, -2525.76]

B1 = [10.0, 20.0, 30.0]            # nT
B2 = [11.0, 19.0, 29.0]
B3 = [10.5, 21.0, 30.5]
B4 = [10.2, 20.3, 29.7]

out = lingradest(B1, B2, B3, B4, R1, R2, R3, R4)
out.div     # ∇·B
out.curl    # ∇×B
out.curv    # field-line curvature 𝐛·∇𝐛
out.R_c     # curvature radius
```

Vectorized over a timeseries (N×3 matrices, time along dim 1):

```julia
res = lingradest(B1m, B2m, B3m, B4m, R1m, R2m, R3m, R4m)
```

With `DimensionalData` time-stamped inputs (auto-synchronizes timestamps via `tsync` if they differ):

```julia
using DimensionalData
out = tlingradest(fields, positions)   # fields, positions :: NTuple{4, AbstractDimMatrix}
out.Bbc                                # DimArray of barycenter field
```

Boundary-normal timing analysis (CVA):

```julia
n, V = CVA(positions, crossing_times)             # normal + speed
n, V, d = CVA(positions, crossing_times, durations)   # + boundary thickness
```

Tetrahedron quality:

```julia
tetrahedron_quality([R1, R2, R3, R4])    # (; det, semiaxes, Qsr, Elongation, Planarity, eigenvectors)
```

## Elsewhere

- [`pyspedas.lingradest`](https://pyspedas.readthedocs.io/en/latest/analysis.html#pyspedas.lingradest) : Calculate magnetic field gradients, divergence, curl, and field line curvature from 4-point observations