# MultiSpacecraftAnalysis

[![DOI](https://zenodo.org/badge/1095241504.svg)](https://doi.org/10.5281/zenodo.17686011)
[![version](https://juliahub.com/docs/General/MultiSpacecraftAnalysis/stable/version.svg)](https://juliahub.com/ui/Packages/General/MultiSpacecraftAnalysis)

Multi-spacecraft analysis methods in Julia.

## Example: spatial gradients from four or more spacecraft

```julia
using Pkg; Pkg.add("MultiSpacecraftAnalysis")
using MultiSpacecraftAnalysis

# One N×3 matrix per spacecraft: rows are time samples, columns x, y, z.
# For 3×N matrices (time along columns), pass `dim = 2`.
# Units carry through: B in nT, R in km ⇒ gradients in nT/km.
N = 100
Bs = Tuple(randn(N, 3) for _ in 1:4)                          # field at each spacecraft
Rs = Tuple(randn(N, 3) .+ 10 .* ((1:3)' .== i) for i in 1:4)  # their positions

res = lingradest(Bs, Rs)  # StructVector, one row per sample
res.curl    # Vector{SVector{3}}; also rbar, vbar, vmag, grad, div, curv, rc
res.grad    # Vector{SMatrix{3,3}}: grad[i, j] = ∂_i B_j
res[1]      # NamedTuple of all fields at sample 1

stack(res.curl; dims = 1)                   # N×3 Matrix
lingradest(Bs, Rs; select = (:div, :curl))  # store only these
```

- A single sample: `Bs`, `Rs` as collections of 3-vectors; returns one NamedTuple.
- More than 4 spacecraft use generalized reciprocal vectors (`reciprocal_vectors(Rs)`).

### Time series with different timestamps

Fields and positions are usually sampled at different times (e.g. MMS burst FGM fields and MEC positions).
`tlingradest` interpolates all of them onto the times of the first field, clipped to their common range, then runs `lingradest`:

```julia
out = tlingradest(fields, positions)  # DimArray, KeyedArray or TimeseriesUtilities.Timed inputs
out.curl                              # same container type and layout as fields[1]
tlingradest(fields, positions; select = (:div, :curl), extrapolation = true)  # other keywords go to `tsync`
```

Each output field keeps the time dimension of `fields[1]`, vectors use its component dimension, and `grad` adds `:∂` for the derivative direction.
`DimArray` outputs are collected in a `DimStack`, others in a `NamedTuple`.
See the [MMS example](https://JuliaSpacePhysics.github.io/MultiSpacecraftAnalysis.jl/dev/mms/), validated against PySPEDAS.

## Elsewhere

- [`pyspedas.lingradest`](https://pyspedas.readthedocs.io/en/latest/analysis.html#pyspedas.lingradest) : Calculate magnetic field gradients, divergence, curl, and field line curvature from 4-point observations