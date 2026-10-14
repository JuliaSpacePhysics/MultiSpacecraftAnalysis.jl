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
res.curl    # Vector{SVector{3}}; also Rbary, Bbc, Bmag, LGBx/y/z, div, curv, R_c
res[1]      # NamedTuple of all fields at sample 1

lingradest(Bs, Rs; flatten = true).curl     # N×3 Matrix
lingradest(Bs, Rs; select = (:div, :curl))  # store only these
```

- A single sample: `Bs`, `Rs` as collections of 3-vectors; returns one NamedTuple.
- More than 4 spacecraft use generalized reciprocal vectors (`reciprocal_vectors(Rs)`).
- With DimensionalData loaded, `DimArray` fields return a `DimStack` over their time dimension.

## Elsewhere

- [`pyspedas.lingradest`](https://pyspedas.readthedocs.io/en/latest/analysis.html#pyspedas.lingradest) : Calculate magnetic field gradients, divergence, curl, and field line curvature from 4-point observations