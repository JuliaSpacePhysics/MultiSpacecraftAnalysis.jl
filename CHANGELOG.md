# Changelog

## [0.2.0] - 2026-10-19

### Changed

- **Breaking**: `lingradest` output fields are field-agnostic: `rbar, vbar, vmag, grad, div, curl, curv, rc`. `grad::SMatrix{3,3}` (`grad[i, j] = ∂_i V_j`) replaces `LGBx`, `LGBy`, `LGBz`; `Rbary → rbar`, `Bbc → vbar`, `Bmag → vmag`, `R_c → rc`.
- **Breaking**: `lingradest` `flatten` keyword removed; the time-series method always returns a `StructVector` (`stack(res.curl; dims = 1)` for an N×3 matrix). `DimArray` inputs to `lingradest` are no longer special-cased: `tlingradest` returns the `DimStack`, after aligning all inputs onto the first field's times.
- **Breaking**: time-series `lingradest` `dim` defaults to each input's own time dimension (e.g. `Ti` of a `DimArray`); plain matrices must pass it.
- **Breaking**: `lingradest(Vs, Rs)` is the only form; the positional 8-argument and SPEDAS 16-argument methods are removed.
- **Breaking**: `reciprocal_vector` removed; use `reciprocal_vectors(rs)[α]`.
- **Breaking**: `curv` is the field-line curvature 𝐯 · ∇𝐯 (`rc` its inverse norm). SPEDAS `curvB` keeps an extra component 𝐯 ∂_∥|𝐕| / |𝐕| along the field.
- DimensionalData output metadata no longer carries units; the estimators keep the units of the inputs.

### Removed

- `ConstantThicknessApproach` / `CTA` and `DiscontinuityAnalyzer` / `DA` stubs.
