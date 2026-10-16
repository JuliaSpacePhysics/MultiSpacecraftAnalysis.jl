# Changelog

## [0.2.0] - Unreleased

### Changed

- **Breaking**: `lingradest` output fields are field-agnostic: `rbar, vbar, vmag, grad, div, curl, curv, rc`. `grad::SMatrix{3,3}` (`grad[i, j] = ∂_i V_j`) replaces `LGBx`, `LGBy`, `LGBz`; `Rbary → rbar`, `Bbc → vbar`, `Bmag → vmag`, `R_c → rc`.
- **Breaking**: `lingradest` `flatten` keyword removed; the time-series method always returns a `StructVector` (`stack(res.curl; dims = 1)` for an N×3 matrix). `DimArray` inputs to `lingradest` are no longer special-cased: `tlingradest` returns the `DimStack`, after aligning all inputs onto the first field's times. Its vector and tensor fields are numeric arrays laid out like the first field and reusing its component dimension (`grad` adds `Dim{:∂}` for the derivative direction).
- **Breaking**: `lingradest(Vs, Rs)` is the only form; the positional 8-argument and SPEDAS 16-argument methods are removed.
- **Breaking**: `reciprocal_vector` removed; use `reciprocal_vectors(rs)[α]`.

### Removed

- `ConstantThicknessApproach` / `CTA` and `DiscontinuityAnalyzer` / `DA` stubs.

### Changed

- DimensionalData output metadata no longer carries units; the estimators keep the units of the inputs.
