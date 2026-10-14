"""
    lingradest(Bs, Rs)

Compute spatial derivatives such as grad, div, curl and curvature using the reciprocal vector technique (linear interpolation),
from field vectors `Bs` measured at positions `Rs` by N ≥ 4 spacecraft.

# Returns
A named tuple containing:
  • Rbary: Barycenter position
  • Bbc: Field at the barycenter
  • Bmag: Field magnitude at the barycenter
  • LGBx, LGBy, LGBz: Linear gradient estimators for each component
  • div: Linear divergence estimator
  • curl: Linear curl estimator
  • curv: Field-line curvature vector: 𝐛 · ∇𝐛
  • R_c: Field-line curvature radius

# References
Based on the method of Chanteur (ISSI, 1998, Ch. 11).
- [lingradest.pro](https://github.com/spedas/bleeding_edge/blob/master/projects/mms/common/curlometer/lingradest.pro)
- [lingradest.py](https://github.com/spedas/pyspedas/blob/master/pyspedas/analysis/lingradest.py#L5)
"""
function lingradest(Bs, Rs)
    Rbary = sum(Rs) / length(Rs)
    ks = reciprocal_vectors(Rs)

    # Interpolation weights 1/N + 𝒌_α ⋅ (𝐫 - Rbary) are all 1/N at the barycenter
    Bbc = sum(Bs) / length(Bs)
    Bmag = norm(Bbc)

    div = sum(map(dot, ks, Bs))
    curl = sum(map(cross, ks, Bs))
    ∇𝐛 = sum(map((k, B) -> k * B', ks, Bs))
    curv = curvature(∇𝐛, Bbc, Bmag)
    R_c = 1 / norm(curv)

    LGBx = ∇𝐛[:, 1]
    LGBy = ∇𝐛[:, 2]
    LGBz = ∇𝐛[:, 3]

    return (;
        Rbary, Bbc, Bmag,
        LGBx, LGBy, LGBz,
        div, curl, curv, R_c,
    )
end

@inline function curvature(∇𝐛, Bbc, Bmag = norm(Bbc))
    return (∇𝐛' * Bbc) / (Bmag^2)
end

# Compose any callable with a NamedTuple projection so
# downstream materialization only includes the selected columns.
# Equivalent to `args... -> NamedTuple{N}(f(args...))`
# but with the reliable inlining / dispatch behavior.
struct PickFields{Names, F}
    f::F
end
PickFields{N}(f::F) where {N, F} = PickFields{N, F}(f)

@inline (p::PickFields{N})(args...) where {N} = NamedTuple{N}(p.f(args...))

# One time sample: `n` field slices then `n` position slices, as SVectors.
struct RowKernel{n, T} end

@inline function (::RowKernel{n, T})(slices...) where {n, T}
    vs = map(SVector{3, T}, slices)
    return @inline lingradest(ntuple(i -> vs[i], Val(n)), ntuple(i -> vs[n + i], Val(n)))
end

const MatrixSet = Union{Tuple{Vararg{AbstractMatrix}}, AbstractVector{<:AbstractMatrix}}

"""
    lingradest(Bs, Rs; dim = 1, flatten = false, select = nothing)

Time-series method: `Bs` and `Rs` are tuples or vectors of matrices, each with time along `dim`. Returns a `StructArray` containing the results.

Set `flatten = true` to flatten the output arrays, making the output shape similar to the input array shape.

Pass `select = (:field1, :field2, ...)` to materialize only the listed columns.
"""
Base.@constprop :aggressive function lingradest(Bs::MatrixSet, Rs::MatrixSet; dim = 1, flatten = false, select = nothing)
    length(Bs) == length(Rs) || throw(DimensionMismatch("$(length(Bs)) fields but $(length(Rs)) positions"))
    row = RowKernel{length(Bs), Base.promote_eltype(Bs..., Rs...)}()
    f = isnothing(select) ? row : PickFields{select}(row)
    s = StructArray(mappedarray(f, map(a -> eachslice(a; dims = dim), (Bs..., Rs...))...))
    flatten || return s
    return map(StructArrays.components(s)) do c
        eltype(c) <: Number && return c
        a = reinterpret(reshape, eltype(eltype(c)), c)
        dim == 1 ? transpose(a) : a
    end
end

Base.@constprop :aggressive lingradest(B1, B2, B3, B4, R1, R2, R3, R4; kw...) =
    lingradest((B1, B2, B3, B4), (R1, R2, R3, R4); kw...)

"""
    lingradest(
        Bx1, Bx2, Bx3, Bx4,
        By1, By2, By3, By4,
        Bz1, Bz2, Bz3, Bz4,
        R1, R2, R3, R4
    )

SPEDAS-argument-compatible version of lingradest.
"""
function lingradest(
        Bx1, Bx2, Bx3, Bx4,
        By1, By2, By3, By4,
        Bz1, Bz2, Bz3, Bz4,
        R1, R2, R3, R4
    )
    Bs = (hcat(Bx1, By1, Bz1), hcat(Bx2, By2, Bz2), hcat(Bx3, By3, Bz3), hcat(Bx4, By4, Bz4))
    return lingradest(Bs, (R1, R2, R3, R4))
end

const lingradest_metadata = Dict(
    :Bmag => Dict(
        :name => "B",
        :desc => "Magnetic field magnitude at the barycenter",
    ),
    :Bbc => Dict(
        :name => "𝐁",
        :desc => "Magnetic field at the barycenter",
    ),
    :div => Dict(
        :name => "∇ ⋅ 𝐁",
        :desc => "Linear divergence estimator",
    ),
    :curl => Dict(
        :name => "∇ × 𝐁",
        :desc => "Linear curl estimator",
    ),
    :curv => Dict(
        :name => "𝐛 · ∇𝐛",
        :desc => "Field-line curvature vector: 𝐛 · ∇𝐛",
    ),
    :R_c => Dict(
        :name => "R_c",
        :desc => "Field-line curvature radius",
    ),
)
