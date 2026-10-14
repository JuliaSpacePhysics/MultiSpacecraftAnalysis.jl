"""
    lingradest(Vs, Rs)

Compute spatial derivatives such as grad, div, curl and curvature using the reciprocal vector technique (linear interpolation),
from vectors `Vs` of any vector field 𝐕 measured at positions `Rs` by N ≥ 4 spacecraft.

# Returns
A named tuple containing:
  • rbar: Barycenter position (mean of `Rs`)
  • vbar: Field at the barycenter (mean of `Vs`)
  • vmag: Field magnitude at the barycenter
  • grad: Linear gradient estimator, `grad[i, j] = ∂_i V_j`
  • div: Linear divergence estimator
  • curl: Linear curl estimator
  • curv: Field-line curvature vector 𝐛 · ∇𝐛, with 𝐛 = 𝐕 / |𝐕|
  • rc: Field-line curvature radius

# References
Based on the method of Chanteur (ISSI, 1998, Ch. 11).
- [lingradest.pro](https://github.com/spedas/bleeding_edge/blob/master/projects/mms/common/curlometer/lingradest.pro)
- [lingradest.py](https://github.com/spedas/pyspedas/blob/master/pyspedas/analysis/lingradest.py#L5)
"""
function lingradest(Vs, Rs)
    rbar = sum(Rs) / length(Rs)
    ks = reciprocal_vectors(Rs)

    # Interpolation weights 1/N + 𝒌_α ⋅ (𝐫 - rbar) are all 1/N at the barycenter
    vbar = sum(Vs) / length(Vs)
    vmag = norm(vbar)

    grad = sum(map((k, V) -> k * V', ks, Vs))
    div = sum(map(dot, ks, Vs))
    curl = sum(map(cross, ks, Vs))
    curv = curvature(grad, vbar, vmag)
    rc = 1 / norm(curv)
    return (; rbar, vbar, vmag, grad, div, curl, curv, rc)
end

@inline function curvature(∇𝐕, 𝐕, V = norm(𝐕))
    return (∇𝐕' * 𝐕) / V^2
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
    lingradest(Vs, Rs; dim = 1, select = nothing)

Time-series method: `Vs` and `Rs` are tuples or vectors of matrices, each with time along `dim`. Returns a `StructVector` with one row per time sample.

Pass `select = (:field1, :field2, ...)` to materialize only the listed columns.
"""
Base.@constprop :aggressive function lingradest(Vs::MatrixSet, Rs::MatrixSet; dim = 1, select = nothing)
    length(Vs) == length(Rs) || throw(DimensionMismatch("$(length(Vs)) fields but $(length(Rs)) positions"))
    row = RowKernel{length(Vs), Base.promote_eltype(Vs..., Rs...)}()
    f = isnothing(select) ? row : PickFields{select}(row)
    return StructArray(mappedarray(f, map(a -> eachslice(a; dims = dim), (Vs..., Rs...))...))
end

const lingradest_metadata = Dict(
    :rbar => Dict(:name => "𝐫̄", :desc => "Barycenter position"),
    :vbar => Dict(:name => "𝐕", :desc => "Field at the barycenter"),
    :vmag => Dict(:name => "|𝐕|", :desc => "Field magnitude at the barycenter"),
    :grad => Dict(:name => "∇𝐕", :desc => "Linear gradient estimator: grad[i, j] = ∂_i V_j"),
    :div => Dict(:name => "∇ ⋅ 𝐕", :desc => "Linear divergence estimator"),
    :curl => Dict(:name => "∇ × 𝐕", :desc => "Linear curl estimator"),
    :curv => Dict(:name => "𝐛 · ∇𝐛", :desc => "Field-line curvature vector"),
    :rc => Dict(:name => "R_c", :desc => "Field-line curvature radius"),
)
