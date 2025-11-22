outer(x, y) = x * y'

"""
    lingradest(B1, B2, B3, B4, R1, R2, R3, R4)

Compute spatial derivatives such as grad, div, curl and curvature using reciprocal vector technique (linear interpolation).

# Arguments
- `B1, B2, B3, B4`: 3-element vectors giving magnetic field measurements at each probe
- `R1, R2, R3, R4`: 3-element vectors giving the probe positions

# Returns
A named tuple containing:
  • Rbary: Barycenter position
  • Bbc: Magnetic field at the barycenter
  • Bmag: Magnetic field magnitude at the barycenter
  • LGBx, LGBy, LGBz: Linear gradient estimators for each component
  • LD: Linear divergence estimator
  • LCB: Linear curl estimator
  • curvature: Field-line curvature vector: 𝐛 · ∇𝐛
  • R_c: Field-line curvature radius

# References
Based on the method of Chanteur (ISSI, 1998, Ch. 11).
- [lingradest.pro](https://github.com/spedas/bleeding_edge/blob/master/projects/mms/common/curlometer/lingradest.pro)
- [lingradest.py](https://github.com/spedas/pyspedas/blob/master/pyspedas/analysis/lingradest.py#L5)
"""
function _lingradest(B1, B2, B3, B4, R1, R2, R3, R4)
    Rs = (R1, R2, R3, R4)
    Bs = (B1, B2, B3, B4)

    # Barycenter of the tetrahedron
    Rbary = (R1 .+ R2 .+ R3 .+ R4) ./ 4
    dRs = map(x -> Rbary - x, Rs)
    ks = reciprocal_vectors(R1, R2, R3, R4)

    # Magnetic field at barycenter
    # μs = @. 1 + dot(ks, dRs)
    Bbc = sum(@. (1 + dot(ks, dRs)) * Bs)
    Bmag = norm(Bbc)

    # Linear estimator
    div = sum(dot.(ks, Bs))
    curl = sum(cross.(ks, Bs))
    ∇𝐛 = sum(outer.(ks, Bs))
    # Field-line curvature
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

function _fast_lingradest(B1, B2, B3, B4, R1, R2, R3, R4)
    T = Base.promote_eltype(B1, B2, B3, B4, R1, R2, R3, R4)
    SV3 = SVector{3, T}
    return _lingradest(SV3(B1), SV3(B2), SV3(B3), SV3(B4), SV3(R1), SV3(R2), SV3(R3), SV3(R4))
end

"""
    lingradest(B1, args...)

Vectorized method for simplified usage. Returns a `StructArray` containing the results.

Set `flatten = true` to flatten the output arrays, making the output shape similar to the input array shape.
"""
Base.@constprop :aggressive function lingradest(args::AbstractMatrix...; dim = 1, flatten = false)
    new_args = map(args) do arg
        eachslice(arg; dims = dim)
    end
    s = StructArray(mappedarray(_fast_lingradest, new_args...))
    # Alternative methods (much slower)
    ## s = StructArray(Broadcast.instantiate(Broadcast.broadcasted(_fast_lingradest, new_args...)))
    ## s = StructArray(Iterators.map(_fast_lingradest, new_args...))
    return flatten ? map(StructArrays.components(s)) do c
            a = flatview(c)
            dim == 1 && ndims(a) == 2 ? a' : a
    end : s
end

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
    # Construct magnetic field vectors
    B1 = hcat(Bx1, By1, Bz1)
    B2 = hcat(Bx2, By2, Bz2)
    B3 = hcat(Bx3, By3, Bz3)
    B4 = hcat(Bx4, By4, Bz4)

    return lingradest(B1, B2, B3, B4, R1, R2, R3, R4)
end


const lingradest_metadata = Dict(
    :Bmag => Dict(
        :name => "B",
        :desc => "Magnetic field magnitude at the barycenter",
        :unit => "nT",
    ),
    :Bbc => Dict(
        :name => "𝐁",
        :desc => "Magnetic field at the barycenter",
        :unit => "nT",
    ),
    :div => Dict(
        :name => "∇ ⋅ 𝐁",
        :desc => "Linear divergence estimator",
        :unit => "nT/km",
    ),
    :curl => Dict(
        :name => "∇ × 𝐁",
        :desc => "Linear curl estimator",
        :unit => "nT/km",
    ),
    :curv => Dict(
        :name => "𝐛 · ∇𝐛",
        :desc => "Field-line curvature vector: 𝐛 · ∇𝐛",
        :unit => "nT/km",
    ),
    :R_c => Dict(
        :name => "R_c",
        :desc => "Field-line curvature radius",
        :unit => "km",
    ),
)
