"""
    reciprocal_vectors(rs)

Reciprocal vectors {``𝒌_α``} of N ≥ 4 positions `rs`. For N = 4 this is the reciprocal base of the tetrahedron,

```math
𝒌_α = \\frac{𝐫_{βγ} × 𝐫_{βλ}}{𝐫_{βα} ⋅ (𝐫_{βγ} × 𝐫_{βλ})}, \\quad 𝐫_{αβ} = 𝐫_β - 𝐫_α
```

for (α, β, γ, λ) a cyclic permutation of (1, 2, 3, 4).
Otherwise the generalized reciprocal vectors ``𝐪_α = 𝐑^{-1} (𝐫_α - 𝐫_b)``, which coincide with ``𝒌_α`` at N = 4.

See also: [`position_tensor`](@ref)

# References
- Vogt, Paschmann & Chanteur (2008), Reciprocal vectors, in Multi-spacecraft analysis methods revisited, ISSI SR-008, Ch. 4.
"""
@inline function reciprocal_vectors(rs)
    # Offsets from one spacecraft are exact for separations ≪ |𝐫|.
    r₀ = first(rs)
    ds = map(r -> r - r₀, rs)
    d̄ = sum(ds) / length(ds)
    ds = map(d -> d - d̄, ds)
    X = reduce(hcat, ds)
    𝐑 = X * X'
    # 𝐑⁻¹ = adj(𝐑) / det(𝐑), cheaper than `inv` or `\`: for symmetric 𝐑 the adjugate's columns are cross products of 𝐑's columns.
    a, b, c = 𝐑[:, 1], 𝐑[:, 2], 𝐑[:, 3]
    A = hcat(cross(b, c), cross(c, a), cross(a, b))
    s = inv(dot(a, A[:, 1]))
    return map(d -> s * (A * d), ds)
end
