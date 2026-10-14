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
function reciprocal_vectors(rs)
    r̄ = sum(rs) / length(rs)
    𝐑 = position_tensor(rs, r̄)
    return map(r -> 𝐑 \ (r - r̄), rs)
end

@inline function reciprocal_vectors((r1, r2, r3, r4)::NTuple{4, Any})
    # Edges from vertex 4 share one triple product; ∑ 𝒌_α = 0 gives 𝒌_4.
    a, b, c = r1 - r4, r2 - r4, r3 - r4
    n1, n2, n3 = cross(b, c), cross(c, a), cross(a, b)
    V = dot(a, n1)
    k1, k2, k3 = n1 / V, n2 / V, n3 / V
    return (k1, k2, k3, -(k1 + k2 + k3))
end
