"""
Compute the set of reciprocal vectors {``𝒌_α``}, which is also called the reciprocal base of the tetrahedron.

See also: [`reciprocal_vector`](@ref)
"""
@inline function reciprocal_vectors(r1, r2, r3, r4)
    # Edges from vertex 4 share one triple product; ∑ 𝒌_α = 0 gives 𝒌_4.
    a, b, c = r1 - r4, r2 - r4, r3 - r4
    n1, n2, n3 = cross(b, c), cross(c, a), cross(a, b)
    V = dot(a, n1)
    k1, k2, k3 = n1 / V, n2 / V, n3 / V
    return (k1, k2, k3, -(k1 + k2 + k3))
end

"""
    reciprocal_vector(r_βα, r_βγ, r_βλ)

Compute the reciprocal vector ``𝒌_α`` for a vertex of a tetrahedron given the relative position vectors.

```math
𝒌_α = \\frac{𝐫_{βγ} × 𝐫_{βλ}}{𝐫_{βα} ⋅ (𝐫_{βγ} × 𝐫_{βλ})}
```

where ``𝐫_{αβ} = r_β - r_α`` are relative position vectors.

# References
- Multi-spacecraft analysis methods revisited : 4.3 Properties of reciprocal vectors
"""
function reciprocal_vector(r_βα, r_βγ, r_βλ)
    numerator = cross(r_βγ, r_βλ)
    return numerator / dot(r_βα, numerator)
end

"""
    reciprocal_vector(rα, rβ, rγ, rλ)

Compute the reciprocal vector ``𝒌_α`` for a vertex of a tetrahedron given the position vectors of all vertices.

The vertices (α, β, γ, λ) must form a cyclic permutation of (1, 2, 3, 4).
"""
function reciprocal_vector(rα, rβ, rγ, rλ)
    rβα = rα - rβ
    rβγ = rγ - rβ
    rβλ = rλ - rβ
    return reciprocal_vector(rβα, rβγ, rβλ)
end


"""
    reciprocal_vector(rα, r0s)

Generalised reciprocal vector for N != 4

```math
𝐪_α = 𝐑^{-1} 𝐫_α
```
See also: [`reciprocal_vector`](@ref), [`position_tensor`](@ref)
"""
function reciprocal_vector(rα, r0s)
    r_all = (rα, r0s...)
    r̄ = mean(r_all)
    𝐑 = position_tensor(r_all, r̄)
    return 𝐑 \ (rα - r̄)
end
