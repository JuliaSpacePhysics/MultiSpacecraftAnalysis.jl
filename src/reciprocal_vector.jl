"""
Compute the set of reciprocal vectors {``𝒌_α``}, which is also called the reciprocal base of the tetrahedron.

See also: [`reciprocal_vector`](@ref)
"""
@inline function reciprocal_vectors(r_α, r_β, r_γ, r_λ)
    return (
        reciprocal_vector(r_α, r_β, r_γ, r_λ),
        reciprocal_vector(r_β, r_γ, r_λ, r_α),
        reciprocal_vector(r_γ, r_λ, r_α, r_β),
        reciprocal_vector(r_λ, r_α, r_β, r_γ),
    )
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
    reciprocal_vector(rβα, rβγ, rβλ)
end


"""
    reciprocal_vector(rα, r0s::AbstractVector{<:AbstractVector})

Generalised reciprocal vector for N != 4

```math
𝐪_α = 𝐑^{-1} 𝐫_α
```
See also: [`reciprocal_vector`](@ref), [`position_tensor`](@ref)
"""
function reciprocal_vector(rα, r0s::AbstractVector{<:AbstractVector})
    r_all = [rα, r0s...]
    𝐑 = position_tensor(r_all)
    inv(𝐑) * (rα - mean(r_all))
end