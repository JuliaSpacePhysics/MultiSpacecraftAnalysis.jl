"""
```math
𝐑 = ∑_α (𝐫_α-𝐫_b) (𝐫_α-𝐫_b)' = ∑_α 𝐫_α 𝐫_α'-𝐫_b 𝐫_b'
```

with ``𝐫_b = ∑_α 𝐫_α / N`` and `N` is the number of positions.

# References
- [paschmannMultispacecraftAnalysisMethods2008](@citet) Paschmann & Daly, 2008. Section 4.7
"""
position_tensor(rs, r̄ = mean(rs)) = sum(r -> (r - r̄) * (r - r̄)', rs)

"""
    volumetric_tensor(rs::AbstractVector{<:AbstractVector})

``\\frac{1}{N} 𝐑'``.

See also: [`position_tensor`](@ref)
"""
volumetric_tensor(rs::AbstractVector{<:AbstractVector}) = position_tensor(rs) / length(rs)

"""Calculate tetrahedron quality factors"""
function tetrahedron_quality(rs::AbstractVector{<:AbstractVector})
    Rvol = volumetric_tensor(rs)
    # Ascending for a positive semidefinite tensor; StaticArrays' eigen has no `sortby`
    F = eigen(Symmetric(Rvol))
    semiaxes = sqrt.(reverse(F.values))
    eigenvectors = F.vectors[:, end:-1:1]
    # Calculate quality parameters
    Qsr = 0.5 * (sum(semiaxes) / semiaxes[1] - 1)
    Elongation = 1 - (semiaxes[2] / semiaxes[1])
    Planarity = 1 - (semiaxes[3] / semiaxes[2])

    return (; det = det(Rvol), semiaxes, Qsr, Elongation, Planarity, eigenvectors)
end
