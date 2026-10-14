# References
# - Harvey, C. C., 1998, Spatial Gradients and the Volumetric Tensor, in Analysis Methods for Multi-Spacecraft Data, edited by G. Paschmann and P. W. Daly, no. SR-001 in ISSI Scientific Reports, chap. 12, pp. 307–322, ESA Publ. Div., Noordwijk, Netherlands.
# - Paschmann, G., & Daly, P. W. (2008). Multi-spacecraft analysis methods revisited.

"Convert slowness vector ``𝐦 = 𝐧/V`` to normal vector and velocity"
function m2nV(m)
    n = m / norm(m)
    V = 1 / norm(m)
    return (; n, V)
end

"""
    ConstantVelocityApproach(positions, times)
    CVA(positions, times)

Constant Velocity Approach (CVA): boundary normal `n` and speed `V` from the crossing `times` at N ≥ 4 spacecraft
`positions` (a collection of 3-vectors), solving the timing equation ``(𝐫_α - 𝐫_1) ⋅ 𝐧 / V = t_α - t_1`` (least squares for N > 4).
"""
function ConstantVelocityApproach(positions, times)
    # Calculate time delays relative to first spacecraft
    Δts = times[2:end] .- times[1]
    # Calculate position differences relative to first spacecraft
    D = reduce(hcat, [r - positions[1] for r in positions[2:end]])'
    m = D \ Δts
    return m2nV(m)
end

"""
    ConstantVelocityApproach(positions, times, durations)

Given `durations` of the boundary crossings, calculate the thickness of the boundary
"""
function ConstantVelocityApproach(positions, times, durations)
    n, V = ConstantVelocityApproach(positions, times)
    d = V .* durations
    return (; n, V, d)
end

"""
Constant Thickness Approach (CTA) for determining boundary normal and velocity.
Based on the method described in Haaland et al., Annales Geophysicae, 2004.
"""
function ConstantThicknessApproach(positions, times, thickness)

end


"""Discontinuity Analyzer (DA) for analyzing properties of discontinuities using multi-spacecraft measurements."""
function DiscontinuityAnalyzer end