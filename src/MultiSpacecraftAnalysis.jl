module MultiSpacecraftAnalysis

using LinearAlgebra
using StructArrays
using MappedArrays: mappedarray
using ArraysOfArrays: flatview
using StaticArrays
using Statistics: mean
using TimeseriesUtilities: tsync, times

export lingradest, tlingradest
export ConstantVelocityApproach, CVA
export reciprocal_vector, reciprocal_vectors
export position_tensor, volumetric_tensor, tetrahedron_quality

include("reciprocal_vector.jl")
include("tetrahedron.jl")
include("lingradest.jl")
include("timing.jl")

function tlingradest end

"""
    tlingradest(fields, positions)

Interpolate and Compute spatial derivatives such as grad, div, curl and curvature using reciprocal vector technique.
"""
function tlingradest(fields, positions; kw...)
    all_data = (fields..., positions...)

    ref_times = times(first(all_data))
    have_same_timestamps = all(d -> times(d) == ref_times, all_data[2:end])

    return have_same_timestamps ?
        lingradest(all_data...; kw...) :
        lingradest(tsync(all_data...)...; kw...)
end


const CVA = ConstantVelocityApproach
const CTA = ConstantThicknessApproach
const DA = DiscontinuityAnalyzer

end
