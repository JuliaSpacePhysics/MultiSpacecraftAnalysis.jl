module MultiSpacecraftAnalysis

using LinearAlgebra
using StructArrays
using MappedArrays: mappedarray
using ArraysOfArrays: flatview
using StaticArrays
using Statistics: mean

export lingradest
export ConstantVelocityApproach, CVA
export reciprocal_vector, reciprocal_vectors
export position_tensor, volumetric_tensor, tetrahedron_quality

include("reciprocal_vector.jl")
include("tetrahedron.jl")
include("lingradest.jl")
include("timing.jl")

const CVA = ConstantVelocityApproach
const CTA = ConstantThicknessApproach
const DA = DiscontinuityAnalyzer

end
