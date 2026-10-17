module MultiSpacecraftAnalysis

using LinearAlgebra
using StructArrays
using MappedArrays: mappedarray
using PrecompileTools: @setup_workload, @compile_workload
using StaticArrays
using Statistics: mean
using TimeseriesUtilities: tsync, unwrap, resolve_dim, axiskeys, rebuild

export lingradest, tlingradest
export ConstantVelocityApproach, CVA
export reciprocal_vectors
export position_tensor, volumetric_tensor, tetrahedron_quality

include("reciprocal_vector.jl")
include("tetrahedron.jl")
include("lingradest.jl")
include("timing.jl")

const CVA = ConstantVelocityApproach

@setup_workload begin
    Bs = ntuple(_ -> ones(2, 3), 4)
    Rs = ntuple(i -> repeat(Float64.((1:3)' .== i), 2), 4)
    @compile_workload lingradest(Bs, Rs)
end

end
