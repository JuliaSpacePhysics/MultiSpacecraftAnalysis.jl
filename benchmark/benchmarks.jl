using BenchmarkTools
using MultiSpacecraftAnalysis
using StaticArrays

const SUITE = BenchmarkGroup()

const R = (
    SVector(64278.6, 17683.5, -2512.4),
    SVector(64276.4, 17704.7, -2514.06),
    SVector(64291.1, 17697.7, -2514.05),
    SVector(64282.6, 17697.4, -2525.76),
)
const B = (
    SVector(10.0, 20.0, 30.0), SVector(11.0, 19.0, 29.0),
    SVector(10.5, 21.0, 30.5), SVector(10.2, 20.3, 29.7),
)

const N = 10_000
_pos(o) = randn(N, 3) .+ reshape(o, 1, 3)
const BM = ntuple(_ -> randn(N, 3), 4)
const RM = (_pos([0.0, 0, 0]), _pos([10.0, 0, 0]), _pos([0.0, 10, 0]), _pos([0.0, 0, 10]))

SUITE["reciprocal_vector"] = @benchmarkable reciprocal_vector($(R)...)
SUITE["tetrahedron_quality"] = @benchmarkable tetrahedron_quality($(collect(R)))
SUITE["CVA"] = @benchmarkable CVA($(collect(R)), $([0.0, 1.0, 2.5, 4.0]))

SUITE["lingradest/single"] = @benchmarkable lingradest($(B)..., $(R)...)
SUITE["lingradest/matrix"] = @benchmarkable lingradest($(BM)..., $(RM)...)
