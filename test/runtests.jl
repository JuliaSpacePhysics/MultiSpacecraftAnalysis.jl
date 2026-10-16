using MultiSpacecraftAnalysis
using Test
using LinearAlgebra
using StaticArrays
using DimensionalData
using Aqua
using PythonCall

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(MultiSpacecraftAnalysis)
end

# skip nan
_is_approx(x, y; kw...) = isnan(x) && isnan(y) || isapprox(x, y; kw...)
function _is_approx(x::AbstractArray, y::AbstractArray; kw...)
    for i in eachindex(x, y)
        _is_approx(x[i], y[i]; kw...) || return false
    end
    return true
end

@testset "Cross validation with PySPEDAS.jl" begin
    # Reference: https://github.com/spedas/pyspedas/blob/master/pyspedas/projects/mms/tests/test_mms_curlometer.py#L4
    # https://github.com/spedas/pyspedas/blob/283bc4ce21a302af3b894451b489a6175dd8fd5d/pyspedas/projects/mms/fgm_tools/mms_lingradest.py#L8
    using PySPEDAS

    @py import pyspedas.projects.mms.tests.test_mms_curlometer: CurlTestCases
    @py import pyspedas.projects.mms.fgm_tools.mms_lingradest: mms_lingradest
    pytest = CurlTestCases()
    pytest.test_lingradest()


    field_names = ("mms1_fgm_b_gse_brst_l2", "mms2_fgm_b_gse_brst_l2_i", "mms3_fgm_b_gse_brst_l2_i", "mms4_fgm_b_gse_brst_l2_i")
    position_names = ("mms1_fgm_r_gse_brst_l2_i", "mms2_fgm_r_gse_brst_l2_i", "mms3_fgm_r_gse_brst_l2_i", "mms4_fgm_r_gse_brst_l2_i")

    func(x) = DimArray(get_data(x))[:, 1:3]
    fields = func.(field_names)
    positions = func.(position_names)

    res = tlingradest(fields, positions)

    @test _is_approx(res.vbar, get_data("Bbc_lingradest")[:, 2:4])
    @test _is_approx(res.vmag, get_data("Bt_lingradest"))
    scale = 1000
    @test _is_approx(res.grad[:, :, 1] .* scale, get_data("gradBx_lingradest"))
    @test _is_approx(res.grad[:, :, 2] .* scale, get_data("gradBy_lingradest"))
    @test _is_approx(res.grad[:, :, 3] .* scale, get_data("gradBz_lingradest"))
    @test _is_approx(res.div .* scale, get_data("divB_nT/1000km_lingradest"))
    @test _is_approx(res.curl .* scale, get_data("curlB_nT/1000km_lingradest")[:, 2:4])
    @test _is_approx(res.curv .* scale, get_data("curvB_lingradest"))
    @test _is_approx(res.rc ./ scale, get_data("Rc_1000km_lingradest"))
    @test _is_approx(res.curl .* scale .* 0.8, get_data("jtotal_lingradest"))
end


@testset "lingradest is exact for linear fields" begin
    G = @SMatrix randn(3, 3)  # G[i, j] = ∂_i B_j
    B0 = @SVector randn(3)
    field(r) = B0 + G' * r
    curlG = SVector(G[2, 3] - G[3, 2], G[3, 1] - G[1, 3], G[1, 2] - G[2, 1])
    for n in (4, 6)
        Rs = ntuple(_ -> 10 * @SVector(randn(3)), n)
        res = lingradest(map(field, Rs), Rs)
        @test res.grad ≈ G
        @test res.div ≈ tr(G)
        @test res.curl ≈ curlG
        @test res.vbar ≈ field(res.rbar)
    end

    N, n = 8, 5
    Rs = [10 .* randn(N, 3) for _ in 1:n]
    Bs = [stack(field, eachrow(R); dims = 1) for R in Rs]
    res = lingradest(Bs, Rs)
    @test all(≈(curlG), res.curl)
    @test all(≈(tr(G)), res.div)
end

@testset "tlingradest aligns time series" begin
    using Dates
    using TimeseriesUtilities: Timed
    G = @SMatrix randn(3, 3)
    B0 = @SVector randn(3)
    R0s = ntuple(i -> 10 * @SVector(randn(3)), 4)
    v = @SVector randn(3)
    # Linear motion through a linear field: linear interpolation is exact.
    R(α, s) = R0s[α] + v * s
    B(α, s) = B0 + G' * R(α, s)
    sec(t) = Dates.value(t - DateTime(2020)) / 1000
    t(offset, dt, n) = DateTime(2020) .+ Millisecond.(offset .+ dt .* (0:(n - 1)))

    comp = Dim{:comp}([:x, :y, :z])
    tVs = ntuple(α -> t(100α, 1000, 60), 4)
    Vs = ntuple(α -> DimArray(stack(s -> B(α, sec(s)), tVs[α]), (comp, Ti(tVs[α]))), 4)
    tR = t(0, 10_000, 7)
    Rs = ntuple(α -> Timed(stack(s -> R(α, sec(s)), tR; dims = 1), tR; dim = 1), 4)

    res = tlingradest(Vs, Rs)
    tout = filter(s -> s >= tVs[4][1] && s <= tR[end], tVs[1])
    @test dims(res.curl) == (comp, Ti(tout))
    @test dims(res.grad) == (Dim{:∂}([:x, :y, :z]), comp, Ti(tout))
    @test all(g -> g ≈ G, eachslice(parent(res.grad); dims = 3))
    @test parent(res.rbar) ≈ stack(s -> sum(α -> R(α, sec(s)), 1:4) / 4, tout)

    Vts = ntuple(α -> Timed(permutedims(parent(Vs[α])), tVs[α]; dim = 1), 4)
    div, ts = tlingradest(Vts, Rs; select = (:div,)).div
    @test ts == tout
    @test div ≈ fill(tr(G), length(tout))
end

@testset "lingradest select kwarg" begin
    N = 64
    BM = ntuple(_ -> randn(N, 3), 4)
    RM = (
        randn(N, 3) .+ reshape([0.0, 0, 0], 1, 3),
        randn(N, 3) .+ reshape([10.0, 0, 0], 1, 3),
        randn(N, 3) .+ reshape([0.0, 10, 0], 1, 3),
        randn(N, 3) .+ reshape([0.0, 0, 10], 1, 3),
    )
    full = lingradest(BM, RM)
    sel = lingradest(BM, RM; select = (:vmag, :div))
    @test propertynames(sel) == (:vmag, :div)
    @test sel.vmag ≈ full.vmag
    @test sel.div ≈ full.div
end

@testset "CVA recovers a planar boundary" begin
    n = normalize(randn(3))
    V = 7.0
    positions = [10 .* randn(3) for _ in 1:5]
    times = [dot(r, n) / V + 3 for r in positions]
    res = CVA(positions, times)
    @test res.n ≈ n
    @test res.V ≈ V
end

@testset "Reciprocal vector" begin
    r1 = [64278.6, 17683.5, -2512.4]
    r2 = [64276.4, 17704.7, -2514.06]
    r3 = [64291.1, 17697.7, -2514.05]
    r4 = [64282.6, 17697.4, -2525.76]
    rv = [-0.022940724814179882, -0.04814997115792451, 0.01788566629102226]
    @test reciprocal_vectors([r1, r2, r3, r4])[1] ≈ rv
end
