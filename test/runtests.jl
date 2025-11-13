using MultiSpacecraftAnalysis
using Test
using Aqua
using PythonCall

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(MultiSpacecraftAnalysis)
end

@testset "Integration test" begin
    # Reference: https://beforerr.github.io/beforerr/docs/courses/epss261/homework/ps4.html

end

# skip nan
_is_approx(x, y; kw...) = isnan(x) && isnan(y) || isapprox(x, y; kw...)
_is_approx(x::AbstractArray, y::AbstractArray; kw...) = all(_is_approx.(x, y; kw...))

@testset "Cross validation with PySPEDAS.jl" begin
    # Reference: https://github.com/spedas/pyspedas/blob/master/pyspedas/projects/mms/tests/test_mms_curlometer.py#L4
    # https://github.com/spedas/pyspedas/blob/283bc4ce21a302af3b894451b489a6175dd8fd5d/pyspedas/projects/mms/fgm_tools/mms_lingradest.py#L8
    using PySPEDAS

    @py import pyspedas.projects.mms.tests.test_mms_curlometer: CurlTestCases
    @py import pyspedas.projects.mms.fgm_tools.mms_lingradest: mms_lingradest
    pytest = CurlTestCases()
    pytest.test_lingradest()


    field_names = ["mms1_fgm_b_gse_brst_l2", "mms2_fgm_b_gse_brst_l2_i", "mms3_fgm_b_gse_brst_l2_i", "mms4_fgm_b_gse_brst_l2_i"]
    position_names = ["mms1_fgm_r_gse_brst_l2_i", "mms2_fgm_r_gse_brst_l2_i", "mms3_fgm_r_gse_brst_l2_i", "mms4_fgm_r_gse_brst_l2_i"]

    func(x) = get_data(x)[:, 1:3]
    fields = func.(field_names)
    positions = func.(position_names)

    res = lingradest(fields..., positions...)

    @test _is_approx(res.Bbc, get_data("Bbc_lingradest")[:, 2:4])
    @test _is_approx(res.Bmag, get_data("Bt_lingradest"))
    scale = 1000
    @test _is_approx(res.LGBx .* scale, get_data("gradBx_lingradest"))
    @test _is_approx(res.LGBy .* scale, get_data("gradBy_lingradest"))
    @test _is_approx(res.LGBz .* scale, get_data("gradBz_lingradest"))
    @test _is_approx(res.div .* scale, get_data("divB_nT/1000km_lingradest"))
    @test _is_approx(res.curl .* scale, get_data("curlB_nT/1000km_lingradest")[:, 2:4])
    @test _is_approx(res.curv .* scale, get_data("curvB_lingradest"))
    @test _is_approx(res.R_c ./ scale, get_data("Rc_1000km_lingradest"))
    @test _is_approx(res.curl .* scale .* 0.8, get_data("jtotal_lingradest"))
end


@testset "Reciprocal vector" begin
    r1 = [64278.6, 17683.5, -2512.4]
    r2 = [64276.4, 17704.7, -2514.06]
    r3 = [64291.1, 17697.7, -2514.05]
    r4 = [64282.6, 17697.4, -2525.76]
    rv = [-0.022940724814179882, -0.04814997115792451, 0.01788566629102226]
    @test reciprocal_vector(r1, r2, r3, r4) == rv
    @test reciprocal_vector(r1, [r2, r3, r4]) == reciprocal_vector(r1, (r2, r3, r4)) ≈ rv
end
