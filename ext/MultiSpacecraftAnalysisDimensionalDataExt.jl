module MultiSpacecraftAnalysisDimensionalDataExt

using DimensionalData
using MultiSpacecraftAnalysis
import MultiSpacecraftAnalysis: lingradest

"""
    lingradest(B1::AbstractDimArray, args...)

Method for handling dimensional arrays. Takes `AbstractDimArray` inputs with a time dimension
and returns a `DimStack` containing all computed quantities.
"""
function MultiSpacecraftAnalysis.lingradest(B1::AbstractDimArray, args...)
    time = dims(B1, Ti)
    out = lingradest(parent(B1), parent.(args)...)
    names = propertynames(out)
    das = map(out, names) do p, name
        DimArray(p, time; name = name)
    end
    return DimStack(das...)
end

end
