module MultiSpacecraftAnalysisDimensionalDataExt

using DimensionalData
using MultiSpacecraftAnalysis
using DimensionalData: TimeDim, dims, hasdim, dimnum
import MultiSpacecraftAnalysis: lingradest, lingradest_metadata

# A no-error version of `dimnum`
_dimnum(x, dim) = hasdim(x, dim) ? dimnum(x, dim) : nothing

"""
    lingradest(B1::AbstractDimArray, args...)

Method for handling dimensional arrays. Takes `AbstractDimArray` inputs with a time dimension
and returns a `DimStack` containing all computed quantities.
"""
function MultiSpacecraftAnalysis.lingradest(B1::AbstractDimMatrix, args::AbstractDimMatrix...; flatten = false, dim = nothing, kw...)
    dim = @something dim _dimnum(B1, TimeDim) 1
    time = dims(B1, dim)
    out = lingradest(parent(B1), parent.(args)...; dim, flatten, kw...)
    names = propertynames(out)
    das = map(names) do name
        data = getproperty(out, name)
        dims = if flatten && ndims(data) > 1
            dim == 1 ? (time, Y(axes(data, 2))) : (Y(axes(data, 1)), time)
        else
            (time,)
        end
        metadata = get(lingradest_metadata, name, Dict())
        DimArray(data, dims; name, metadata)
    end
    return DimStack(das...)
end

end
