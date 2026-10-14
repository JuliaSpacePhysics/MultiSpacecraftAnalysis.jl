module MultiSpacecraftAnalysisDimensionalDataExt

using DimensionalData
using MultiSpacecraftAnalysis
using DimensionalData: TimeDim, dims, hasdim, dimnum
import MultiSpacecraftAnalysis: lingradest, lingradest_metadata, MatrixSet

# A no-error version of `dimnum`
_dimnum(x, dim) = hasdim(x, dim) ? dimnum(x, dim) : nothing

const DimMatrixSet = Union{Tuple{Vararg{AbstractDimMatrix}}, AbstractVector{<:AbstractDimMatrix}}

"""
    lingradest(Bs::DimMatrixSet, Rs; dim = nothing, kw...)

`Bs` as `AbstractDimMatrix`es: `dim` defaults to their time dimension, and the result is a `DimStack`.
"""
function MultiSpacecraftAnalysis.lingradest(Bs::DimMatrixSet, Rs::MatrixSet; flatten = false, dim = nothing, kw...)
    B1 = first(Bs)
    dim = @something dim _dimnum(B1, TimeDim) 1
    time = dims(B1, dim)
    out = lingradest(map(parent, Bs), map(parent, Rs); dim, flatten, kw...)
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
