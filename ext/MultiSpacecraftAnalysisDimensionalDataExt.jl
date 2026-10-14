module MultiSpacecraftAnalysisDimensionalDataExt

using DimensionalData
using MultiSpacecraftAnalysis
using DimensionalData: TimeDim, dims, hasdim, dimnum, otherdims, lookup
import MultiSpacecraftAnalysis: lingradest, lingradest_metadata, MatrixSet

# A no-error version of `dimnum`
_dimnum(x, dim) = hasdim(x, dim) ? dimnum(x, dim) : nothing

const DimMatrixSet = Union{Tuple{Vararg{AbstractDimMatrix}}, AbstractVector{<:AbstractDimMatrix}}

# Component dims of an output element type, reusing the input's component dim `comp`;
# `grad[i, j] = ∂_i V_j` needs a distinct dim for the derivative direction i.
_compdims(T, comp) = T <: Number ? () : ndims(T) == 1 ? (comp,) : (Dim{:∂}(lookup(comp)), comp)

"""
    lingradest(Vs, Rs; dim = nothing, kw...)

`Vs` as tuples or vectors of `AbstractDimMatrix`: `dim` defaults to their time dimension. Returns a `DimStack`
of numeric arrays, time along `dim`.
"""
function MultiSpacecraftAnalysis.lingradest(Vs::DimMatrixSet, Rs::MatrixSet; dim = nothing, kw...)
    V1 = first(Vs)
    dim = @something dim _dimnum(V1, TimeDim) 1
    time = dims(V1, dim)
    comp = only(otherdims(V1, time))
    out = lingradest(map(parent, Vs), map(parent, Rs); dim, kw...)
    das = map(propertynames(out)) do name
        col = getproperty(out, name)
        comps = _compdims(eltype(col), comp)
        data, dims = if isempty(comps)
            col, (time,)
        elseif dim == 1
            stack(col; dims = 1), (time, comps...)
        else
            stack(col), (comps..., time)
        end
        DimArray(data, dims; name, metadata = get(lingradest_metadata, name, Dict()))
    end
    return DimStack(das...)
end

end
