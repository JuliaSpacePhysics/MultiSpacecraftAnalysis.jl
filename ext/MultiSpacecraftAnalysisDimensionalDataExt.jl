module MultiSpacecraftAnalysisDimensionalDataExt

using DimensionalData
using MultiSpacecraftAnalysis
import MultiSpacecraftAnalysis: _keyed, _keyed_fields, lingradest_metadata

# `rebuild` keeps the name and metadata (e.g. units) of the first field, which describe none of the outputs.
function _keyed(V1::AbstractDimArray, out, d)
    fields = _keyed_fields(V1, out, d)
    return DimStack(NamedTuple{keys(fields)}(map(keys(fields), values(fields)) do name, A
        DimensionalData.rebuild(A; name, metadata = lingradest_metadata[name])
    end))
end

end
