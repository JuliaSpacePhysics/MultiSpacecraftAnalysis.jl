module MultiSpacecraftAnalysisDimensionalDataExt

using DimensionalData
using MultiSpacecraftAnalysis
import MultiSpacecraftAnalysis: _keyed, _keyed_fields, lingradest_metadata

# `rebuild` keeps the name and metadata (e.g. units) of the first field, which describe none of the outputs.
function _keyed(V1::AbstractDimArray, out)
    fields = _keyed_fields(V1, out)
    return DimStack(map(keys(fields), values(fields)) do name, A
        DimensionalData.rebuild(A; name, metadata = get(lingradest_metadata, name, Dict()))
    end...)
end

end
