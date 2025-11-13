"""
Transform matrix-like `A` to `n×m` shape
"""
function ensure_nxm(A, n, m)
    if size(A, 1) == m && size(A, 2) == n
        return permutedims(A)  # Convert from m×n to n×m
    elseif size(A, 1) == n && size(A, 2) == m
        return A  # Already in n×m format
    else
        throw(ArgumentError("A must be either n×m or m×n matrices, but got size $(size(A))"))
    end
end
