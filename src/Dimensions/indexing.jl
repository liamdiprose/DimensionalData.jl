#### dims2indices ####

"""
    dims2indices(dim::Dimension, I) => NTuple{Union{Colon,AbstractArray,Int}}

Convert a `Dimension` or `Selector` `I` to indices of `Int`, `AbstractArray` or `Colon`.
"""
@inline dims2indices(dim::Dimension, I::StandardIndices) = I
@inline dims2indices(dim::Dimension, I) = _dims2indices(dim, I)

@inline dims2indices(x, I) = dims2indices(dims(x), I)
@inline dims2indices(::Nothing, I) = _dimsnotdefinederror()
@inline dims2indices(dims::DimTuple, I) = dims2indices(dims, (I,))
# Standard array indices are simply returned
@inline dims2indices(dims::DimTuple, I::Tuple{Vararg{<:StandardIndices}}) = I
@inline dims2indices(dims::DimTuple, I::Tuple{<:CartesianIndex}) = I
@inline dims2indices(dims::DimTuple, sel::Tuple) = 
    LookupArrays.selectindices(lookup(dims), sel)
@inline dims2indices(dims::DimTuple, ::Tuple{}) = ()
# Otherwise attempt to convert dims to indices
@inline function dims2indices(dims::DimTuple, I::DimTuple)
    extradims = otherdims(I, dims)
    length(extradims) > 0 && _warnextradims(extradims)
    _dims2indices(lookup(dims), dims, sortdims(I, dims))
end

# Handle tuples with @generated
@inline _dims2indices(::Tuple{}, dims::Tuple{}, ::Tuple{}) = ()
@generated function _dims2indices(lookups::Tuple, dims::Tuple, I::Tuple)
    # We separate out Aligned and Unaligned lookups as
    # Unaligned must be selected in groups e.g. X and Y together.
    unalligned = Expr(:tuple)
    uaI = Expr(:tuple)
    alligned = Expr(:tuple)
    dimmerge = Expr(:tuple)
    a_count = ua_count = 0
    for (i, lkup) in enumerate(lookups.parameters)
        if lkup <: Unaligned
            ua_count += 1
            push!(unalligned.args, :(dims[$i]))
            push!(uaI.args, :(I[$i]))
            push!(dimmerge.args, :(uadims[$ua_count]))
        else
            a_count += 1
            push!(alligned.args, :(_dims2indices(dims[$i], I[$i])))
            # Update the merged tuple
            push!(dimmerge.args, :(adims[$a_count]))
        end
    end

    if length(unalligned.args) > 1
        # Output the `dimmerge` tuple, which will
        # combine uadimsand adims in the right order
        quote
             adims = $alligned
             # Unaligned dims have to be run together as a set
             uadims = unalligned_dims2indices($unalligned, map(_unwrapdim, $uaI))
             $dimmerge
        end
    else
        alligned
    end
end

@inline function unalligned_dims2indices(dims::DimTuple, sel::Tuple)
    map(sel) do s
        s isa Selector && throw(ArgumentError("Unalligned dims: use selectors for all $(join(map(string ∘ dim2key, dims), ", ")) dims, or none of them"))
        isnothing(s) ? Colon() : s
    end
end
@inline function unalligned_dims2indices(dims::DimTuple, sel::Tuple{<:Selector,Vararg{<:Selector}})
    LookupArrays.select_unalligned_indices(lookup(dims), sel)
end

_unwrapdim(dim::Dimension) = val(dim)
_unwrapdim(x) = x

# Single dim methods
# A Dimension type always means Colon(), as if it was constructed with the default value.
@inline _dims2indices(dim::Dimension, ::Type{<:Dimension}) = Colon()
# Nothing means nothing was passed for this dimension
@inline _dims2indices(dim::Dimension, ::Nothing) = Colon()
# Simply unwrap dimensions
@inline _dims2indices(dim::Dimension, seldim::Dimension) = 
    LookupArrays.selectindices(val(dim), val(seldim))
