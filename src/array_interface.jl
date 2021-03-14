ArrayInterface.dimnames(A::AbstractDimArray) = 
    map(ArrayInterface.StaticSymbol ∘ dim2key, dims(A))
