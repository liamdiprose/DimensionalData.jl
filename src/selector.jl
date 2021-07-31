"""
    Selector

Abstract supertype for all selectors.

Selectors are wrappers that indicate that passed values are not the array indices,
but values to be selected from the dimension index, such as `DateTime` objects for
a `Ti` dimension.

Selectors provided in DimensionalData are:

- [`At`](@ref)
- [`Between`](@ref)
- [`Near`](@ref)
- [`Where`](@ref)
- [`Contains`](@ref)

"""
abstract type Selector{T} end

val(sel::Selector) = sel.val
rebuild(sel::Selector, val) = basetypeof(sel)(val)

@inline maybeselector(I...) = maybeselector(I)
@inline maybeselector(I::Tuple) = map(maybeselector, I)
# Int AbstractArray and Colon do normal indexing
@inline maybeselector(i::StandardIndices) = i
# Selectors are allready selectors
@inline maybeselector(i::Selector) = i
# Anything else becomes `At`
@inline maybeselector(i) = At(i)
@inline maybeselector(f::Fix2Selector) = Between(f)

"""
    At <: Selector

    At(x, atol, rtol)
    At(x; atol=nothing, rtol=nothing)

Selector that exactly matches the value on the passed-in dimensions, or throws an error.
For ranges and arrays, every intermediate value must match an existing value -
not just the end points.

`x` can be any value or `Vector` of values.

`atol` and `rtol` are passed to `isapprox`.
For `Number` `rtol` will be set to `Base.rtoldefault`, otherwise `nothing`,
and wont be used.

## Example

```jldoctest
using DimensionalData

A = DimArray([1 2 3; 4 5 6], (X(10:10:20), Y(5:7)))
A[X(At(20)), Y(At(6))]

# output

5
```
"""
struct At{T,A,R} <: Selector{T}
    val::T
    atol::A
    rtol::R
end
At(val; atol=nothing, rtol=nothing) =
    At{typeof.((val, atol, rtol))...}(val, atol, rtol)

atol(sel::At) = sel.atol
rtol(sel::At) = sel.rtol

"""
    Near <: Selector

    Near(x)

Selector that selects the nearest index to `x`.

With [`Points`](@ref) this is simply the index values nearest to the `x`,
however with [`Intervals`](@ref) it is the interval _center_ nearest to `x`.
This will be offset from the index value for `Start` and
[`End`](@ref) loci.

## Example

```jldoctest
using DimensionalData

A = DimArray([1 2 3; 4 5 6], (X(10:10:20), Y(5:7)))
A[X(Near(23)), Y(Near(5.1))]

# output

4
```
"""
struct Near{T} <: Selector{T}
    val::T
end

"""
    Contains <: Selector

    Contains(x)

Selector that selects the interval the value is contained by. If the
interval is not present in the index, an error will be thrown.

Can only be used for [`Intervals`](@ref) or [`Categorical`](@ref).

## Example

```jldoctest
using DimensionalData

dims_ = X(10:10:20; mode=Sampled(sampling=Intervals())),
        Y(5:7; mode=Sampled(sampling=Intervals()))
A = DimArray([1 2 3; 4 5 6], dims_)
A[X(Contains(8)), Y(Contains(6.8))]

# output

3
```
"""
struct Contains{T} <: Selector{T}
    val::T
end

struct PosInf end

Base.isless(::PosInf, x) = false
Base.isless(x, ::PosInf) = true

struct NegInf end

Base.isless(::NegInf, x) = true
Base.isless(x, ::NegInf) = false

"""
    Between <: Selector

    Between(a, b)

Selector that retreive all indices located between 2 values, evaluated
with `>=` for the lower value `a`, and `<` for the upper value `b`.
This means the same value will not be counted twice in 2 adjacent
`Between` selections.

For [`Intervals`](@ref) the whole interval must be lie between the
values. For [`Points`](@ref) the points must fall between
the values. Different [`Sampling`](@ref) types may give different
results with the same input - this is the intended behaviour.

`Between` for [`Irregular`](@ref) intervals is a little complicated. The
interval is the distance between a value and the next (for `Start` locus)
or previous (for [`End`](@ref) locus) value.

For [`Center`](@ref), we take the mid point between two index values
as the start and end of each interval. This may or may not make sense for
the values in your indes, so use `Between` with `Irregular` `Intervals(Center())`
with caution.

## Example

```jldoctest
using DimensionalData

A = DimArray([1 2 3; 4 5 6], (X(10:10:20), Y(5:7)))
A[X(Between(15, 25)), Y(Between(4, 6.5))]

# output

1×2 DimArray{Int64,2} with dimensions:
  X: 20:10:20 (Sampled - Ordered Regular Points),
  Y: 5:6 (Sampled - Ordered Regular Points)
 4  5
```
"""
struct Between{T<:Union{Tuple{Any,Any},Nothing}} <: Selector{T}
    val::T
end
Between(a, b) = Between((a, b))

# Between using Fix2
const LowerFix2Selector = Union{Fix2{typeof(>=),<:Any},Fix2{typeof(>),<:Any}}
const UpperFix2Selector = Union{Fix2{typeof(<=),<:Any},Fix2{typeof(<),<:Any}}
const Fix2Selector = Union{LowerFix2Selector,UpperFix2Selector}

# Between(t::Tuple{Any,Any}) = Between((>=(t[1]), <=(t[2])))
# Between(t::Tuple{Fix2Selector,Fix2Selector}) = Between{typeof(t)}((t[1], t[2]))
# Between(a::Fix2Selector, b::Fix2Selector) = Between(_comparison_order(a, b))
Between(a::Fix2Selector) = Between(_to_inf(a))

Base.first(sel::Between) = first(val(sel))
Base.last(sel::Between) = last(val(sel))

_comparison_order(a::LowerFix2Selector, b::UpperFix2Selector) = (a, b)
_comparison_order(a::UpperFix2Selector, b::LowerFix2Selector) = (b, a)
_comparison_order(a::Fix2Selector, b::Fix2Selector) =
    throw(Argumen("Use less than (<=,<) and/or greater than (>=,>) symbols, but not two of either"))

_to_inf(f::LowerFix2Selector) = (f, PosInf())
_to_inf(f::UpperFix2Selector) = (NegInf(), f)

# Dimension constructor passed Fix2 uses Between
(::Type{D})(f1::Fix2Selector, f2::Fix2Selector) where D<:Dimension = D(Between(f1, f2))

"""
    Where <: Selector

    Where(f::Function)

Selector that filters a dimension by any function that accepts
a single value from the index and returns a `Bool`.

## Example

```jldoctest
using DimensionalData

A = DimArray([1 2 3; 4 5 6], (X(10:10:20), Y(19:21)))
A[X(Where(x -> x > 15)), Y(Where(x -> x in (19, 21)))]

# output

1×2 DimArray{Int64,2} with dimensions:
  X: Int64[20] (Sampled - Ordered Regular Points),
  Y: Int64[19, 21] (Sampled - Ordered Regular Points)
 4  6
```
"""
struct Where{T} <: Selector{T}
    f::T
end

val(sel::Where) = sel.f


# sel2indices ==========================================================================

# Converts Selectors to regular indices

@inline sel2indices(x, ls...) = sel2indices(dims(x), ls...)
@inline sel2indices(dims::Tuple, l1, ls...) = sel2indices(dims, (l1, ls...))
@inline sel2indices(dims::Tuple, lookup::Tuple) =
    map((d, l) -> sel2indices(d, l), dims, lookup)
@inline sel2indices(dims::Tuple, lookup::Tuple{}) = ()
@inline sel2indices(dim::Dimension, sel) = _sel2indices(dim, maybeselector(sel))

# First filter based on rough selector properties -----------------

# Standard indices are just returned.
@inline _sel2indices(::Dimension, sel::StandardIndices) = sel
# Vectors are mapped
@inline _sel2indices(dim::Dimension, sel::Selector{<:AbstractVector}) =
    [_sel2indices(mode(dim), dim, rebuild(sel, v)) for v in val(sel)]
@inline _sel2indices(dim::Dimension, sel::Selector) = _sel2indices(mode(dim), dim, sel)

# Where selector ==============================
# Yes this is everything. Where doesn't need mode specialisation
@inline _sel2indices(dim::Dimension, sel::Where) =
    [i for (i, v) in enumerate(index(dim)) if sel.f(v)]

# Then dispatch based on IndexMode -----------------
# Selectors can have varied behaviours depending on the index mode.

# Noindex Contains just converts the selector to standard indices. Implemented
# so the Selectors actually work, not because what they do is useful or interesting.
@inline _sel2indices(mode::NoIndex, dim::Dimension, sel::Union{At,Near,Contains}) = val(sel)
@inline _sel2indices(mode::NoIndex, dim::Dimension, sel::Between{<:Tuple{Integer,Integer}}) =
    val(sel)[1]:val(sel)[2]
@inline _sel2indices(mode::Categorical, dim::Dimension, sel::Selector) =
    if sel isa Union{Contains,Near}
        _sel2indices(Points(), mode, dim, At(val(sel)))
    else
        _sel2indices(Points(), mode, dim, sel)
    end
@inline _sel2indices(mode::AbstractSampled, dim::Dimension, sel::Selector) =
    _sel2indices(sampling(mode), mode, dim, sel)

# For Sampled filter based on sampling type and selector -----------------

@inline _sel2indices(sampling::Sampling, mode::IndexMode, dim::Dimension, sel::At) =
    at(sampling, mode, dim, sel)
@inline _sel2indices(sampling::Sampling, mode::IndexMode, dim::Dimension, sel::Near) = begin
    span(mode) isa Irregular && locus(mode) isa Union{Start,End} && _nearirregularerror()
    near(sampling, mode, dim, sel)
end
@inline _sel2indices(sampling::Points, mode::IndexMode, dim::Dimension, sel::Contains) =
    _containspointserror()
@inline _sel2indices(sampling::Intervals, mode::IndexMode, dim::Dimension, sel::Contains) =
    contains(sampling, mode, dim, sel)
@inline _sel2indices(sampling::Sampling, mode::IndexMode, dim::Dimension, sel::Between{<:Tuple}) =
    between(sampling, mode, dim, sel)

@noinline _nearirregularerror() = throw(ArgumentError("Near is not implemented for Irregular with Start or End loci. Use Contains"))
@noinline _containspointserror() = throw(ArgumentError("`Contains` has no meaning with `Points`. Use `Near`"))


# Unaligned IndexMode ------------------------------------------

# unalligned2indices is callled directly from dims2indices

# We use the transformation from the first Transformed dim.
# In practice the others could be empty.
@inline unalligned2indices(dims::DimTuple, sel::Tuple) = sel
@inline unalligned2indices(dims::DimTuple, sel::Tuple{<:Dimension,Vararg{<:Dimension}}) =
    unalligned2indices(dims, map(val, sel))
@inline unalligned2indices(dims::DimTuple, sel::Tuple{<:Selector,Vararg{<:Selector}}) = begin
    coords = [map(val, sel)...]
    transformed = transformfunc(mode(dims[1]))(coords)
    map(_transform2int, sel, transformed)
end

_transform2int(::At, x) = convert(Int, x)
_transform2int(::Near, x) = round(Int, x)


# Internal use structs --------------------------------------

struct _Upper end
struct _Lower end

# Selector methods

# at =============================================================================

at(dim::Dimension, sel::At) = at(sampling(mode(dim)), mode(dim), dim, sel)
function at(::Sampling, mode::IndexMode, dim::Dimension, sel::At)
    relate(dim, at(indexorder(dim), val(dim), dim, val(sel), atol(sel), rtol(sel)))
end
#
function at(o::IndexOrder, a::AbstractArray{<:Union{Number,Dates.TimeType}}, dim::Dimension, selval, atol, rtol::Nothing)
    x = unwrap(selval)
    i = searchsortedlast(a, x; order=_ordering(o))
    _isat(x, a[i], atol, 1) || _selvalnotfound(dim, selval)
    return i
end
# range optimisation - we leave bounds checking to the array, but maybe this is not cool
function at(o::IndexOrder, a::AbstractRange{<:Number}, dim::Dimension, selval, atol, rtol::Nothing)
    x = unwrap(selval)
    f = (x - first(a)) / step(a)
    i = unsafe_trunc(Int, f)
    _isat(i, f, atol, step(a)) && checkbounds(Bool, a, i + 1) || _selvalnotfound(dim, selval)
    return i + 1
end

@inline _isat(i, f, atol::Nothing, scalar) = i == f
@inline _isat(i, f, atol, scalar) = abs(i - f) <= atol / scalar

# catch-all for an unordered or non-number index
function at(
    ::IndexOrder, a, dim::Dimension, selval, atol, rtol::Nothing
)
    i = findfirst(x -> x == unwrap(selval), index(dim))
    i == nothing && _selvalnotfound(dim, selval)
    return i
end
# compile-time indexing
@generated function at(
    ::IndexOrder, ::Val{Index}, dim::Dimension, selval::Val{X}, atol, rtol::Nothing
) where {Index,X}
    i = findfirst(x -> x == X, Index)
    if i == nothing
        :(_selvalnotfound(dim, selval))
    else
        return i
    end
end

@noinline _selvalnotfound(dim, selval) = throw(ArgumentError("$selval not found in $(name(dim))"))


# near ===========================================================================

# Finds the nearest point in the index, adjusting for locus if necessary.
# In Intevals we are finding the nearest point to the center of the interval.

near(dim::Dimension, sel::Near) = near(sampling(mode(dim)), mode(dim), dim, sel)
function near(::Sampling, mode::IndexMode, dim::Dimension, sel::Near)
    order = indexorder(dim)
    order isa UnorderedIndex && _nearunorderederror()
    locus = DD.locus(dim)

    lower_dimside = _dimside(_Lower(), order, dim)
    v = _locus_adjust(locus, unwrap(val(sel)), dim)
    i = _inbounds(_near_search_func(order)(order, dim, v), dim)

    islower = if order isa ForwardIndex
        i <= lower_dimside
    else
        i >= lower_dimside
    end

    i = if islower
        lower_dimside
    else
        previ = _prevind(order, i)
        vl, vi = map(abs, (dim[previ] - v, dim[i] - v))
        # We have to use the right >/>= for Start/End locus
        _lt(locus)(vl, vi) ? previ : i
    end
    relate(dim, i)
end

_near_search_func(::ForwardIndex) = _searchfirst
_near_search_func(::ReverseIndex) = _searchlast

_locus_adjust(locus::Center, v, dim) = v
_locus_adjust(locus::Start, v, dim) = v - abs(step(dim)) / 2
_locus_adjust(locus::End, v, dim) = v + abs(step(dim)) / 2
_locus_adjust(locus::Start, v::DateTime, dim) = v - (v - (v - abs(step(dim)))) / 2
_locus_adjust(locus::End, v::DateTime, dim) = v + (v + abs(step(dim)) - v) / 2

@noinline _nearunorderederror() = throw(ArgumentError("`Near` has no meaning in an `Unordered` index"))

# contains ================================================================================

# Finds which interval contains a point

function contains(dim::Dimension, sel::Contains)
    contains(sampling(mode(dim)), mode(dim), dim, sel)
end
# Points --------------------------------------
@noinline function contains(::Points, ::IndexMode, dim::Dimension, sel::Contains)
    throw(ArgumentError("Points IndexMode cannot use 'Contains', use 'Near' instead."))
end
# Intervals -----------------------------------
function contains(sampling::Intervals, mode::IndexMode, dim::Dimension, sel::Contains)
    relate(dim, contains(span(mode), sampling, indexorder(mode), locus(mode), dim, sel))
end
# Regular Intervals ---------------------------
function contains(span::Regular, ::Intervals, order, locus, dim::Dimension, sel::Contains)
    v = val(sel); s = abs(val(span))
    _locus_checkbounds(locus, bounds(dim), v)
    i = _contains_search_func(locus, order)(order, dim, _maybeaddhalf(locus, s, v))
    # Check the value is in this cell.
    # It is always for AbstractRange but might not be for Val tuple or Vector.
    if !(val(dim) isa AbstractRange) 
        _lt(locus)(v, dim[i] + s) || _notcontainederror(v)
    end
    i
end
# Explicit Intervals ---------------------------
function contains(span::Explicit, ::Intervals, order, locus, dim::Dimension, sel::Contains)
    x = val(sel)
    i = searchsortedlast(view(val(span), 1, :), x; order=_ordering(order))
    if i <= 0 || val(span)[2, i] < x 
        _notcontainederror(x)
    end
    i
end
# Irregular Intervals -------------------------
function contains(span::Irregular, ::Intervals, order::IndexOrder, locus::Locus, dim::Dimension, sel::Contains)
    _locus_checkbounds(locus, bounds(span), val(sel))
    _contains_search_func(locus, order)(order, dim, val(sel))
end
function contains(span::Irregular, ::Intervals, order::IndexOrder, locus::Center, dim::Dimension, sel::Contains)
    v = val(sel)
    _locus_checkbounds(locus, bounds(span), v)
    i = _searchfirst(order, dim, v)
    i <= firstindex(dim) && return firstindex(dim)
    i > lastindex(dim) && return lastindex(dim)

    interval = abs(dim[i] - dim[i - 1])
    distance = abs(dim[i] - v)
    _contains_lt(order)(interval / 2, distance) ? i - 1 : i
end

@noinline _notcontainederror(v) = throw(ArgumentError("No interval contains $(v)"))

_contains_search_func(::Locus, ::ForwardIndex) = _searchlast
_contains_search_func(::Locus, ::ReverseIndex) = _searchfirst
_contains_search_func(::End, ::ForwardIndex) = _searchfirst
_contains_search_func(::End, ::ReverseIndex) = _searchlast

_maybeaddhalf(::Locus, s, v) = v
_maybeaddhalf(::Center, s, v) = v + s / 2

_contains_lt(::ForwardIndex) = (<)
_contains_lt(::ReverseIndex) = (<=)


# between ================================================================================

# between => AbstractRange
# Finds all values between two points, adjusted for locus where necessary
between(dim::Dimension, sel::Between) = between(sampling(mode(dim)), mode(dim), dim, sel)
function between(sampling::Sampling, mode::IndexMode, dim::Dimension, sel::Between)
    sel = _to_fix2(locus(dim), sel)
    order = indexorder(dim)
    # Only accept an ordered index
    order isa UnorderedIndex && throw(ArgumentError("Cannot use `Between` with UnorderedIndex"))
    # Get lower and higher indices of `dim` between `sel` values
    return between(sampling, order, mode, dim, sel)
end

#= Points ------------------------------------
Points are selected

# Example with Start locus
Bounds ┃        ┃
Lower  |  |  |  |
Upper     |  |  |  
=#
function between(sampling::Points, o::IndexOrder, ::IndexMode, d::Dimension, sel::Between)
    # Search for the upper and lower indices
    I = map(val(sel)) do f
        _search(o, d, f)
    end
    return _finish_between(o, d, I)
end

#= Intervals -------------------------
Select whole intervals falling between two bounds, 
using the given less-than and greater-than operators
=#
function between(sampling::Intervals, o::IndexOrder, mode::IndexMode, dim::Dimension, sel::Between)
    between(span(mode), sampling, o, mode, dim, sel)
end

#= Regular Intervals -------------------------
Bounds come from the index so that the upper bound of an interval is the 
lower bound of the next interval. The interval size is regular.

# Example with Start locus
Bounds ┃          ┃
Lower  |  |  |  |
Upper     |  |  |  
Index  1- 2- 3- 4-
=#
function between(span::Regular, ::Intervals, o::IndexOrder, mode::IndexMode, d::Dimension, sel::Between)
    # Search for the upper and lower indices
    I = map(val(sel), _interval_bounds_diffs(mode)) do f, ibs
        # Rebuild Fix2 with value adjusted to be the interval bound
        f = f.f(f.x + ibs)
        # Search in dim
        _search(o, d, f)
    end
    return _finish_between(o, d, I)
end

#= Irregular Intervals -----------------------
Bounds come from the index so that the upper bound of an interval is the 
lower bound of the next interval. The interval size may be irregular.

Example with Start locus:
       
Bounds ┃           ┃
Lower  |  |  |   |  
Upper     |  |   |   
Index  1- 2- 3-  4-
=#
function between(span::Irregular, ::Intervals, o::IndexOrder, mode::IndexMode, d::Dimension, sel::Between)
    I = map(val(sel), bounds(span)) do f, b
        # Check if the search value is within the bounds
        if f(b)
            _dimside(f, o, d)
        else
            # In bounds so search for the index
            _irreg_search(locus(mode), o, d, f, b)
        end
    end
    return _finish_between(o, d, I)
end

#= Explicit Intervals -------------------------
Lower and upper interval bounds come from an explicit bounds matrix. They usually line up exactly 
and are simply duplicated with an offset of 1 in lower and upper rows, but they may not be:

Example with Start locus and gap in intervals:

Bounds ┃      ┋    ┃
Lower  |  |   |  | 
Upper     |  |   | |  
Index  1- 2-  3- 4- 
We store the intervals in a matrix and take views of it to get the lower and upper bounds vectors
=#
function between(span::Explicit, ::Intervals, o::IndexOrder, mode::IndexMode, d::Dimension, sel::Between)
    # Use the bounds matrix as the vector index for both lower and upper bounds values
    lowerbounds = rebuild(d, view(val(span), 1, :))
    upperbounds = rebuild(d, view(val(span), 2, :))
    f1, f2 = val(sel)
    I = _search(o, lowerbounds, f1), _search(o, upperbounds, f2)
    return _finish_between(o, d, I)
end

function _finish_between(o, d, I)
    # Maybe reverse the indices for a reverse index
    l, h = _maybereverse(o, I)
    # Make sure indices are in bounds
    return _maybeflip(relation(d), d, l:h)
end

function _to_fix2(locus, sel::Between) 
    a = _lower_fix2(locus, first(sel))
    b = _upper_fix2(locus, last(sel)) 
    Between(a, b)
end
# We always include the index value in the interval by 
# using <= or >=. `Center` preferences the lower mid-point.
_lower_fix2(::Locus, val::Fix2) = val 
_lower_fix2(::Locus, val) = >=(val)
_upper_fix2(::Locus, val::Fix2) = val
_upper_fix2(::Locus, val) = <=(val)

# _interval_bounds_diffs => Tuple{<:Any,<:Any}
# Calculate the amounts differences between an index value to lower and upper interval bounds
_interval_bounds_diffs(mode) = _interval_bounds_diffs(locus(mode), abs(step(span(mode))))
_interval_bounds_diffs(locus::Start, step) = zero(step), -step
_interval_bounds_diffs(locus::Center, step) = step/2, -step/2
_interval_bounds_diffs(locus::End, step) = step, zero(step)

_search(o::ForwardIndex, dim, f::LowerFix2Selector) = _searchfirst(o, dim, f.x, _fix2_lt(f))
_search(o::ForwardIndex, dim, f::UpperFix2Selector) = _searchlast(o, dim, f.x, _fix2_lt(f))
_search(o::ReverseIndex, dim, f::LowerFix2Selector) = _searchlast(o, dim, f.x, _fix2_lt(f))
_search(o::ReverseIndex, dim, f::UpperFix2Selector) = _searchfirst(o, dim, f.x, _fix2_lt(f))

# We need to convert our Fix2 functions to less-than as used in searchsorted functions
# searchsortedfirst lt=< means first value >= ; order=_ordering(o), lt=<
_fix2_lt(::Fix2{typeof(>=)}) = <
_fix2_lt(::Fix2{typeof(>)}) = <=
# searchsortedlast why do we swap this
_fix2_lt(::Fix2{typeof(<=)}) = <
_fix2_lt(::Fix2{typeof(<)}) = <=

# _irreg_search
# Searching an irregular index for interval bounds is a little more difficult as we
# have to calculate them from adjactent values or the midpoint between two values.
# This has to be done after we find the closest index, rather than beforehand.
function _irreg_search(locus::Center, o::IndexOrder, d::Dimension, f::LowerFix2Selector, bound)
    i = _search(o, d, f)
    os = _order_shift(o)
    ilowerbound = i - os 
    interval_bound = if ilowerbound < firstindex(d) || ilowerbound > lastindex(d)
        # Use the array boundary as the interval boundary
        bound
    else
        # Use the half way point between index values as the interval boundary
        d[i] + (d[ilowerbound] - d[i]) / 2 
    end
    interval_is_inside_bound = f(interval_bound)
    return interval_is_inside_bound ? i : i + os
end
function _irreg_search(locus::Center, o::IndexOrder, d::Dimension, f::UpperFix2Selector, bound)
    i = _search(o, d, f)
    os = _order_shift(o)
    iupperbound = i + os 
    interval_bound = if iupperbound > lastindex(d) || iupperbound < firstindex(d)
        # Use the array boundary as the interval boundary
        bound
    else
        # Use the half way point between index values as the interval boundary
        d[i] + (d[iupperbound] - d[i]) / 2
    end
    interval_is_inside_bound = f(interval_bound)
    return interval_is_inside_bound ? i : i - os
end
function _irreg_search(locus, o::IndexOrder, d::Dimension, f::Fix2Selector, bound)
    _search(o, d, f) + _order_shift(o) * (_locus_shift(locus) + _side_shift(f))
end

# Index shuffling for irregular interval search results
_locus_shift(::Start) = 0
_locus_shift(::End) = 1
_side_shift(::LowerFix2Selector) = 0
_side_shift(::UpperFix2Selector) = -1
_order_shift(::ForwardIndex) = 1
_order_shift(::ReverseIndex) = -1

_lt(::typeof(<)) = (<)
_lt(::typeof(<=)) = (<=)
_lt(::typeof(>)) = (<=)
_lt(::typeof(>=)) = (<)
_lt(::_Lower) = (<)
_lt(::_Upper) = (<=)


# Shared utils ============================================================================

# _searchlast => Int
# searchsortedlast for a dimension
_searchlast(o::IndexOrder, dim::Dimension, v, lt=<) =
    searchsortedlast(index(dim), unwrap(v); order=_ordering(o), lt=lt)
_searchlast(o::IndexOrder, dim::Dimension{<:Val{Index}}, v, lt=<) where Index =
    searchsortedlast(Index, unwrap(v); order=_ordering(o), lt=lt)

# _searchfirst => Int
# searchsortedfirst for a dimension
_searchfirst(o::IndexOrder, dim::Dimension, v, lt=<) =
    searchsortedfirst(index(dim), unwrap(v); order=_ordering(o), lt=lt)
_searchfirst(o::IndexOrder, dim::Dimension{<:Val{Index}}, v, lt=<) where Index =
    searchsortedfirst(Index, unwrap(v); order=_ordering(o), lt=lt)

_asfunc(::Type{typeof(<)}) = <
_asfunc(::Type{typeof(<=)}) = <=

# Return an inbounds index
_inbounds(is::Tuple, dim::Dimension) = map(i -> _inbounds(i, dim), is)
_inbounds(i::Int, dim::Dimension) =
    if i > lastindex(dim)
        lastindex(dim)
    elseif i <= firstindex(dim)
        firstindex(dim)
    else
        i
    end

_sorttuple(sel::Between) = _sorttuple(val(sel))
_sorttuple((a, b)) = a < b ? (a, b) : (b, a)

_maybereverse(o::ForwardIndex, (a, b)) = (a, b)
_maybereverse(o::ReverseIndex, (a, b)) = (b, a)

_lt(::Locus) = (<)
_lt(::End) = (<=)
_gt(::Locus) = (>=)
_gt(::End) = (>)

_locus_checkbounds(loc, (l, h), v) =
    (_lt(loc)(v, l) || _gt(loc)(v, h)) && throw(BoundsError())

_prevind(::ForwardIndex, i) = i - 1
_prevind(::ReverseIndex, i) = i + 1

_dimside(::Union{_Lower,LowerFix2Selector}, o::ForwardIndex, d) = firstindex(d)
_dimside(::Union{_Lower,LowerFix2Selector}, o::ReverseIndex, d) = lastindex(d)
_dimside(::Union{_Upper,UpperFix2Selector}, o::ForwardIndex, d) = lastindex(d)
_dimside(::Union{_Upper,UpperFix2Selector}, o::ReverseIndex, d) = firstindex(d)

_ordering(::ForwardIndex) = Base.Order.ForwardOrdering()
_ordering(::ReverseIndex) = Base.Order.ReverseOrdering()
