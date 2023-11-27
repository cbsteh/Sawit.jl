
"""
    Afgen{ITP <: AbstractArray}

Table of x vs. y values, with linear interpolation utility for y given x.

# Fields
- `x::Vector{Float64}`: list of x values
- `y::Vector{Float64}`: list of y values
- `itp<:AbstractArray`: linear interpolation function for y given x

# Notes
- No extrapolation.
"""
@with_kw struct Afgen{ITP <: AbstractArray}
    x::Vector{Float64} = []
    y::Vector{Float64} = []
    itp::ITP = []
end


"""
    Afgen(dct::AbstractDict)

Constructor for `Agfen` object.

# Arguments
- `dct::AbstractDict`: (x,y) data, where `x' are the dictionary keys,
                       and `y` are the dictonary values

# Returns
- `Afgen`: instantiated Afgen object
"""
function Afgen(dct::AbstractDict)
    lst = [(parse(Float64, k), v) for (k, v) ∈ pairs(dct)]
    sort!(lst)  # order by `x` for linear interpolation
    x, y = first.(lst), last.(lst)
    itp = interpolate((x,), y, Gridded(Linear()))   # use linear interpolation
    Afgen(x=x, y=y, itp=itp)
 end


"""
    (ag::Afgen)(x)

Function object to return the interpolated y for x.

# Arguments
- `x`: x value

# Notes
- No extrapolation.

# Returns
- `<:Real`: interpolated y for x
"""
(ag::Afgen)(x) = ag.itp(x)
