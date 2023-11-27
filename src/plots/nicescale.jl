# Internal use
function nicenum(rng, to_round::Bool)
    exponent = floor(log10(rng))
    fraction = rng / 10^exponent
    fn = to_round ? [<(1.5), <(3), <(7)] : [<=(1), <=(2), <=(5)]
    nicefraction = fn[1](fraction) ? 1 :
                   fn[2](fraction) ? 2 :
                   fn[3](fraction) ? 5 : 10
    nicefraction * 10^exponent
end


# Internal use: Nice no. of decimal places for ticks
function nicedecimal(tickspacing)
    logtickx = log10(tickspacing)
    (logtickx >= 0) ? 0 : Int(abs(floor(logtickx)))
end


"""
    niceticks(minpt, maxpt, maxticks::Int=10)

Return an array of ticks for an axis.

# Arguments
- `minpt`: smallest value to be plotted
- `maxpt`: largest value to be plotted
- `maxticks::Int=10`: no. of ticks

# Returns
- `Vector{Float64}`: array of ticks. Length of array is `maxticks`.
"""
function niceticks(minpt, maxpt, maxticks::Int=10)
    maxticks = max(2, maxticks)
    rng = nicenum(maxpt - minpt, false)
    step = nicenum(rng / (maxticks - 1), true)      # tick spacing
    nicemin = floor(minpt / step) * step
    nicemax = ceil(maxpt / step) * step
    step = (nicemax - nicemin) / (maxticks - 1)     # nice tick spacing
    ar = [nicemin + i * step for i ∈ 0:maxticks-1]
    ar = round.(ar, digits=nicedecimal(step))

    # after rounding, tick spacings may not be equal
    # if not equal, retry fit but with one fewer tick count:
    Δ = ar[2] - ar[1]
    res = map(i -> isapprox(ar[i] - ar[i-1], Δ), 3:length(ar))
    if !all(res)
        ar = niceticks(minpt, maxpt, maxticks-1)
    end

    ar
end


"""
    niceticks(vals::AbstractVector, maxticks::Int=10)

Return an array of ticks for an axis.

# Arguments
- `vals::AbstractVector`: values to be plotted
- `maxticks::Int=10`: no. of ticks

# Returns
- `Vector{Float64}`: array of ticks. Length of array is `maxticks`.
"""
function niceticks(vals::AbstractVector, maxticks::Int=10)
    allvals = collect(Iterators.flatten(vals))
    minpt = minimum(allvals)
    maxpt = maximum(allvals)
    niceticks(minpt, maxpt, maxticks)
end
