"""
    BinXY

One bin (group) of x and y values.

# Fields
- `xkeys::Vector{Float64}`: start and end values of x for bin
- `xvals::Vector{Float64}`: list of x values in bin
- `yvals::Vector{Float64}`: list of y values in bin
"""
@with_kw mutable struct BinXY
    xkeys::Vector{Float64} = []
    xvals::Vector{Float64} = []
    yvals::Vector{Float64} = []
end


"""
    BinXYZ

One bin (group) of x, y, and z values.

# Fields
- `xkeys::Vector{Float64}`: start and end values of x for bin
- `ykeys::Vector{Float64}`: start and end values of y for bin
- `xvals::Vector{Float64}`: list of x values in bin
- `yvals::Vector{Float64}`: list of y values in bin
- `zvals::Vector{Float64}`: list of z values in bin
"""
@with_kw mutable struct BinXYZ
    xkeys::Vector{Float64} = []
    ykeys::Vector{Float64} = []
    xvals::Vector{Float64} = []
    yvals::Vector{Float64} = []
    zvals::Vector{Float64} = []
end


"""
     binXYZ(x::AbstractVector, y::AbstractVector, z::AbstractVector)

Digitize continuous 3D data (x, y, and z) into bins (or groups).

# Arguments
- `x::AbstractVector`: list of x values
- `y::AbstractVector`: list of y values
- `z::AbstractVector`: list of z values

# Returns
- `Vector{BinXYZ}`: bins (groups) of x, y, and z values.
"""
function binXYZ(x::AbstractVector, y::AbstractVector, z::AbstractVector)

    function in_group(bin, xval, yval)
        x1, x2 = bin.xkeys[1], bin.xkeys[2]
        y1, y2 = bin.ykeys[1], bin.ykeys[2]
        (x1 <= xval < x2) && (y1 <= yval < y2)
     end

    n = length(x)
    min_x, max_x = extrema(x)
    min_y, max_y = extrema(y)

    # mutual information binning:
    p2 = cor(x, y) ^ 2
    lx = (max_x - min_x) / std(x)
    ly = (max_y - min_y) / std(y)
    L = (n * p2) / (12 * (1 - p2)) * (lx^2 + ly^2)
    k = round(Int, 0.5 + 0.5 * sqrt(1 + 4 * sqrt(L)))

    Δx = (max_x - min_x) / k
    Δy = (max_y - min_y) / k
    δx = Δx / 2
    δy = Δy / 2

    bins = BinXYZ[]
    for xv ∈ min_x:Δx:max_x, yv ∈ min_y:Δy:max_y
        push!(bins, BinXYZ(xkeys=[xv-δx, xv+δx], ykeys=[yv-δy, yv+δy]))
    end

    foreach(zip(x, y, z)) do xyz
        idx = findfirst(b -> in_group(b, xyz[1], xyz[2]), bins)
        if !isnothing(idx)
            b = bins[idx]
            push!(b.xvals, xyz[1])
            push!(b.yvals, xyz[2])
            push!(b.zvals, xyz[3])
        end
    end

    bins
end


"""
     binXY(x::AbstractVector, y::AbstractVector)

Digitize continuous 2D data (x and y) into bins (or groups).

# Arguments
- `x::AbstractVector`: list of x values
- `y::AbstractVector`: list of y values

# Returns
- `Vector{BinXY}`: bins (groups) of x and y values.
"""
function binXY(x::AbstractVector, y::AbstractVector)
    n = length(x)
    min_x, max_x = extrema(x)
    min_y, max_y = extrema(y)

    # mutual information binning:
    p2 = cor(x, y) ^ 2
    lx = (max_x - min_x) / std(x)
    ly = (max_y - min_y) / std(y)
    L = (n * p2) / (12 * (1 - p2)) * (lx^2 + ly^2)
    k = round(Int, 0.5 + 0.5 * sqrt(1 + 4 * sqrt(L)))

    Δ = (max_x - min_x) / k
    δ = Δ / 2
    bins = map(v -> BinXY(xkeys=[v-δ, v+δ]), min_x:Δ:max_x)

    for xy ∈ zip(x, y)
        idx = findfirst(b -> b.xkeys[1] <= xy[1] < b.xkeys[2], bins)
        if !isnothing(idx)
            push!(bins[idx].xvals, xy[1])
            push!(bins[idx].yvals, xy[2])
        end
    end

    bins
end


"""
    binstats(df::AbstractDataFrame, x::Symbol, y::Symbol)

Collect statistics (mean, SD, and SE) for x and y data in each bin.

# Arguments
- `df::AbstractDataFrame`: `DataFrame` containing the x and y values
- `x::Symbol`: x parameter/field in `df`
- `y::Symbol`: y parameter/field  in `df`

# Notes
- Every x and y groups will have their SD and SE determined
- Estimation of histogram bin count based on mutual information:
     Hacine-Gharbi, Deriche, M., Ravier, P., Harba, R., & Mohamadi, T. (2013).
     A new histogram-based estimation technique of entropy and mutual
     information using mean squared error minimization. Computers and
     Electrical Engineering, 39, 918-933.

# Returns
- `DataFrame`: mean, SD, and SE for x and y data in each bin.
"""
function binstats(df::AbstractDataFrame, x::Symbol, y::Symbol)
    xvals = getproperty(df, x)
    yvals = getproperty(df, y)
    bins = binXY(xvals, yvals)

    mean_x = map(b -> mean(b.xvals), bins)
    sd_x = map(b -> std(b.xvals), bins)
    se_x = map(b -> std(b.xvals) / sqrt(length(b.xvals)), bins)
    mean_y = map(b -> mean(b.yvals), bins)
    sd_y = map(b -> std(b.yvals), bins)
    se_y = map(b -> std(b.yvals) / sqrt(length(b.yvals)), bins)

    paramx, paramy = String(x), String(y)
    vals = [mean_x, sd_x, se_x, mean_y, sd_y, se_y]
    lbls = map(s -> ["mean_", "sd_", "se_"] .* s, [paramx, paramy])
    lbls = collect(Iterators.flatten(lbls))

    df = DataFrame([z[1] => z[2] for z ∈ zip(lbls, vals)])
    subset!(df, All() .=> ByRow(!isnan))
    df
end


"""
    @tictoc(expr) macro

Print and return elapsed time for a given code expression run.

# Arguments
- `expr`: expression to time

# Returns
- `Dates.Milliseconds`: elapsed time (ms)
"""
macro tictoc(expr)
    quote
        tic = now()
        $(esc(expr))
        toc = now()
        Δ = toc - tic
        e = Dates.canonicalize(Dates.CompoundPeriod(Δ))
        println("Completed in $e.")
        Δ
    end
end


"""
    gauss(N::Int, a, b)

Abscissas and weights for N-point Gauss-Legendre integration over an interval (a,b).

# Arguments
- `N::Int`: number of integration points
- `a`: lower integration limit
- `b`: upper integration limit

# Notes
- Code from: rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature#Julia

# Returns
- `Tuple{Float64, Float64}`: abscissas and weights
"""
function gauss(N::Int, a, b)
    λ, Q = eigen(SymTridiagonal(zeros(N), [n / sqrt(4n^2 - 1) for n ∈ 1:N-1]))
    @. (λ + 1) * (b - a) / 2 + a, [2Q[1, i]^2 for i ∈ 1:N] * (b - a) / 2
end


"""
    integrate(N::Int, lower, upper)

Numerical integration using N-point Gauss–Legendre quadrature method.

# Arguments
- `N::Int`: number of points for integration
- `lower`: lower boundary
- `upper`: upper boundary

# Notes
- Returns a function object for integration. This function object has the
  following signature:

      itg(func, args...; fields)

  where
     - `func` function to integrate that must return a named tuple with keys/fields
              listed in keyword `fields`
     - `args...`: arguments to pass into function `func`
     - `fields`: `Tuple` fields to return from the function `func`

   This function object returns a `Tuple` of integrated values, as listed
   in the argument `fields`.

# Returns
- Function object for integration (see Notes)
"""
function integrate(N::Int, lower, upper)
    x, w = gauss(N, lower, upper)   # abscissas and weights over (lower, upper) interval

    function itg(func, args...; fields)
        ret = map(a -> func(a, args...), x)
        map(f -> sum(map(r -> getproperty(r, f), ret) .* w), fields)
    end
end


"""
    mergepath(workdir::AbstractString, fname::AbstractString,
              defaultdir::AbstractString)

Merge the working folder path with file name path.

# Arguments
- `workdir::AbstractString`: working folder path (see Notes)
- `fname::AbstractString`: file name path
- `defaultdir::AbstractString`: default working folder path (see Notes)

# Notes
- If `workdir` is "@@", the default working folder path will be used.
  This default folder path is usually the location of the model input `JSON` file.
- Inserts '/' between the working folder path and file name path, if required.
- Replaces all back slashes to forward slashes.

# Returns
- `String`: merged path name
"""
function mergepath(workdir::AbstractString, fname::AbstractString,
                   defaultdir::AbstractString)
    function clean(arr, pos)
        # replace back slash with forward slash at specified character position
        ch = arr[pos:pos]   # get the character at the specified position
        if ch == "\\"
            arr = collect(arr)  # String is immutable, so convert to vector of Char
            arr[pos] = '/'      # replace the back slash with forward slash
            arr = join(arr)     # now convert back to String
            ch = "/"            # update the character
        end
        ch, arr
    end

    pth = fname
    if workdir != ""
        if workdir == "@@"
            # @@ code given, so use the provided default location
            workdir = defaultdir
        end

        ch1, workdir = clean(workdir, length(workdir))
        ch2, fname = clean(fname, 1)

        if ch1 == "/" && ch2 == "/"
            pth = wkd * fname[2:end]
        elseif (ch1 != "/" && ch2 == "/") || (ch1 == "/" && ch2 != "/")
            pth = "$(workdir)$(fname)"
        else
            pth = "$(workdir)/$(fname)"
        end
    end
    replace(pth, "\\" => "/")   # all back slashes to forward slashes
end


"""
    fraction_year(date::Date)

Convert `Date` into a fraction of its year.

# Arguments
- `date::Date`: date to be converted

# Returns
- `Float64`: equivalent decimal date of its year
"""
fraction_year(date::Date) = year(date) + (dayofyear(date) - 1) / daysinyear(date)


"""
    has_ext(fname::AbstractString, ext::AbstractVector{<:AbstractString})

Determine if a given filename has the provided file extension(s).

# Arguments
- `fname::AbstractString`: filename to check for the required extension
- `ext::AbstractVector{<:AbstractString}`: list fo extensions to find

# Returns
- `Bool`: `true` if the file name has the extension, else `false`
"""
function has_ext(fname::AbstractString, ext::AbstractVector{<:AbstractString})
    fnm = uppercase(fname)
    lst = uppercase.(ext)
    # find the characters after the last dot, but exclude the dot from text capture
    m = match(r"\.([^.]+$)", fnm)
    !isnothing(m) && (m[1] ∈ lst)
end


"""
    is_xl(fname::AbstractString)

Determine if a given filename has an `Excel` file extension.

# Arguments
- `fname::AbstractString`: filename to check

# Notes
- Files extensions "xlsx", "xlsm", or "xltx" are considered as `Excel`

# Returns
- `Bool`: `true` if the file name has an `Excel` file extension, else `false`
"""
is_xl(fname::AbstractString) = has_ext(fname, ["xlsx", "xlsm", "xltx"])


"""
    is_csv(fname::AbstractString)

Determine if a given filename has a CSV file extension.

# Arguments
- `fname::AbstractString`: filename to check

# Notes
- File extensions "csv", "txt", or "dat" are considered as `CSV`

# Returns
- `Bool`: `true` if the file name has a CSV file extension, else `false`
"""
is_csv(fname::AbstractString) = has_ext(fname, ["csv", "txt", "dat"])


"""
    round_df!(dfs::AbstractVector{<:AbstractDataFrame}, digits)

Round all float columns of a `DataFrame` to a specified number of digits.

# Arguments
- `dfs::AbstractVector{<:AbstractDataFrame}`: list of `DataFrame`
- `digits`: number of digits to round

# Returns
- nothing
"""
function round_df!(dfs::AbstractVector{<:AbstractDataFrame}, digits)
    is_float(col) = eltype(col) <: AbstractFloat

    transform!.(dfs,
                All() .=> c -> is_float(c) ? round.(c, digits=digits) : c,
                renamecols=false)
end
