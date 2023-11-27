"""
    ParDct

Glossary of model output parameters and their respective units.

# Fields
- `par::Dict{String, Vector{String}}`: `Dict` holding parameter names and units
"""
@with_kw struct ParDct
    par::Dict{String, Vector{String}} = Dict(
        "yap" => ["YAP", ""],
        "yield" => ["FFB", "kg DM palm^{-1} yr^{-1}"],
        "vdmgro" => ["VDM", "kg DM palm^{-1} yr^{-1}"],
        "tdmgro" => ["TDM", "kg DM palm^{-1} yr^{-1}"],
        "parts_pinnae_weight" => ["pinnae\\,weight", "kg DM palm^{-1}"],
        "parts_rachis_weight" => ["rachis\\,weight", "kg DM palm^{-1}"],
        "frdwgt" => ["frond\\,weight", "kg DM palm^{-1}"],
        "parts_trunk_weight" => ["trunk\\,weight", "kg DM palm^{-1}"],
        "lai" => ["LAI", "m^2 m^{-2}"],
        "trunkhgt" => ["trunk\\,height", "m"],
    )
end


"""
    getpar(pd::ParDct, par)

Retrieve name of parameter and text for chart axis labeling.

# Arguments
- `pd::ParDct`: glossary of all parameters and their units
- `par`: name of parameter

# Notes
- Returns the name of given parameter and without unit if given parameter
  is not in the glossary.
- The parameter label is for Latex labeling

# Returns
- `NamedTuple{(:param, :label), Tuple{String, String}}`: name of parameter
      and parameter text for chart labeling
"""
function getpar(pd::ParDct, par)
    param, unit = get(pd.par, par, [par, ""])
    label = param * (isempty(unit) ? "" : " ($unit)")
    unit = replace(unit, " " => "\\,")    # add latex space
    label = replace(label, " " => "\\,")
    (; param, unit, label)
end


"""
    add_inset_axis!(fig::Figure, r::Int, c::Int, fontpt::Int=20)

Add an inset axis to an exisitng `Figure`.

# Arguments
- `fig::Figure`: chart `Figure`
- `r::Int`: row number in the chart layout
- `c::Int`: column number in the chart layout
- `fontpt::Int`: font size (points) for the y labels (default: 20)

# Returns
- `Axis`: new inset axis
"""
function add_inset_axis!(fig::Figure, r::Int, c::Int, fontpt::Int=20)
    ax = Axis(
        fig[r, c],
        width=Relative(0.3), height=Relative(0.3),
        halign=0.15, valign=0.95,
        xlabelsize=fontpt, yticklabelsize=fontpt, ytickalign=1,
    )
    hidexdecorations!(ax; label=false)
    hideydecorations!(ax; ticklabels=false, ticks=false)
    ax
end


"""
    plot_validation(simdf::AbstractDataFrame,
                    pairdf::AbstractDataFrame,
                    gofdf::AbstractDataFrame,
                    title::AbstractString,
                    saveplot::Bool,
                    plotfname::AbstractString)

Plot agreement between observations and estimations (simulation results).

# Arguments:
- `simdf::AbstractDataFrame` : all simulated values
- `pairdf::AbstractDataFramw`: observed and simulated paired values, arranged in such
                               a way that for each paramater, its observed column
                               is adjacent to its simulated column
- `gofdf::AbstractDataFrame`: Goodness-of-fit tests results
- `title::AbstractString`: chart title (default: empty string)
- `saveplot::Bool`: set `true` to save the plots, else `false` for no saves
- `plotfname::AbstractString`: file name to save the plots

# Notes
- INTERNAL USE.
- Should not be directly used, rather called by `validation` function,
  because `unidf` and `gofdf` must be arranged in a specific order/manner
  (see `validate`).

# Returns
- `Figure`: chart figure
"""
function plot_validation(simdf::AbstractDataFrame,
                         pairdf::AbstractDataFrame,
                         gofdf::AbstractDataFrame,
                         title::AbstractString,
                         saveplot::Bool,
                         plotfname::AbstractString)
    N = size(gofdf)[1]
    ncols = min(4, N * 2)    # no. of charts per row we want
    nrows = Int(ceil(N * 2 / ncols))   # no. of charts to plot
    chtsz = 2    # chart size in inches
    rowsz, colsz = nrows * chtsz, ncols * chtsz

    fig, colors = create_figure((colsz, rowsz); fontpt=30)
    axs = Axis[]
    Label(fig[0, 1:ncols], title)   # chart title

    pairfields = names(pairdf)
    x = pairfields[1]        # key name (e.g., "YAP")
    obsx = pairdf[!, x]      # key values for observed
    simx = simdf[!, x]       # key values for simulated

    pardct = ParDct()

    goffields = map(p -> getpar(pardct, p).param, names(gofdf))    # names of GoF indexes
    goflbls = join(goffields[2:end], ", ") * "\n"   # text box for all GoF names

    row = col = 0
    r = 0   # track current col. no. in `pairdf` `DataFrame`

    for n ∈ 2:2:(length(pairfields)-1)
        r += 1
        param = gofdf[r, 1]     # current parameter being drawn
        simy = simdf[!, pairfields[n]]   # all simulated values
        o = pairdf[!, pairfields[n]]     # observed and simulated pair
        s = pairdf[!, pairfields[n+1]]

        # observed values may have missing values:
        keep = map(!ismissing, o)
        ox = collect(Float64, obsx[keep])
        o = collect(Float64, o[keep])
        s = collect(Float64, s[keep])

        parstr = getpar(pardct, param)
        parx = getpar(pardct, x)
        xparam = parx.param
        yparam = parstr.param
        xunit = isempty(parx.unit) ? "" : "($(parx.unit))"
        yunit = isempty(parstr.unit) ? "" : "($(parstr.unit))"

        # 1. scatter and line charts
        row, col = add_axis!(axs, fig, row, col,
                             L"\mathrm{%$xparam}\,%$xunit",
                             L"\mathrm{%$yparam}\,%$yunit"; maxcols=ncols)
        l = lines!(axs[end], simx, simy, color=colors[1], linewidth=3)
        scat = scatter!(axs[end], ox, o, color=colors[1], markersize=15)
        (n == 2) && add_legend(axs[end], [l, scat], ["sim.", "obs."])

        allx = unique(vcat(simx, ox))
        axs[end].xticks = niceticks(allx)

        # 2. obs vs. sim with bisector (1:1 line / line of agreement)
        paramdf = gofdf[r, 2:end]
        paramvals = round.((i -> paramdf[i][1]).(eachindex(paramdf)), digits=2)
        txt = (n == 2 ? goflbls : "") * join(string.(paramvals), ", ")

        row, col = add_axis!(axs, fig, row, col,
                             L"\mathrm{sim.\,(%$yparam)}",
                             L"\mathrm{obs.\,(%$yparam)}"; maxcols=ncols)
        scatter!(axs[end], s, o, color=colors[1], markersize=15)
        dummy = scatter!(axs[end], [s[1]], [o[1]], color=:transparent, markersize=1)
        ablines!(axs[end], 0, 1, color=colors[1], linestyle=:dash, linewidth=3)
        add_legend(axs[end], [dummy], [txt])

        # draw a 1:1 chart, with equal x- and y-axis limits
        xy = [s; o]
        mn = minimum(xy)
        mx = maximum(xy)
        Δ = (mx - mn) * 0.1     # leave a 10% margin
        mn -= Δ
        mx += Δ
        limits!(axs[end], (mn, mx), (mn, mx))

        # 3. box plot
        err = s .- o
        inset_ax = add_inset_axis!(fig, row, col)
        xs = fill(1, length(err))
        violin!(inset_ax, xs, err; color=:lightgray, strokewidth=1)
        boxplot!(inset_ax, xs, err;
                 width=0.4, strokewidth=1, color=(colors[1], 0.3),
                 outliercolor=colors[1], markersize=10)
        xlims!(inset_ax, (0.5, 1.5))
        inset_ax.yticks = niceticks(err, 4)
        inset_ax.xlabel = "sim. - obs."
        inset_ax.backgroundcolor = :transparent
    end

    saveplot && save(plotfname, fig)
    display(fig)
    fig
end


"""
    validate(obsfname::AbstractString,
             simdf::AbstractDataFrame,
             mi::Input,
             x::Symbol=:yap,
             gof::AbstractVector{Function}=[GoF.NMBE, GoF.NMAE, GoF.KGE],
             plot::Bool=true,
             title::AbstractString="",
             saveplot::Bool=false,
             plotfname::AbstractString="validation.png")

Validate estimations (simulation results) against observations using specified
Goodness-of-fit indexes.

# Arguments
- `obsfname::AbstractString`: file name for observation (measured) data
- `simdf::AbstractDataFrame`: model estimation/simulation results
- `mi::Input`: model input
- `x::Symbol`: key (default is :yap)

# Keywords
- `gof::::AbstractVector{<:Function}`: list of goodness-of-fit `GoF` functions to use
                                       (default: [NMBE, NMAE, KGE])
- `plot:Bool`: `true` to plot validation charts (default), else `false` for no plots
- `title::AbstractString`: chart title (default: empty string)
- `saveplot:Bool`: `true` to save the plots, else `false` for no saves (default)
- `plotfname::AbstractString`: plot file name (default: "validation.png")

# Notes
- Validation charts can only be saved if they are plotted. In other words,
  argument `saveplot` is only effective when argument `plot` is `true`. If `plot`
  is `false`, whatever value `saveplot` has will be ignored.

# Returns
- `NamedTuple{(:unidf, :gofdf), Tuple{DataFrame, DataFrame}`:
     1) `DataFrame` with adjacent (paired) columns of observed and simulated values
        for each parameter
     2) `DataFrame` containing the goodness-of-fit tests
"""
function validate(obsfname::AbstractString,
                  simdf::AbstractDataFrame,
                  mi::Input,
                  x::Symbol=:yap;
                  gof::AbstractVector{<:Function}=[GoF.NMBE, GoF.NMAE, GoF.KGE],
                  plot::Bool=true,
                  title::AbstractString="",
                  saveplot::Bool=false,
                  plotfname::AbstractString="validation.png")
    # read in the observed data and convert them into a `DataFrame`:
    obsfname = mergepath(mi.datadir, obsfname, mi.jsondir)
    obsdf = CSV.File(obsfname; comment="#", missingstring=["-999"]) |> DataFrame

    # create a unified (combined) `DataFrame` for both observed and simulated data:
    pairdf = innerjoin(obsdf, simdf, on=x, makeunique=true)

    ofields = names(obsdf)      # obs
    sfields = names(simdf)      # sim
    cfields = ofields ∩ sfields   # select only common parameters in both obs and sim
    strx = String(x)
    deleteat!(cfields, findall(==(strx), cfields))  # remove key from GoF comparisons

    lst = (f -> (f, "$(f)_1")).(cfields)  # pair obs vs sim columns for each parameter
    renlst = [strx]     # add key here
    append!(renlst, lst...)
    select!(pairdf, renlst)  # trim columns to only those with matching obs vs sim
    rename!(f -> replace(f, r"_1$" => "_sim"), pairdf)   # rename all sim parameters

    goflst = []
    for field ∈ cfields
        o = pairdf[!, field]
        s = pairdf[!, "$(field)_sim"]
        o, s = GoF.filter_data(o, s)   # exclude missing data or those with 0 value
        fit = (g -> g(o, s)).(gof)     # do every GoF test on current parameter
        nt = ((g, f) -> Symbol(g)=>f).(gof, fit)
        append!(goflst, ((param=field, nt...), ))
    end

    # some GoF indexes return a `NamedTuple`, so parse its members
    #    into their individual columns:
    gofdf = DataFrame(goflst)
    cols = [c for c in names(gofdf) if eltype(getproperty(gofdf, c)) <: NamedTuple]
    transform!(gofdf, cols .=> AsTable; renamecols=false)
    select!(gofdf, Not(cols))

    plot && plot_validation(simdf, pairdf, gofdf, title, saveplot, plotfname)
    (; pairdf, gofdf)
end


"""
    function plot_gof(param::AbstractVector{<:AbstractString},
                     series::AbstractVector{<:AbstractString},
                     x::AbstractVector{<:AbstractVector},
                     y::AbstractVector{<:AbstractVector},
                     z::AbstractVector{<:AbstractVector};
                     zmin, zmax,
                     xopt, yopt,
                     xlabel::AbstractString,
                     ylabel::AbstractString,
                     zlabel::AbstractString,
                     flabel=1.0,
                     fdist=1.0,
                     mesh_on::Bool=false,
                     link_axes::Bool=false,
                     saveplot::Bool=false,
                     plotfname::AbstractString="gof.svg",
                     ncolors::Int=5,
                     set_seed::Bool=true)

Plot three goodness-of-fit (GoF) indexes as a 2D scatter chart, with
non-overlapping label placement (EXPERIMENTAL).

# Arguments
- `param::AbstractVector{<:AbstractString}`: list of model parameter names
- `series::AbstractVector{<:AbstractString}`: list of series
- `x::AbstractVector{<:AbstractVector}`: list of first GoF index values
- `y::AbstractVector{<:AbstractVector}`: list of second GoF index values
- `z::AbstractVector{<:AbstractVector}`: list of third GoF index values

# Keywords
- `zmin`: lowest value (worst fit) for the third GoF index, `z`
- `zmax`: highest value (best fit) for the third GoF index, `z`
- `xopt`: highest value (best fit) for the first GoF index, `x`
- `yopt`: highest value (best fit) for the second GoF index, `y`
- `xlabel::AbstractString`: name of the first GoF index
- `ylabel::AbstractString`: name of the second GoF index
- `zlabel::AbstractString`: name of the third GoF index
- `flabel`: fraction by which to increase facet size (default: 1.0)
- `fdist`: inflate the max distance between facet and its anchor point
           before a line will be drawn between them (default: 1.0)
- `mesh_on::Bool`: set `true` to draw a mesh over the plot, else `false` for no
                   mesh overlay. Used for testing if facet size is appropriate
                   (default: `false`)
- `link_axes::Bool`: set `true` to link all x and y axes (default: `false`)
- `saveplot:Bool`: `true` to save the plots, else `false` for no saves (default)
- `plotfname::AbstractString`: plot file name (default: "gof.png")
- `ncolors::Int`: no. of intervals in the color scheme for the `z` GoF index (default: 5)
- `set_seed::Bool`: set `true` to always randomize label placement attempts
                    (default: `true`), else `false` for deterministic attempts

# Notes
- A chart is drawn for each model parameter, and each chart comprises
  one or more series.
- The third GoF index, `z`, will be plotted as color values.
- `jet` color scheme is used. It is reversed so that low values are in red
  and high values in blue.

# Returns
- `Figure`: chart figure
"""
function plot_gof(param::AbstractVector{<:AbstractString},
                  series::AbstractVector{<:AbstractString},
                  x::AbstractVector{<:AbstractVector},
                  y::AbstractVector{<:AbstractVector},
                  z::AbstractVector{<:AbstractVector};
                  zmin, zmax,
                  xopt, yopt,
                  xlabel::AbstractString,
                  ylabel::AbstractString,
                  zlabel::AbstractString,
                  flabel=1.0,
                  fdist=1.0,
                  mesh_on::Bool=false,
                  link_axes::Bool=false,
                  saveplot::Bool=false,
                  plotfname::AbstractString="gof.png",
                  ncolors::Int=5,
                  set_seed::Bool=true)
    N = length(param)
    ncols = min(3, N)
    nrows = max(1, Int(ceil(N / ncols)))
    chtsz = 2    # chart size in inches
    rowsz, colsz = nrows * chtsz, ncols * chtsz

    dpi = 300
    res = (chtsz * dpi, chtsz * dpi)
    fig, _ = create_figure((colsz, rowsz); fontpt=30, dpi=dpi)
    axs = Axis[]
    row = col = 0
    fontpt = 20

    # adjust as necessary:
    # fdist=0.5
    # fontpt=20
    # flabel=2.5

    cmap = cgrad(:jet, ncolors, rev=true, categorical=true)
    clims = (zmin, zmax)

    lastx = last_row_per_col(ncols, nrows, length(param))

    pardct = ParDct()
    p = map(p -> getpar(pardct, p).param, param)

    if link_axes
        xi_opt = vcat(collect(Iterators.flatten(x)), xopt)
        yi_opt = vcat(collect(Iterators.flatten(y)), yopt)
        mesh = create_mesh(extrema(xi_opt), extrema(yi_opt), series;
                           resolution=res, flabel=flabel, fontpt=fontpt)
    end

    for i ∈ eachindex(p)
        xi, yi = x[i], y[i]

        xlbl = i in lastx ? xlabel : ""
        ylbl = ((i-1) % ncols == 0) ? ylabel : ""

        row, col = add_axis!(axs, fig, row, col, xlbl, ylbl; maxcols=ncols)
        axs[end].title = L"\mathrm{%$(p[i])}"
        scatter!(axs[end], xi, yi, marker=:circle, markersize=20,
                 color=z[i], colormap=cmap, colorrange=clims)

        if !link_axes
            xi_opt = vcat(xi, xopt)
            yi_opt = vcat(yi, yopt)
            mesh = create_mesh(extrema(xi_opt), extrema(yi_opt), series;
                               resolution=res, flabel=flabel, fontpt=fontpt)
        end

        set_seed && Random.seed!(Dates.value(now()))

        fill_mesh!(mesh, xi, yi, series; niter=100, max_order=5)
        mesh_on && draw_mesh!(axs[end], mesh)
        draw_labels!(axs[end], mesh, fdist, fontpt)

        scatter!(axs[end], [xopt], [yopt], marker='⊗', markersize=25, color=:black)

        # add dummy points to have a fixed scale:
        scatter!(axs[end], [0.0, 0.6], [-0.25, 0.6], markersize=1, color=:white)
    end

    link_axes && linkaxes!(axs...)

    ival = (zmax - zmin) / ncolors
    Colorbar(fig[nrows+1, 1:max(1, ncols-1)], colormap=cmap, limits=clims,
             ticks=zmin:ival:zmax, label=zlabel,
             labelsize=20, ticklabelsize=20, vertical=false)

    saveplot && save(plotfname, fig)

    display(fig)
    fig
end


"""
    last_row_per_col(ncols::Int, nrows::Int, n::Int)

Return the chart numbers where these charts are in the last row of each column.
(INTERNAL USE).

# Arguments
- `ncols::Int`: number of columns in grid
- `nrows::Int`: number of rows in grid
- `n::Int`: number of charts

# Notes
- Charts in the last row are not necessarily always in the last row of each column,
  because some columns in the last row may not have any charts, e.g., plotting
  10 charts in a 3 by 4 grid layout, where the last two columns in the last
  row will be empty. In this case, this function will return [8, 9, 10] to
  indicate that the charts in row 4 for column  1, row 2 for column 2, and row 2
  for column 3 are the charts in the last row.

# Returns
- `Vector{Int}`: the last row number for each column
"""
function last_row_per_col(ncols::Int, nrows::Int, n::Int)
    lastrow = nrows
    while (lastrow - 1) * ncols >= n
        lastrow -= 1
    end
    N = lastrow * ncols
    d = N - n
    k = max(1, (lastrow - 1) * ncols + 1 - d)
    [k:n...]
end


"""
    plot_pairs(obs::AbstractDataFrame, sim::AbstractDataFrame;
               gof::AbstractVector{<:Function}=[GoF.NMBE, GoF.NMAE, GoF.KGE],
               title::AbstractString="",
               saveplot::Bool=false,
               plotfname::AbstractString="pairs.png")

Plot agreement between pairs of observations and estimations for various parameters.

# Arguments
- `obs::AbstractDataFrame` : observation values
- `sim::AbstractDataFramw`: simulated values

# Keywords
- `gof::::AbstractVector{<:Function}`: list of goodness-of-fit `GoF` functions to use
                                       (default: [NMBE, NMAE, KGE])
- `title::AbstractString`: chart title (default: empty string)
- `saveplot::Bool`: set `true` to save the plots, else `false` for no saves
- `plotfname::AbstractString`: file name to save the plots (default: "pairs.png")

# Notes
- `obs` and `sim` must have the same `DataFrame` column names.

# Returns
- `Figure`: chart figure
"""
function plot_pairs(obs::AbstractDataFrame, sim::AbstractDataFrame;
                    gof::AbstractVector{<:Function}=[GoF.NMBE, GoF.NMAE, GoF.KGE],
                    title::AbstractString="",
                    saveplot::Bool=false,
                    plotfname::AbstractString="pairs.png")
    par_obs = names(obs)
    N = length(par_obs)

    ncols = min(3, N)
    nrows = Int(ceil(N / ncols))
    chtsz = 2    # chart size in inches
    rowsz, colsz = nrows * chtsz, ncols * chtsz

    lastx = last_row_per_col(ncols, nrows, N)
    pardct = ParDct()
    goflbls = join(map(g -> "$g", gof), ", ") * "\n"

    fig, colors = create_figure((colsz, rowsz); fontpt=30)
    axs = Axis[]
    Label(fig[0, 1:ncols], title)   # chart title

    row = col = 0

    for i ∈ eachindex(par_obs)
        po = par_obs[i]
        o = obs[:, po]
        s = sim[:, po]

        idxs = findall(v -> !ismissing(v), o)
        o = o[idxs]
        s = s[idxs]

        fit = round.(map(g -> g(o, s), gof); digits=2)
        txt = (i == 1 ? goflbls : "") * join(string.(fit), ", ")

        xlabel = i in lastx ? "sim." : ""
        ylabel = ((i-1) % ncols == 0) ? "obs." : ""
        row, col = add_axis!(axs, fig, row, col, xlabel, ylabel; maxcols=ncols)
        scatter!(axs[end], s, o, color=colors[1], markersize=12)
        dummy = scatter!(axs[end], [s[1]], [o[1]], color=:transparent, markersize=1)
        ablines!(axs[end], 0, 1, color=colors[1], linestyle=:dash, linewidth=3)
        add_legend(axs[end], [dummy], [txt])

        par = getpar(pardct, po)
        axs[end].title = L"\mathrm{%$(par.param)}\,(%$(par.unit))"

        err = s .- o
        inset_ax = add_inset_axis!(fig, row, col)
        xs = fill(1, length(err))
        violin!(inset_ax, xs, err; color=:lightgray, strokewidth=1)
        boxplot!(inset_ax, xs, err;
                 width=0.4, strokewidth=1, color=(colors[1], 0.3),
                 outliercolor=colors[1], markersize=8)
        xlims!(inset_ax, (0.5, 1.5))
        inset_ax.yticks = niceticks(err, 5)
        inset_ax.xlabel = "sim. - obs."
        inset_ax.backgroundcolor = :transparent

        xy = [s; o]
        mn = minimum(xy)
        mx = maximum(xy)
        Δ = (mx - mn) * 0.08     # leave a little margin to avoid clipping
        mn -= Δ
        mx += Δ
        limits!(axs[end], (mn, mx), (mn, mx))
    end

    saveplot && save(plotfname, fig)

    display(fig)
    fig
end
