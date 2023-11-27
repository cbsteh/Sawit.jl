"""
    add_axis!(axs::Vector{Axis},
              fig::Figure,
              r::Int, c::Int,
              xlabel::AbstractString, ylabel::AbstractString;
              maxcols::Int=4)

Add a new `Axis` to the chart `Figure`.

# Arguments
- `axs::Vector{Axis}`: array containing the `Figure` axes
- `fig::Figure`: chart `Figure`
- `r::Int`: row number in the chart layout
- `c::Int`: column number in the chart layout
- `xlabel::AbstractString`: label name of the x-axis
- `ylabel::AbstractString`: label name of the y-axis

# Keywords
- `maxcols::Int`: no. of columns per row (default=4)

# Notes
- Modifies the axes array `axs` by adding the new axis for the current `Figure`.

# Returns
- `Tuple{Int, Int}`: row and column number of the new axis
"""
function add_axis!(axs::Vector{Axis},
                   fig::Figure,
                   r::Int, c::Int,
                   xlabel::AbstractString, ylabel::AbstractString;
                   maxcols::Int=4)
    c = max(1, (c+1) % (maxcols+1))
    r = (c == 1) ? (r + 1) : r
    push!(axs, Axis(fig[r, c], xlabel=xlabel, ylabel=ylabel,
                    xgridvisible=false, ygridvisible=false))
    r, c
end


"""
    add_legend(ax::Axis, series, labels; position=:rb)

Add a `Legend` to current chart.

# Arguments
- `ax::Axis`: chart `Axis`
- `series`: chart `FigureAxisPlot`
- `labels`: series name/label

# Keywords
- `position`: position of legend withij axis (default is `:rb` for right bottom)

# Returns
- `Legend`
"""
function add_legend(ax::Axis, series, labels; position=:rb)
    axislegend(
        ax, series, labels,
        padding=(3,3,3,3),
        colgap=3,
        orientation=:horizontal,
        position=position,
        bgcolor=:transparent,
    )
end


"""
    create_figure(sz_in::Tuple{T, R};
                  fontpt::Int=20, dpi::Int=300) where {T<:Real, R<:Real}

Create a `Figure` for the chart(s), using default values.

# Arguments
- `sz_in::Tuple{T, R}`: width and height of figure in inches

# Keywords
- `fontpt::Int`: font size (points) (default: 20)
- `dpi::Int`: dots per inch for the figure size (default: 300)

# Returns
- `Tuple{Figure, Vector{Symbol}}`: the new chart `Figure` and color scheme
"""
function create_figure(sz_in::Tuple{T, R};
                       fontpt::Int=20, dpi::Int=300) where {T<:Real, R<:Real}
    sz_px = sz_in .* dpi
    fig = Figure(resolution=sz_px, font="Arial", fontsize=fontpt)

    # qualitative color scheme
    colors = [
        :grey0,
        :red,
        :green3,
        :orange,
        :magenta,
        :blue
    ]

    fig, colors
end


"""
    add_lines!(ax::Axis, x, ys, ylabels, colors; position=:rb)

Add multiple line series to a given chart axis.

# Arguments
- `ax::Axis`: chart `Axis`
- `x`: x values
- `ys`: y values for each series
- `ylabels`: y series labels
- `colors`: color to use for each line series

# Keywords
- `position`: position of legend withij axis (default is `:rb` for right bottom)

# Returns
- nothing
"""
function add_lines!(ax::Axis, x, ys, ylabels, colors; position=:rb)
    l = []
    foreach(eachindex(ys)) do i
        push!(l, lines!(ax, x, ys[i], color=colors[i], linewidth=3))
    end
    add_legend(ax, l, ylabels; position=position)

    yticks = niceticks(ys, 8)
    ax.yticks = yticks
    δ = yticks[2] - yticks[1]
    if position == :rt || position == :lt
        ylims!(ax, yticks[1]-0.5δ, yticks[end]+δ)
    elseif position == :rb || position == :lb
        ylims!(ax, yticks[1]-δ, yticks[end]+0.5δ)
    end
end


"""
    plot_annual(df::AbstractDataFrame, mi::Input;
                saveplot::Bool=false, plotfname::AbstractString="annual.png")

Plot annual output results.

# Arguments
- `df::AbstractDataFrame`: annual model output results
- `mi::Input`: model input

# Keywords
- `saveplot:Bool`: `true` to save the plots, else `false` for no saves (default)
- `plotfname::AbstractString`: plot file name (default: "annual.png")

# Returns
- plot figure
"""
function plot_annual(df::AbstractDataFrame, mi::Input;
                     saveplot::Bool=false, plotfname::AbstractString="annual.png")
    chtsz = 2 # inches
    nrows, ncols = 5, 5
    rowsz, colsz = nrows * chtsz, ncols * chtsz

    fig, colors = create_figure((colsz, rowsz); fontpt=30)
    axs = Axis[]

    yap = df.yap    # x values
    nlayers = mi.nlayers
    row = col = 0
    mxc = 5

    # yield:
    row, col = add_axis!(axs, fig, row, col, "YAP", "yield"; maxcols=mxc)
    lines!(axs[end], yap, df.yield, color=colors[1], linewidth=3)

    # LAI:
    row, col = add_axis!(axs, fig, row, col, "YAP", "LAI"; maxcols=mxc)
    lines!(axs[end], yap, df.lai, color=colors[1], linewidth=3)

    # pinnae, rachis, and roots weights:
    row, col = add_axis!(axs, fig, row, col, "YAP", "weight"; maxcols=mxc)
    ys = [df.parts_pinnae_weight, df.parts_rachis_weight, df.parts_roots_weight]
    ylabels = ["pinnae", "rachis", "roots"]
    add_lines!(axs[end], yap, ys, ylabels, colors; position=:lt)

    # trunk weight:
    row, col = add_axis!(axs, fig, row, col, "YAP", "trunk weight"; maxcols=mxc)
    lines!(axs[end], yap, df.parts_trunk_weight, color=colors[1], linewidth=3)

    # flowers and bunches weight:
    row, col = add_axis!(axs, fig, row, col, "YAP", "weight"; maxcols=mxc)
    ys = [df.parts_maleflo_weight, df.parts_femaflo_weight, df.parts_bunches_weight]
    ylabels = ["male", "female", "bunches"]
    add_lines!(axs[end], yap, ys, ylabels, colors; position=:lt)

    # assimilates and maintenance, growth, and generative respiration:
    row, col = add_axis!(axs, fig, row, col, "YAP", "weight"; maxcols=mxc)
    ys = [df.assimilates, df.assim4maint, df.assim4growth, df.assim4gen]
    ylabels = ["asm", "mnt", "gro", "gen"]
    add_lines!(axs[end], yap, ys, ylabels, colors; position=:lt)

    # VDM and TDM:
    row, col = add_axis!(axs, fig, row, col, "YAP", "weight"; maxcols=mxc)
    ys = [df.vdmgro, df.tdmgro]
    ylabels = ["VDM", "TDM"]
    add_lines!(axs[end], yap, ys, ylabels, colors; position=:lt)

    # trunk height:
    row, col = add_axis!(axs, fig, row, col, "YAP", "trunk height"; maxcols=mxc)
    lines!(axs[end], yap, df.trunkhgt, color=colors[1], linewidth=3)

    # rooting depth:
    row, col = add_axis!(axs, fig, row, col, "YAP", "roots depth"; maxcols=mxc)
    lines!(axs[end], yap, df.roots_depth, color=colors[1], linewidth=3)

    # gross and net rainfall:
    row, col = add_axis!(axs, fig, row, col, "YAP", "rain"; maxcols=mxc)
    ys = [df.rain, df.netrain]
    ylabels = ["gross", "net"]
    add_lines!(axs[end], yap, ys, ylabels, colors; position=:lt)

    # root zone VWC:
    row, col = add_axis!(axs, fig, row, col, "YAP", "roots VWC"; maxcols=mxc)
    lines!(axs[end], yap, df.roots_vwc, color=colors[1], linewidth=3)

    # root zone WC:
    row, col = add_axis!(axs, fig, row, col, "YAP", "roots WC"; maxcols=mxc)
    lines!(axs[end], yap, df.roots_wc, color=colors[1], linewidth=3)

    # deep percolation:
    row, col = add_axis!(axs, fig, row, col, "YAP", "percolation"; maxcols=mxc)
    lines!(axs[end], yap, df.dp, color=colors[1], linewidth=3)

    # ET:
    row, col = add_axis!(axs, fig, row, col, "YAP", "ET"; maxcols=mxc)
    ys = [df.aet_soil, df.aet_crop, df.aet_soil .+ df.aet_crop]
    ylabels = ["E", "T", "ET"]
    add_lines!(axs[end], yap, ys, ylabels, colors; position=:lt)

    # water deficit:
    row, col = add_axis!(axs, fig, row, col, "YAP", "water deficit"; maxcols=mxc)
    lines!(axs[end], yap, df.waterdeficit, color=colors[1], linewidth=3)

    # influx:
    row, col = add_axis!(axs, fig, row, col, "YAP", "influx"; maxcols=mxc)
    ys = map(i -> getproperty(df, Symbol("layers$(i)_influx")), 1:nlayers)
    ylabels = string.([1:nlayers...])
    add_lines!(axs[end], yap, ys, ylabels, colors; position=:lt)

    # outflux:
    row, col = add_axis!(axs, fig, row, col, "YAP", "outflux"; maxcols=mxc)
    ys = map(i -> getproperty(df, Symbol("layers$(i)_outflux")), 1:nlayers)
    add_lines!(axs[end], yap, ys, ylabels, colors; position=:lt)

    # netflux:
    row, col = add_axis!(axs, fig, row, col, "YAP", "netflux"; maxcols=mxc)
    ys = map(i -> getproperty(df, Symbol("layers$(i)_netflux")), 1:nlayers)
    add_lines!(axs[end], yap, ys, ylabels, colors; position=:lt)

    # T:
    row, col = add_axis!(axs, fig, row, col, "YAP", "actual T"; maxcols=mxc)
    ys = map(i -> getproperty(df, Symbol("layers$(i)_t")), 1:nlayers)
    add_lines!(axs[end], yap, ys, ylabels, colors; position=:lt)

    # min. and max. air temperatures:
    row, col = add_axis!(axs, fig, row, col, "YAP", "air temp."; maxcols=mxc)
    ys = [df.tmin, df.tmax]
    ylabels = [L"T_{min}", L"T_{max}"]
    add_lines!(axs[end], yap, ys, ylabels, colors; position=:lt)

    # solar irradiance:
    row, col = add_axis!(axs, fig, row, col, "YAP", "radiation"; maxcols=mxc)
    ys = [(df.totrad ./ 1000), (df.drrad ./ 1000), (df.dfrad ./ 1000)]
    ylabels = [L"I_{total}", L"I_{dr}", L"I_{df}"]
    add_lines!(axs[end], yap, ys, ylabels, colors; position=:lt)

    # wind speed:
    row, col = add_axis!(axs, fig, row, col, "YAP", "wind speed"; maxcols=mxc)
    lines!(axs[end], yap, df.wind, color=colors[1], linewidth=3)

    row, col = add_axis!(axs, fig, row, col, "YAP", "CGR"; maxcols=mxc)
    lines!(axs[end], yap, df.cgr, color=colors[1], linewidth=3)

    row, col = add_axis!(axs, fig, row, col, "YAP", "BI"; maxcols=mxc)
    lines!(axs[end], yap, df.bi, color=colors[1], linewidth=3)

    row, col = add_axis!(axs, fig, row, col, "YAP", "NAR"; maxcols=mxc)
    lines!(axs[end], yap, df.nar, color=colors[1], linewidth=3)

    xticks = niceticks(yap)
    foreach(ax -> ax.xticks = xticks, axs)

    saveplot && save(plotfname, fig)
    display(fig)
end


"""
    plot_daily(df::AbstractDataFrame, mi::Input;
               saveplot::Bool=false, plotfname::AbstractString="daily.png")

Plot daily output results.

# Arguments
- `df::AbstractDataFrame`: daily model output results
- `mi::Input`: model input

# Keywords
- `saveplot:Bool`: `true` to save the plots, else `false` for no saves (default)
- `plotfname::AbstractString`: plot file name (default: "daily.png")

# Returns
- nothing
"""
function plot_daily(df::AbstractDataFrame, mi::Input;
                    saveplot::Bool=false, plotfname::AbstractString="daily.png")
    chtsz = 2 # inches
    nrows, ncols = 5, 4
    rowsz, colsz = nrows * chtsz, ncols * chtsz

    fig, colors = create_figure((colsz, rowsz); fontpt=30)
    axs = Axis[]

    dates = fraction_year.(df.date)
    nlayers = mi.nlayers
    row = col = 0

    row, col = add_axis!(axs, fig, row, col, "year", "water stress")
    lines!(axs[end], dates, df.stress_t, color=colors[1], linewidth=3)

    row, col = add_axis!(axs, fig, row, col, "year", "water deficit")
    lines!(axs[end], dates, df.et_crop .- df.aet_crop, color=colors[1], linewidth=3)

    row, col = add_axis!(axs, fig, row, col, "year", "percolation")
    lines!(axs[end], dates, df.dp, color=colors[1], linewidth=3)

    row, col = add_axis!(axs, fig, row, col, "year", "roots VWC")
    lines!(axs[end], dates, df.roots_vwc, color=colors[1], linewidth=3)

    row, col = add_axis!(axs, fig, row, col, "year", "roots WC")
    lines!(axs[end], dates, df.roots_wc, color=colors[1], linewidth=3)

    row, col = add_axis!(axs, fig, row, col, "year", "actual ET")
    ys = [df.aet_soil, df.aet_crop]
    ylabels = ["soil", "crop"]
    add_lines!(axs[end], dates, ys, ylabels, colors; position=:lt)

    row, col = add_axis!(axs, fig, row, col, "year", "actual total ET")
    lines!(axs[end], dates, df.aet_soil .+ df.aet_crop, color=colors[1], linewidth=3)

    row, col = add_axis!(axs, fig, row, col, "year", "H")
    ys = [df.h_soil, df.h_crop]
    add_lines!(axs[end], dates, ys, ylabels, colors; position=:lt)

    row, col = add_axis!(axs, fig, row, col, "year", "total H")
    lines!(axs[end], dates, df.h_soil .+ df.h_crop, color=colors[1], linewidth=3)

    row, col = add_axis!(axs, fig, row, col, "year", "layer T")
    l = lines!(axs[end], dates, df.layers1_t, color=colors[1], linewidth=3)
    add_legend(axs[end], [l], ["1"])

    row, col = add_axis!(axs, fig, row, col, "year", "layer T")
    ys = map(i -> getproperty(df, Symbol("layers$(i)_t")), 2:nlayers)
    ylabels = string.([2:nlayers...])
    add_lines!(axs[end], dates, ys, ylabels, colors; position=:lt)

    row, col = add_axis!(axs, fig, row, col, "year", "layer VWC")
    l = lines!(axs[end], dates, df.layers1_vwc, color=colors[1], linewidth=3)
    add_legend(axs[end], [l], ["1"])

    row, col = add_axis!(axs, fig, row, col, "year", "layer VWC")
    ys = map(i -> getproperty(df, Symbol("layers$(i)_vwc")), 2:nlayers)
    add_lines!(axs[end], dates, ys, ylabels, colors; position=:lt)

    row, col = add_axis!(axs, fig, row, col, "year", "layer WC")
    l = lines!(axs[end], dates, df.layers1_wc, color=colors[1], linewidth=3)
    add_legend(axs[end], [l], ["1"])

    row, col = add_axis!(axs, fig, row, col, "year", "layer WC")
    lyr = []
    ys = map(i -> getproperty(df, Symbol("layers$(i)_wc")), 2:nlayers)
    ylabels = string.([2:nlayers...])
    add_lines!(axs[end], dates, ys, ylabels, colors; position=:lt)

    # assimilates and maintenance, growth, and generative respiration:
    row, col = add_axis!(axs, fig, row, col, "year", "weight")
    ys = [df.assimilates, df.assim4maint, df.assim4growth, df.assim4gen]
    ylabels = ["asm", "mnt", "gro", "gen"]
    add_lines!(axs[end], dates, ys, ylabels, colors; position=:lt)

    # VDM and TDM:
    row, col = add_axis!(axs, fig, row, col, "year", "weight")
    ys = [df.vdmgro, df.tdmgro]
    ylabels = ["VDM", "TDM"]
    add_lines!(axs[end], dates, ys, ylabels, colors; position=:lt)

    # gross and net rainfall:
    row, col = add_axis!(axs, fig, row, col, "year", "rain")
    ys = [df.rain, df.netrain]
    ylabels = ["gross", "net"]
    add_lines!(axs[end], dates, ys, ylabels, colors; position=:lt)

    # min. and max. air temperatures:
    row, col = add_axis!(axs, fig, row, col, "year", "air temp.")
    ys = [df.tmin, df.tmax]
    ylabels = [L"T_{min}", L"T_{max}"]
    add_lines!(axs[end], dates, ys, ylabels, colors; position=:lt)

    # solar irradiance:
    row, col = add_axis!(axs, fig, row, col, "year", "radiation")
    ys = [df.totrad, df.drrad, df.dfrad]
    ylabels = [L"I_{total}", L"I_{dr}", L"I_{df}"]
    add_lines!(axs[end], dates, ys, ylabels, colors; position=:lt)

    # wind speed:
    row, col = add_axis!(axs, fig, row, col, "year", "wind speed")
    lines!(axs[end], dates, df.wind, color=colors[1], linewidth=3)

    xticks = niceticks(dates)
    foreach(ax -> ax.xticks = xticks, axs)

    saveplot && save(plotfname, fig)
    display(fig)
end


"""
    plot_hourly(df::AbstractDataFrame, mo::Output,
                saveplot::Bool=false, plotfname::AbstractString="hourly.png")

Plot hourly output results for a given day.

# Arguments
- `df::AbstractDataFrame`: annual model output results
- `mo::Output`: model output
- `mi::Input`: model input

# Keywords
- `saveplot:Bool`: `true` to save the plots, else `false` for no saves (default)
- `plotfname::AbstractString`: plot file name (default: "hourly.png")

# Returns
- plot figure
"""
function plot_hourly(df::AbstractDataFrame, mo::Output;
                     saveplot::Bool=false, plotfname::AbstractString="hourly.png")
    chtsz = 2 # inches
    nrows, ncols = 5, 4
    rowsz, colsz = nrows * chtsz, ncols * chtsz

    fig, colors = create_figure((colsz, rowsz); fontpt=30)
    axs = Axis[]

    th = df.solarhour
    # some plots only for between sunrise and sunset hour
    tsr, tss = mo.op.sunrise+1.5, mo.op.sunset-1.5
    th_f = th[tsr .< th .< tss]
    df_f = filter(:solarhour => (t -> tsr < t < tss), df)

    row = col = 0

    # air and canopy temperatures:
    row, col = add_axis!(axs, fig, row, col, "hour", "temp.")
    ys = [df.air_temp, df.canopy_temp]
    ylabels = ["air", "canopy"]
    add_lines!(axs[end], th, ys, ylabels, colors; position=:lt)

    # air vapor pressures:
    row, col = add_axis!(axs, fig, row, col, "hour", "vap. pressure")
    ys = [df.svp, df.vp, df.vpd, df.vpd0]
    ylabels = ["SVP", "VP", "VPD", L"VPD_0"]
    add_lines!(axs[end], th, ys, ylabels, colors; position=:lt)

    # wind speeds:
    row, col = add_axis!(axs, fig, row, col, "hour", "wind speed")
    ys = [df.wind, df.utreehgt]
    ylabels = ["above", "canopy top"]
    add_lines!(axs[end], th, ys, ylabels, colors; position=:lt)

    # solar radiation components:
    row, col = add_axis!(axs, fig, row, col, "hour", "radiation")
    ys = [df.totrad, df.drrad, df.dfrad, df.netrad]
    ylabels = ["total", "direct", "diffuse", "net"]
    add_lines!(axs[end], th, ys, ylabels, colors; position=:lt)

    # LAI components:
    row, col = add_axis!(axs, fig, row, col, "hour", "LAI")
    ys = [df.lsl, df.lsh]
    ylabels = ["sunlit", "shaded"]
    add_lines!(axs[end], th, ys, ylabels, colors; position=:lt)

    # PAR irradiance components within and outside canopies:
    row, col = add_axis!(axs, fig, row, col, "hour", "PAR")
    ys = [df.par_outdr, df.par_outdf, df.par_indr, df.par_inscatter, df.par_indf]
    ylabels = [L"out_{dr}", L"out_{df}", L"in_{dr}", L"in_{scatter}", L"in_{df}"]
    add_lines!(axs[end], th, ys, ylabels, colors; position=:lt)

    # PAR absorption components:
    row, col = add_axis!(axs, fig, row, col, "hour", "PAR abs.")
    ys = [df.par_abssunlit, df.par_absshaded]
    ylabels = ["sunlit", "shaded"]
    add_lines!(axs[end], th, ys, ylabels, colors; position=:lt)

    # CO2 levels:
    row, col = add_axis!(axs, fig, row, col, "hour", "CO2 level")
    ys = [df.co2_ambient, df.co2_internal]
    ylabels = ["ambient", "internal"]
    add_lines!(axs[end], th, ys, ylabels, colors; position=:lt)

    # leaf assimilation of CO2 components:
    row, col = add_axis!(axs, fig, row, col, "hour", "CO2 assim.")
    vq = df_f.leafassim_vqsl .+ df_f.leafassim_vqsh
    ys = [df_f.leafassim_vc, df_f.leafassim_vqsl, df_f.leafassim_vqsh, vq]
    ylabels = [L"v_c ", L"v_{q,sl}", L"v_{q,sh}", L"v_{q,total}"]
    add_lines!(axs[end], th_f, ys, ylabels, colors; position=:lt)

    # sunlit and shaded leaf and canopy assimilation of CO2:
    row, col = add_axis!(axs, fig, row, col, "hour", "CO2 assim.")
    ys = [df.leafassim_sunlit, df.leafassim_shaded, df.canopyassim]
    ylabels = ["sunlit", "shaded", "canopy"]
    add_lines!(axs[end], th, ys, ylabels, colors; position=:lt)

    # available energy components:
    row, col = add_axis!(axs, fig, row, col, "hour", "avail. energy")
    ys = [df.A_total, df.A_soil, df.A_crop, df.A_net, df.A_g]
    ylabels = ["total", "soil", "crop", L"R_{net}", "G"]
    add_lines!(axs[end], th, ys, ylabels, colors; position=:lt)

    # stomatal stresses:
    row, col = add_axis!(axs, fig, row, col, "hour", "stresses")
    ys = [df.stresses_water, df.stresses_vpd, df.stresses_par]
    ylabels = ["water", "VPD", "PAR"]
    add_lines!(axs[end], th, ys, ylabels, colors; position=:lt)

    # aerodynamic resistances:
    row, col = add_axis!(axs, fig, row, col, "hour", "resistance")
    ys = [df.R_rsa, df.R_rca, df.R_raa]
    ylabels = [L"r_s^a", L"r_c^a", L"r_a^a"]
    add_lines!(axs[end], th, ys, ylabels, colors; position=:lt)

    # stomatal and canopy resistances:
    row, col = add_axis!(axs, fig, row, col, "hour", "resistance")
    ys = [df_f.R_rst, df_f.R_rcs]
    ylabels = [L"r_{st}", L"r_c^s"]
    add_lines!(axs[end], th_f, ys, ylabels, colors; position=:lt)

    # soil surface resistance:
    row, col = add_axis!(axs, fig, row, col, "hour", "soil resistance")
    lines!(axs[end], th, df.R_rss, color=colors[1], linewidth=3)

    # latent heat flux components:
    row, col = add_axis!(axs, fig, row, col, "year", "ET")
    ys = [df.et_total, df.et_soil, df.et_crop]
    ylabels = ["total", "soil", "crop"]
    add_lines!(axs[end], th, ys, ylabels, colors; position=:lt)

    # sensible heat flux components:
    row, col = add_axis!(axs, fig, row, col, "year", "H")
    ys = [df.h_total, df.h_soil, df.h_crop]
    ylabels = ["total", "soil", "crop"]
    add_lines!(axs[end], th, ys, ylabels, colors; position=:lt)

    # canopy extinction coefficient components:
    row, col = add_axis!(axs, fig, row, col, "hour", "ext. coef.")
    ys = [df_f.kdr, df_f.kdf]
    ylabels = [L"k_{dr}", L"k_{df}"]
    add_lines!(axs[end], th_f, ys, ylabels, colors; position=:lt)

    # solar height (elevation):
    row, col = add_axis!(axs, fig, row, col, "hour", "solar hgt.")
    lines!(axs[end], th_f, rad2deg.(df_f.solarhgt), color=colors[1], linewidth=3)

    # canopy clump factor:
    row, col = add_axis!(axs, fig, row, col, "hour", "clump factor")
    lines!(axs[end], th_f, df_f.clump, color=colors[1], linewidth=3)

    saveplot && save(plotfname, fig)
    display(fig)
    fig
end


"""
    explore(xys::AbstractVector{T};
            title::AbstractString="",
            saveplot::Bool=false,
            plotfname::AbstractString="hourly.png") where {T}

Plot a series of x against y for exploratory analysis.

# Arguments
- `xys::AbstractVector{T}`: array of (x, y) values/columns and labels to plot

# Keywords
- `title::AbstractString`: plot title (default: empty string/no title)
- `saveplot:Bool`: `true` to save the plots, else `false` for no saves (default)
- `plotfname::AbstractString`: plot file name (default: "explore.png")

# Returns
- plot figure
"""
function explore(xys::AbstractVector{T};
                 title::AbstractString="",
                 saveplot::Bool=false,
                 plotfname::AbstractString="explore.png") where {T}
    sz = length(xys)
    ncols = min(4, sz)
    n, r = divrem(sz, 4)
    nrows = (r > 0) ? n + 1 : n

    fig, colors = create_figure((2 * ncols, 2 * nrows))
    axs = Axis[]
    Label(fig[0, 1:ncols], title)   # chart title

    row = col = 0
    for xy ∈ xys
        xlbl, ylbl = xy[3][1], xy[3][2]
        row, col = add_axis!(axs, fig, row, col, xlbl, ylbl)
        scatter!(axs[end], xy[1], xy[2], color=colors[1], markersize=10)
    end

    saveplot && save(plotfname, fig)
    display(fig)
    fig
end
