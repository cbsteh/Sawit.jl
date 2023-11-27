"""
    save_results(fname::AbstractString,
                 daly::AbstractDataFrame, annl::AbstractDataFrame)

Save the daily and annual simulation results into `CSV` (default) or `Excel` files.

# Arguments
- `fname::AbstractString`: file name to save the daily and annual model output (see Notes)
- `daly::AbstractDataFrame`: daily simulation results
- `annl::AbstractDataFrame`: annual simulation results

# Notes
- If the file extension in `fname` is `XLSX`, the model output results will be saved
  as an `Excel` workbook with two worksheets, named "daily" and "annual" for daily and
  annual results, respectively.
- If not `Excel`, then it is assumed, by default, that the model output results will
  be saved as `CSV`-delimited text files.
- If saved to `CSV` files, the words "-daily" and "-annual" will be appended to the
  file name `fname` for daily and annual model output results, respectively.

# Returns
- nothing
"""
function save_results(fname::AbstractString,
                      daly::AbstractDataFrame, annl::AbstractDataFrame)
    if is_xl(fname)
        XLSX.openxlsx(fname, mode="w") do xf
            while XLSX.sheetcount(xf) < 2
                XLSX.addsheet!(xf)
            end
            sht1 = xf[1]
            XLSX.rename!(sht1, "daily")
            XLSX.writetable!(sht1, collect(eachcol(daly)), names(daly))
            sht2 = xf[2]
            XLSX.rename!(sht2, "annual")
            XLSX.writetable!(sht2, collect(eachcol(annl)), names(annl))
        end
    else
        m = findlast(".", fname)
        pre = fname[1:m[1]-1]
        pos = fname[m[1]:end]
        dalyfname = "$(pre)-daily$(pos)"
        CSV.write(dalyfname, daly)
        annlfname = "$(pre)-annual$(pos)"
        CSV.write(annlfname, annl)
    end
end


"""
    summarize_annual(df::AbstractDataFrame, mi::Input)

Summarize daily model output into annual model output.

# Arguments
- `df::AbstractDataFrame`: daily model output results as `DataFrame`
- `mi::Input`: model input

# Returns
- `DataFrame`: annual model output
"""
function summarize_annual(df::AbstractDataFrame, mi::Input)
    # fields to sum
    sum_fields = [
        :yield,
        :totrad, :drrad, :dfrad,
        :rain, :netrain,
        :assimilates, :assim4maint, :assim4growth, :assim4gen,
        :vdmgro, :tdmgro,
        :et_crop, :et_soil, :aet_crop, :aet_soil,
        :dp
    ]
    # fields to get their final values at end of every YAP period
    lastval_fields = [
        :lai,
        :parts_pinnae_weight, :parts_rachis_weight, :frdwgt,
        :parts_trunk_weight, :parts_roots_weight,
        :parts_maleflo_weight, :parts_femaflo_weight, :parts_bunches_weight,
        :trunkhgt, :roots_depth,
        :roots_vwc, :roots_wc
    ]
    # fields to average
    mean_fields = [
        :tmin, :tmax, :wind
    ]
    # soil layer field names (as strings to make it easier to append later)
    layer_fields = [
        "t", "influx", "outflux", "netflux",
        "vwc", "wc"
    ]

    nlayers = mi.nlayers
    f = ["layers$(i)_$(f)" for i ∈ 1:nlayers, f ∈ layer_fields] |> vec .|> Symbol
    m = nlayers * 4
    n = nlayers * 6
    append!(sum_fields, f[1:m])
    append!(lastval_fields, f[m+1:n])

    ann = combine(
        groupby(df, :yap),
        sum_fields .=> sum,
        lastval_fields .=> f -> f[end],
        mean_fields .=> Statistics.mean,
        renamecols=false
    )

    # crop growth rate (t/ha/yr)
    ann.cgr = (ann.vdmgro .+ ann.yield) .* mi. plantdens / 1000
    # bunch index
    ann.bi = ann.yield ./ (ann.yield .+ ann.vdmgro)
    # net assimilation rate (g/m2 leaf/week)
    ann.nar = (ann.cgr ./ ann.lai) * (100 / 52)
    # water deficit (mm/year)
    ann.waterdeficit = ann.et_crop .- ann.aet_crop

    ann
end


"""
    fluxes2mm!(sl::SoilLayer)

Convert water fluxes from unit m to mm.

# Arguments
- `sl::SoilLayer`: Water fluxes in soil layer to which to convert to unit mm

# Returns
- nothing
"""
function fluxes2mm!(sl::SoilLayer)
    sl.k *= 1000
    sl.t *= 1000
    sl.e *= 1000
    sl.influx *= 1000
    sl.outflux *= 1000
    sl.netflux *= 1000
end


"""
    init(jsonfname::AbstractString, dct::Union{AbstractDict, Nothing}=nothing)

Read the model input file, make custom changes to the model input,
then initialize the model. Use before running the model.

# Arguments
- `jsonfname::AbstractString`: model input JSON file name
- `dct::Union{AbstractDict, Nothing}`: model scenario (custom changes to model input)
                                       (default: nothing, meaning no scenario)

# Notes
- Model inputs that have been read from the model input file can still be
  modified by specifiying the changes in the `dct` argument. For instance,
  to change the planting density to 148 palms/ha, write:

     dct = Dict(:plantdens => 148)
     mo, mi = init(jsonfname, dct)

  where `:plantdens` is the `Symbol` for the `MI` field for planting
  density. Note: only field names (case sensitive) for `MI` struct argument
  are recognized. Any unrecognized field names will be ignored.

# Returns
- `NamedTuple{(:mo, :mi), Tuple{Output, Input}}`: model output and model input
"""
function init(jsonfname::AbstractString, dct::Union{AbstractDict, Nothing}=nothing)
    # read in the model input from JSON file:
    mi = read_model_input(jsonfname)
    if !isnothing(dct)
        input_fields = fieldnames(Input)
        kys = [k for k ∈ keys(dct) if k ∈ input_fields]
        foreach(k -> setproperty!(mi, k, dct[k]), kys)
    end

    # initialize the model:
    if isempty(mi.boxcarfname)
        bc = BoxCar()
    else
        dct = JSON.parsefile(mi.boxcarfname; dicttype=Dict)
        bc = BoxCar(
            maleflo = dct["maleflo"],
            femaflo = dct["femaflo"],
            bunches = dct["bunches"]
        )
    end

    mo = Output(op=OP(date=mi.date), bc=bc)
    init_crop!(mo, mi)
    init_soilwater!(mo, mi)

    (; mo, mi)
end


"""
    runsim!(mo::Output, mi::Input)

Simulation run.

# Arguments
- `mo::Output`: model output results
- `mi::Input`: model input

# Returns
- `NamedTuple{(:daly, :annl, :mi), Tuple{DataFrame, DataFrame, Input}}`:
     daily and annual model output results, and the model input
"""
function runsim!(mo::Output, mi::Input)
    (mi.seed > 0) && Random.seed!(mi.seed)
    numdays = mi.numdays
    mo.op.date -= (numdays > 0) ? Day(1) : Day(0)
    results = Vector{OP}(undef, numdays)     # store only `OP`

    @inbounds for nday ∈ 1:numdays
        mo.op.date += Day(1)
        update_wthr!(mo, mi)       # order of execution is important
        update_photosyn!(mo, mi)
        update_soilwater!(mo, mi)
        update_crop!(mo, mi)
        results[nday] = clone(mo.op)  # deep copy model results
    end

    foreach(r -> fluxes2mm!.(r.layers), results)
    daly = tabulate(results)
    annl = summarize_annual(daly, mi)
    round_df!([daly, annl], 3)
    save_results(mi.outputfname, daly, annl)

    (; daly, annl, mo, mi)
end


"""
    runsimhour!(jsonfname::AbstractString;
                saveoutput::Bool=true,
                outfname::AbstractString="hour.csv",
                plot::Bool=false,
                saveplot::Bool=false,
                plotfname::AbstractString="hourly.png")

Plot hourly output results for a given day.

# Arguments
- `jsonfname::AbstractString`: model input JSON file name

# Keywords
- `saveoutput::Bool`: `true` to save results to file (default), else `false` for no save
- `outfname::AbstractString`: name of output file (default: "hour.csv")
- `plot::Bool`: `true` to plot results, else `false` for no plot (default)
- `saveplot:Bool`: `true` to save validation charts, else `false` for no saves (default)
- `plotfname::AbstractString`: plot file name (default: "hourly.png")

# Returns
- plot figure
"""
function runsimhour!(jsonfname::AbstractString;
                     saveoutput::Bool=true,
                     outfname::AbstractString="hour.csv",
                     plot::Bool=false,
                     saveplot::Bool=false,
                     plotfname::AbstractString="hourly.png")
    mo, mi = init(jsonfname)
    (mi.seed > 0) && Random.seed!(mi.seed)

    update_wthr!(mo, mi)
    spas = soil_plant_atmos!(mo, mi)
    results = Vector{OPHour}(undef, 24)
    foreach(th -> results[th+1] = spas(th).ophr, 0:23)
    hrly = tabulate(results)

    saveoutput && CSV.write(outfname, hrly)
    plot && plot_hourly(hrly, mo; saveplot=saveplot, plotfname=plotfname)
    (; hrly, mo, mi)
end


"""
    savestate(model_fname::AbstractString, boxcar_fname::AbstractString,
              mo::Output, mi::Input)

Save the model object and Box car state to file. Files are saved in JSON format.

# Arguments
- `model_fname::AbstractString`: file name to save the model object state
- `boxcar_fname::AbstractString`: file name to save the Box car object state
- `mo::Output`: model output results
- `mi::Input`: model input

# Notes
- Example of use:

    jsonfname = "data/demo/input.json"
    res = start(jsonfname)
    savestate("data/demo/model.json", "data/demo/bc.json", res.mo, res.mi)

# Returns
- nothing
"""
function savestate(model_fname::AbstractString, boxcar_fname::AbstractString,
                   mo::Output, mi::Input)
    function to_dict(af::Afgen)
        OrderedDict(map(xy -> string(xy[1]) => xy[2], zip(af.x, af.y)))
    end

    @unpack op, bc = mo

    # save Box car state:
    bc_dct = OrderedDict(map(f -> f => getproperty(bc, f), fieldnames(BoxCar)))
    open(boxcar_fname, "w") do f
        JSON.print(f, bc_dct, 4)
    end

    parts = op.parts
    fields = [:thick, :vwc, :clay, :sand]
    layers = map(l -> OrderedDict(map(f -> f => getproperty(l, f), fields)), op.layers)

    midct = OrderedDict(
        "random generator seed" => mi.seed,
        "no. of simulation days" => mi.numdays,

        "data folder" => "",
        "boxcar filename" => mi.boxcarfname,
        "weather filename" => mi.wthrfname,
        "output filename" => mi.outputfname,

        "site latitude" => rad2deg(mi.lat),
        "weather station hgt" => mi.methgt,

        "year" => Year(op.date).value,
        "month" => Month(op.date).value,
        "day" => Day(op.date).value,

        "co2_path" => mi.ccfn.co2_path,
        "ta_path" => mi.ccfn.Δta_path,
        "wind_path" => mi.ccfn.Δwind_path,
        "rain_path" => mi.ccfn.Δrain_path,

        "tree age" => op.treeage,
        "planting density" => op.plantdens,
        "trunk height" => op.trunkhgt,
        "thinning planting density" => mi.thinplantdens,
        "thinning tree age" => mi.thinage,
        "female prob" => mi.femaleprob,
        "pinnae wgt" => parts.pinnae.weight,
        "rachis wgt" => parts.rachis.weight,
        "trunk wgt" => parts.trunk.weight,
        "roots wgt" => parts.roots.weight,
        "male flowers wgt" => parts.maleflo.weight,
        "female flowers wgt" => parts.femaflo.weight,
        "bunches wgt" => parts.bunches.weight,

        "no. of intervals" => mi.nintervals,
        "rooting depth" => op.roots.depth,
        "max. rate rooting depth" => mi.maxrate_rootdepth * 1000,
        "has watertable" => mi.has_watertable,
        "no. of soil layers" => mi.nlayers,
        "layers" => layers,

        "sla" => to_dict(mi.sla_table),
        "pinnae n" => to_dict(mi.pinnae_n_table),
        "pinnae m" => to_dict(mi.pinnae_m_table),
        "rachis n" => to_dict(mi.rachis_n_table),
        "rachis m" => to_dict(mi.rachis_m_table),
        "roots n" => to_dict(mi.roots_n_table),
        "roots m" => to_dict(mi.roots_m_table),
        "trunk n" => to_dict(mi.trunk_n_table),
        "trunk m" => to_dict(mi.trunk_m_table)
    )

    open(model_fname, "w") do f
        JSON.print(f, midct, 4)
    end
end


"""
    start(jsonfname::AbstractString, scenario::AbstractVector{<:AbstractDict})

Read the model input file, initialize the model with modeling scenarios,
then start simulation runs for each modeling scenario.

# Arguments
- `jsonfname::AbstractString`: model input JSON file name
- `scenario::AbstractVector{<:AbstractDict}`: list of modeling scenarios

# Returns
- `Vector{NamedTuple{(:daly, :annl. :mi), Tuple{DataFrame, DataFrame, Input}}`:
     list of daily and annual model output results and the model input
"""
function start(jsonfname::AbstractString, scenario::AbstractVector{<:AbstractDict})
    mo_mi = map(dct -> init(jsonfname, dct), scenario)
    Folds.map(m -> runsim!(m.mo, m.mi), mo_mi)
end


"""
    start(jsonfname::AbstractString)

Read the model input file, initialize the model, then start simulation runs.

# Arguments
- `jsonfname::AbstractString`: model input JSON file name

# Notes
- This method initializes the model with the values read from the model input
  file, then immediately runs the model, and returns the results.
- In most cases, this method is the entry point of the model. However, to make
  changes to the model input, after their values have been read from the
  mode input file, use method:

     init(jsonfname::AbstractString,
          dct::AbstractDict{Symbol, <:Any}=Dict{Symbol, Int}())

  to specify the changes to the model input, then use method:

     run(mo::Output, mi::Input)

  to run the model.

# Returns
- `NamedTuple{(:daly, :annl, :mi), Tuple{DataFrame, DataFrame, Input}`:
     daily and annual model output results, and the model input
"""
function start(jsonfname::AbstractString)
    mo, mi = init(jsonfname)
    daly, annl = runsim!(mo, mi)
    (; daly, annl, mo, mi)
end
