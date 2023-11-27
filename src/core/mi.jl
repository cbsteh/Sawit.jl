include("mi.afgen.jl")
include("mi.cc.jl")


"""
    Prop

Soil layer physical properties.

# Fields
- `thick::Float64`: layer thickness (m)
- `clay::Float64`: clay (%)
- `sand::Float64`: sand (%)
- `vwc::Float64`: volumetric soil water content VWC (m3/m3)
"""
@with_kw mutable struct Prop
    thick::Float64 = 0.0
    clay::Float64 = 0.0
    sand::Float64 = 0.0
    vwc::Float64 = 0.0
end


"""
    load_wthr(fname::AbstractString)

Read in (load) weather data from file (`Excel` or`CSV` format).

# Arguments
- `fname::AbstractString`: weather file name

# Notes
- The following fields must be defined (case-sensitive) in the weather
  data file (field order is unimportant):

     year, month, day, tmin, tmax, wind, rain

- If weather file is an `Excel` workbook, the weather data must be stored
  in a worksheet named "data".
- Error is displayed if file name has an unrecognized file extension.

# Returns
- `DataFrame`: weather data (with additional column for `Date`)
"""
function load_wthr(fname::AbstractString)
    isfile(fname) || error("ERROR: File \"$fname\" not found.")

    kys = [:year, :month, :day, :tmin, :tmax, :wind, :rain]
    types = [Int, Int, Int, Float64, Float64, Float64, Float64]
    fields = Dict([zip(kys, types)]...)

    if is_csv(fname)
        # multi-threading for small to medium CSV files gives inconsistent read
        #    speeds, so turn off multi-threading by setting `ntasks` to 1 thread
        df = CSV.File(fname, comment="#", ntasks=1, types=fields, validate=false) |>
             DataFrame
    elseif is_xl(fname)
        df = Dataset(XLSX.readtable(fname, "data")...)
        for (k, t) ∈ fields
            df[!, k] .= t.(df[!, k])
        end
    else
        error("Given file type not applicable.")
    end

    transform!(df, kys[1:3] => ByRow(Date) => :date)     # create a Date column
end


"""
    Input

Model input parameters.

# Fields
- `jsondir::String`: model input JSON file path
- `datadir::String`: data folder path
- `boxcarfname::String`: Box car file name
- `outputfname::String`: model output file name
- `wthrfname::String`: weather file name
- `seed::Int`: seed for the random number generator
- `numdays::Int`: number of simulation days (days)
- `lat::Float64`: site latitude (radians)
- `methgt::Float64`: weather station height (m)
- `date::Date`: date (i.e., simulation start date)
- `wthr::DataFrame`: `DataFrame` to hold the weather data
- `treeage::Int`: age of tree (days)
- `plantdens::Int`: planting density (palms/ha)
- `thinplantdens::Int`: thinning planting density (palms/ha) (-1 for no thinning)
- `thinage::Int`: age of tree for thinning (days)
- `femaleprob::Float64`: probability getting female flowers (0 - 1)
- `pinnae_wgt::Float64`: pinnae weights (kg DM/palm)
- `rachis_wgt::Float64`: rachis weight (kg DM/palm)
- `trunk_wgt::Float64`: trunk weight (kg DM/palm)
- `roots_wgt::Float64`: roots weight (kg DM/palm)
- `maleflo_wgt::Float64`: male flowers weight (kg DM/palm)
- `femaflo_wgt::Float64`: female flowers weight (kg DM/palm)
- `bunches_wgt::Float64`: bunches weight (kg DM/palm)
- `sla_table::Afgen`: table of tree age vs. specific leaf area (m2 leaf/kg DM leaf)
- `pinnae_n_table::Afgen`: table of tree age vs. pinnae N content (fraction)
- `pinnae_m_table::Afgen`: table of tree age vs. pinnae mineral content (fraction)
- `rachis_n_table::Afgen`: table of tree age vs. rachis N content (fraction)
- `rachis_m_table::Afgen`: table of tree age vs. rachis mineral content (fraction)
- `roots_n_table::Afgen`: table of tree age vs. roots N content (fraction)
- `roots_m_table::Afgen`: table of tree age vs. roots mineral content (fraction)
- `trunk_n_table::Afgen`: table of tree age vs. trunk N content (fraction)
- `trunk_m_table::Afgen`: table of tree age vs. trunk mineral content (fraction)
- `nintervals::Int`: no. of subintervals per day for daily water integration
- `rootdepth::Float64`: rooting depth (m)
- `maxrate_rootdepth::Float64`: max. rate of increase in rooting depth (m/day)
- `has_watertable::Bool`: `true` for water table, else `false` for none
- `nlayers::Int`: no. of soil layers
- `ccfn::CCfn`: climate change functions

# Notes
- `maxrate_rootdepth` stored as m/day but entered as mm/day in input file
"""
@with_kw mutable struct Input
    jsondir::String = ""
    datadir::String = ""
    boxcarfname::String = ""
    outputfname::String = ""
    wthrfname::String = ""

    seed::Int = -1
    numdays::Int = 0

    lat::Float64 = 0.0
    methgt::Float64 = 0.0
    date::Date = Date(1970, 1, 1)
    wthr::DataFrame = DataFrame()

    treeage::Int = 365
    plantdens::Int = 148
    trunkhgt::Float64 = 0.0
    thinplantdens::Int = -1
    thinage::Int = 3650
    femaleprob::Float64 = 0.5

    pinnae_wgt::Float64 = 0.4
    rachis_wgt::Float64 = 0.7
    trunk_wgt::Float64 = 0.1
    roots_wgt::Float64 = 0.2
    maleflo_wgt::Float64 = 0.0
    femaflo_wgt::Float64 = 0.0
    bunches_wgt::Float64 = 0.0

    sla_table::Afgen = Afgen()
    pinnae_n_table::Afgen = Afgen()
    pinnae_m_table::Afgen = Afgen()
    rachis_n_table::Afgen = Afgen()
    rachis_m_table::Afgen = Afgen()
    roots_n_table::Afgen = Afgen()
    roots_m_table::Afgen = Afgen()
    trunk_n_table::Afgen = Afgen()
    trunk_m_table::Afgen = Afgen()

    nintervals::Int = 0
    rootdepth::Float64 = 0.0
    maxrate_rootdepth::Float64 = 0.0
    has_watertable::Bool = false
    nlayers::Int = 0
    layers::Vector{Prop} = []

    ccfn::CCfn = CCfn()
end


"""
    read_model_input(jsonfname::AbstractString)

Read the model input file into the `Input` object.

# Arguments
- `jsonfname::AbstractString`: model input JSON file name

# Returns
- `Input`: model input
"""
function read_model_input(jsonfname::AbstractString)
    isfile(jsonfname) || error("File \"$jsonfname\" not found.")

    dct = JSON.parsefile(jsonfname; dicttype=Dict)
    jsondir = dirname(jsonfname)   # reference to model input file path
    datadir = dct["data folder"]

    boxcarfname = dct["boxcar filename"]
    if !isempty(boxcarfname)
        boxcarfname = mergepath(datadir, boxcarfname, jsondir)
    end

    outputfname = mergepath(datadir, dct["output filename"], jsondir)
    wthrfname = mergepath(datadir, dct["weather filename"], jsondir)

    nlayers = dct["no. of soil layers"]
    layers = Vector{Prop}(undef, nlayers)
    jl = dct["layers"]
    for i ∈ 1:nlayers
        thick = jl[i]["thick"]
        clay = jl[i]["clay"]
        sand = jl[i]["sand"]
        vwc = jl[i]["vwc"]
        layers[i] = Prop(thick, clay, sand, vwc)
    end

    Input(
        jsondir = jsondir,
        datadir = datadir,
        boxcarfname = boxcarfname,
        outputfname = outputfname,
        wthrfname = wthrfname,
        seed = dct["random generator seed"],
        numdays = dct["no. of simulation days"],
        lat = deg2rad(dct["site latitude"]),   # always work in radians
        methgt = dct["weather station hgt"],
        date = Date(dct["year"], dct["month"], dct["day"]),
        wthr = load_wthr(wthrfname),  # read in entire weather data
        treeage = dct["tree age"],
        plantdens = dct["planting density"],
        trunkhgt = dct["trunk height"],
        thinplantdens = dct["thinning planting density"],
        thinage = dct["thinning tree age"],
        femaleprob = dct["female prob"],
        pinnae_wgt = dct["pinnae wgt"],
        rachis_wgt = dct["rachis wgt"],
        trunk_wgt = dct["trunk wgt"],
        roots_wgt = dct["roots wgt"],
        maleflo_wgt = dct["male flowers wgt"],
        femaflo_wgt = dct["female flowers wgt"],
        bunches_wgt = dct["bunches wgt"],
        sla_table = Afgen(dct["sla"]),
        pinnae_n_table = Afgen(dct["pinnae n"]),
        pinnae_m_table = Afgen(dct["pinnae m"]),
        rachis_n_table = Afgen(dct["rachis n"]),
        rachis_m_table = Afgen(dct["rachis m"]),
        roots_n_table = Afgen(dct["roots n"]),
        roots_m_table = Afgen(dct["roots m"]),
        trunk_n_table = Afgen(dct["trunk n"]),
        trunk_m_table = Afgen(dct["trunk m"]),
        nintervals = dct["no. of intervals"],
        rootdepth = dct["rooting depth"],
        maxrate_rootdepth = dct["max. rate rooting depth"] / 1000,  # store as m/day
        has_watertable = dct["has watertable"],
        nlayers = nlayers,
        layers = layers,
        ccfn = CCfn(dct)
    )
end
