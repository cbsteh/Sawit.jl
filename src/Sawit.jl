module Sawit

using CairoMakie
using CSV
using DataFrames
using Dates
using Folds
using Interpolations
using JSON
using LinearAlgebra
using OrderedCollections
using Parameters
using Random
using Statistics
using XLSX

include("utils/utils.jl")
include("utils/GoF.jl")

include("core/mi.jl")
include("core/mo.jl")
include("core/meteo.jl")
include("core/layer.jl")
include("core/soilwater.jl")
include("core/crop.jl")
include("core/energybal.jl")
include("core/photosyn.jl")

include("plots/nicescale.jl")
include("plots/plotsim.jl")
include("plots/cmesh2d.jl")
include("plots/algolabels.jl")
include("plots/validation.jl")

include("runsim.jl")

# utils.jl
export @tictoc
export BinXY, BinXYZ
export binXY, binXYZ, binstats, integrate, mergepath, fraction_year,
       is_xl, is_csv, has_ext, gauss, round_df!

# GoF.jl
using .GoF
export GoF

# mi.jl
export Prop, Input, Afgen, CCfn
export load_wthr, read_model_input

# mo.jl
export BoxCar, StomatalStresses, AvailEnergy, Resistances, HeatFluxes, SoilLayer,
       RootsLayer, OP, Output, OPHour, Properties, Parts
export tabulate, clone, vegparts, vegparts!, allparts, allparts!

# meteo.jl
export solar_declination, solar_constant, et_solar_radiation, saturated_vapor_pressure,
       vapor_pressure, vapor_pressure_deficit, relative_humidity, svp_fn,
       twilight_hours, day_et_solar_radiation, day_solar_radiation, solar_position,
       solar_radiation, net_solar_radiation, air_temperature, svp_slope_fn,
       wind_speed, read_wthr, est_hour_wthr, update_wthr!

# layer.jl
export to_dg, est_airentry, est_psd, est_critical, to_wc, to_vwc, est_swc,
       ksat, heads_k, new_layer, roots_water!, new_roots_layer

# soilwater.jl
export init_soilwater!, update_soilwater!

# crop.jl
export dm_weights, lai_maximum, vdm_maximum, init_crop!, update_crop!

# energybal.jl
export update_energybal!

# photosyn.jl
export update_photosyn!

# runsim.jl
export save_results, summarize_annual, init, runsim, runsimhour!, savestate, start

# nicescale
export niceticks

# plotsim.jl
export plot_annual, plot_daily, plot_hourly, explore,
       create_figure, add_axis!, add_legend, add_lines!

# cmesh2d.jl
export Pt, Facet, CMesh2D
export distance, mean_dist, min_dist, mpos, mindex, get_neighbors, centroid

# algolabels.jl
export is_available, has_label, reset_place!, save_places, load_places!,
       swap_place!, find_places!, labelsize, segments_intersect,
       create_mesh, fill_mesh!, draw_mesh!, draw_labels!, viz

# validation.jl
export ParDct
export add_inset_axis!, getpar, validate, plot_gof, plot_pairs

end   # module
