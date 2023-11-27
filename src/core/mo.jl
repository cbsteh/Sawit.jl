include("mo.parts.jl")
include("mo.hour.jl")


"""
    BoxCar

Boxcar method to represent the growth cycle of male and female flowers and bunches.

# Fields
- `maleflo::Vector{Float64}`: male flowers life cycle
- `femaflo::Vector{Float64}`: female flowers life cycle
- `bunches::Vector{Float64}`: bunches life cycle

# Notes
- Expensive to be deep copied.
"""
@with_kw struct BoxCar
    maleflo::Vector{Float64} = zeros(210)  # 210 days full cycle
    femaflo::Vector{Float64} = zeros(210)  # same as male flowers
    bunches::Vector{Float64} = zeros(150)  # 150 days full cycle
end


"""
    SoilLayer

Soil layer physical and water properties.

# Fields
- `thick::Float64`: layer thickness (m)
- `clay::Float64`: clay (%)
- `sand::Float64`: sand (%)
- `dg::Float64`: geometric particle size distribution (µm)
- `accthick::Float64`: accumulated thickness of current and all soil layers above it (m)
- `depth::Float64`: depth of layer from soil surface (m)
- `sat::Float64`: saturation (m3/m3)
- `fc::Float64`: field capacity (m3/m3)
- `pwp::Float64`: permanent wilting point (m3/m3)
- `psd::Float64`: pore size distribution (slope of the log. suction-water curve)
- `porosity::Float64`: soil porosity (fraction)
- `airentry::Float64`: air entry suction (bubbling presssure) (kPa)
- `critical::Float64`: critical VWC, below which plant water stress occurs (m3/m3)
- `ksat::Float64`: saturated hydraulic conductivity (m/day)
- `vwc::Float64`: volumetric soil water content VWC (m3/m3)
- `wc::Float64`: soil water content (mm)
- `k::Float64`: unsaturated hydraulic conductivity (m/day)
- `matric::Float64`: matric head (m)
- `gravity::Float64`: gravity head (m)
- `tothead::Float64`: total (matric+gravity) head (m)
- `t::Float64`: plant water uptake (transpiration) (m/day)
- `e::Float64`: soil evaporation (m/day)
- `influx::Float64`: water entry into layer (m/day)
- `outflux::Float64`: water exit out of layer (m/day)
- `netflux::Float64`: difference between influx and outflux (m/day)
"""
@with_kw mutable struct SoilLayer
    thick::Float64 = 0.0
    clay::Float64 = 0.0
    sand::Float64 = 0.0

    dg::Float64 = 0.0
    accthick::Float64 = 0.0
    depth::Float64 = 0.0
    sat::Float64 = 0.0
    fc::Float64 = 0.0
    pwp::Float64 = 0.0
    psd::Float64 = 0.0
    porosity::Float64 = 0.0
    airentry::Float64 = 0.0
    critical::Float64 = 0.0
    ksat::Float64 = 0.0

    vwc::Float64 = 0.0
    wc::Float64 = 0.0
    k::Float64 = 0.0
    matric::Float64 = 0.0
    gravity::Float64 = 0.0
    tothead::Float64 = 0.0

    t::Float64 = 0.0
    e::Float64 = 0.0
    influx::Float64 = 0.0
    outflux::Float64 = 0.0
    netflux::Float64 = 0.0
end


"""
    RootsLayer

Soil root zone physical parameters.

# Fields
- `depth`: depth of soil layer from surface (m)
- `clay`: amount of clay (%)
- `sand`: amount of sand (%)
- `vwc::Float64`: volumetric soil water content VWC (m3/m3)
- `wc::Float64`: soil water content (mm)
- `sat::Float64`: saturation (m3/m3)
- `fc::Float64`: field capacity (m3/m3)
- `pwp::Float64`: permanent wilting point (m3/m3)
- `critical::Float64`: critical VWC, below which plant water stress occurs (m3/m3)
- `last_layer:Int`: the last layer number (index) that holds the roots
"""
@with_kw mutable struct RootsLayer
    depth::Float64 = 0.0
    clay::Float64 = 0.0
    sand::Float64 = 0.0
    vwc::Float64 = 0.0
    wc::Float64 = 0.0
    sat::Float64 = 0.0
    fc::Float64 = 0.0
    pwp::Float64 = 0.0
    critical::Float64 = 0.0
    last_layer::Int = 0
end


"""
    OP

Daily model output results.

# Fields
- `yap::Int`: years after planting (years)
- `doy::Int`: day of year
- `solardec::Float64`: solar declination (radians)
- `sunrise::Float64`: time of sunrise (hours)
- `sunset::Float64`: time of sunset (hours)
- `etrad::Float64`: daily ET solar radiation (MJ/m2)
- `totrad::Float64`: daily total solar radiation (MJ/m2)
- `drrad::Float64`: daily direct solar radiation (MJ/m2)
- `dfrad::Float64`: daily diffuse solar radiation (MJ/m2)
- `tmin::Float64`: daily min. air temperature (°C)
- `tmax::Float64`: daily max. air temperature (°C)
- `wind::Float64`: daily mean wind speed (m/s)
- `rain::Float64`: daily rain (mm)

- `treeage::Int`: tree age (days)
- `plantdens::Int`: planting density (palms/ha)
- `parts::Parts`: tree part properties
- `treehgt::Float64`: tree height (m)
- `trunkhgt::Float64`: trunk height (m)
- `assim4maint::Float64`: assimilates for maintenance (kg CH2O/palm/day)
- `assim4growth::Float64`: assimilates for growth (kg CH2O/palm/day)
- `assim4gen::Float64`: assimilates for generative growth (kg CH2O/palm/day)
- `vdmreq::Float64`: vegetative DM requirement (kg DM/palm/day)
- `vdmwgt::Float64`: vegetative DM weight (kg DM/palm)
- `tdmwgt::Float64`: total DM weight (kg DM/palm)
- `flowgt::Float64`: male and female flowers DM weight (kg DM/palm)
- `frdwgt::Float64`: fronds (pinnae + rachis) weight (kg DM/palm)
- `vdmgro::Float64`: vegetative growth rate (kg DM/palm/day)
- `tdmgro::Float64`: total growth rate (kg DM/palm/day)
- `yield::Float64`: FFB yield (kg DM/palm/day)
- `sla::Float64`: specific leaf area (m2 leaf/kg DM leaf)
- `lai::Float64`: leaf area index (m2 leaf/m2 ground)
- `laimax::Float64`: max. leaf area index (m2 leaf/m2 ground)
- `vdmmax::Float64`: max. annual vegetative DM (kg/palm/year)

- `layers::Vector{SoilLayer}`: soil layers
- `roots::RootsLayer`: root zone soil layer
- `stress_t::Float64`: reduction in potential plant transpiration (0-1)
- `stress_e::Float64`: reduction in potential soil evaporation (0-1)
- `netrain::Float64`: net rainfall (mm)
- `aet_soil::Float64`: actual soil evaporation (mm/day)
- `aet_crop::Float64`: actual plant transpiration (mm/day)
- `dp::Float64`: deep peroclation (mm/day)

- `refhgt::Float64`: reference height (m)
- `d::Float64`: zero plane displacement (m)
- `z0::Float64`: crop roughness length (m)
- `kwind::Float64`: wind extinction coefficient (unitless)
- `keddy::Float64`: eddy diffusivity extinction coefficient (unitless)
- `leaflength::Float64`: pinnae length (m)
- `leafwidth::Float64`: pinnae mean width (m)
- `efflai::Float64`: effective leaf area index (m2 leaf/m2 ground)

- `assimilates::Float64`: daily assimilates from photosynthesis (kg CH2O/palm/day)
- `et::HeatFluxes`: daily latent heat fluxes ((mm water/day)
- `h::HeatFluxes`: daily sensible heat fluxes (MJ/m2/day)
"""
@with_kw mutable struct OP
    yap::Int = 0
    date::Date = Date(1970, 1, 1)
    doy::Int = 0
    solardec::Float64 = 0.0
    sunrise::Float64 = 0.0
    sunset::Float64 = 0.0
    etrad::Float64 = 0.0
    totrad::Float64 = 0.0
    drrad::Float64 = 0.0
    dfrad::Float64 = 0.0
    tmin::Float64 = 0.0
    tmax::Float64 = 0.0
    wind::Float64 = 0.0
    rain::Float64 = 0.0

    treeage::Int = 0
    plantdens::Int = 0
    parts::Parts = Parts()
    treehgt::Float64 = 0.0
    trunkhgt::Float64 = 0.0
    assim4maint::Float64 = 0.0
    assim4growth::Float64 = 0.0
    assim4gen::Float64 = 0.0
    vdmreq::Float64 = 0.0
    tdmwgt::Float64 = 0.0
    vdmwgt::Float64 = 0.0
    flowgt::Float64 = 0.0
    frdwgt::Float64 = 0.0
    vdmgro::Float64 = 0.0
    tdmgro::Float64 = 0.0
    yield::Float64 = 0.0
    sla::Float64 = 0.0
    lai::Float64 = 0.0
    laimax::Float64 = 0.0
    vdmmax::Float64 = 0.0

    layers::Vector{SoilLayer} = []
    roots::RootsLayer = RootsLayer()
    stress_t::Float64 = 0.0
    stress_e::Float64 = 0.0
    netrain::Float64 = 0.0
    aet_soil::Float64 = 0.0
    aet_crop::Float64 = 0.0
    dp::Float64 = 0.0

    refhgt::Float64 = 0.0
    d::Float64 = 0.0
    z0::Float64 = 0.0
    kwind::Float64 = 0.0
    keddy::Float64 = 0.0
    leaflength::Float64 = 0.0
    leafwidth::Float64 = 0.0
    efflai::Float64 = 0.0

    assimilates::Float64 = 0.0
    et::HeatFluxes = HeatFluxes()
    h::HeatFluxes = HeatFluxes()
end


"""
    Output

Holder for daily model output results.

# Fields
- `op::OP`: daily model output results (can be deep copied)
- `bc::BoxCar`: flower and bunch development cycle (should not be deep copied)

# Notes
- This composite holds the parts of the model output results that can be deep copied
  and parts that should not be deep copied. Deep copying the `BoxCar` object, for
  instance, is expensive and will slow down execution; thus, this object is placed
  separately from the other model output results.
"""
@with_kw struct Output
    op::OP = OP()           # can be deep copied
    bc::BoxCar = BoxCar()   # should not be deep copied (expensive)
end

include("mo.clone.jl")
include("mo.tabulate.jl")
