"""
   PAR

Components of photosynthetically active radiation (PAR) within and outside the canopies.

# Notes
- Units are µmol photons on a per ground or leaf area basis, depending on parameter.

# Fields
- `outdr::Float64`: direct PAR outside the canopies (µmol photons/m2 ground/s)
- `outdf::Float64`: diffuse PAR outside the canopies (µmol photons/m2 ground/s)
- `indrscatter::Float64`: direct PAR and scatter within the canopies
                          (µmol photons/m2 ground/s)
- `indr::Float64`: direct PAR within the canopies (µmol photons/m2 ground/s)
- `inscatter::Float64`: scatter within the canopies (µmol photons/m2 leaf/s)
- `indf::Float64`: diffuse PAR within the canopies (µmol photons/m2 leaf/s)
- `abssunlit::Float64`: PAR absorbed by sunlit leaves (µmol photons/m2 leaf/s)
- `absshaded::Float64`: PAR absorbed by shaded leaves (µmol photons/m2 leaf/s)
"""
@with_kw struct PAR
    outdr::Float64 = 0.0
    outdf::Float64 = 0.0
    indrscatter::Float64 = 0.0
    indr::Float64 = 0.0
    inscatter::Float64 = 0.0
    indf::Float64 = 0.0
    abssunlit::Float64 = 0.0
    absshaded::Float64 = 0.0
end


"""
   AssimCoef

Photosynthesis coefficients.

# Fields
- `mmco2::Float64`: Michaelis-Menten constant for CO2 (µmol/mol)
- `mmo2::Float64`: Michaelis-Menten constant for O2 (µmol/mol)
- `specificity::Float64`: CO2 / O2 specificity factor (unitless)
- `vcmax::Float64`: Rubisco maximum capacity rate (µmol CO2/m2 leaf/s)
- `co2_pt::Float64`: CO2 compensation point (µmol CO2/mol)
"""
@with_kw struct AssimCoef
    mmco2::Float64 = 0.0
    mmo2::Float64 = 0.0
    specificity::Float64 = 0.0
    vcmax::Float64 = 0.0
    co2_pt::Float64 = 0.0
end


"""
   LeafAssim

Leaf CO2 assimilation components.

# Notes
- All units in µmol CO2/m2 leaf/s

# Fields
- `vc::Float64`: Rubisco-limited assimilation
- `vqsl::Float64`: light-limited assimilation by sunlit leaves
- `vqsh::Float64`: light-limited assimilation by shaded leaves
- `vs::Float64`: sink-limited assimilation
- `sunlit::Float64`: total assimilation by sunlit leaves
- `shaded::Float64`: total assimilation by shaded leaves
"""
@with_kw struct LeafAssim
    vc::Float64 = 0.0
    vqsl::Float64 = 0.0
    vqsh::Float64 = 0.0
    vs::Float64 = 0.0
    sunlit::Float64 = 0.0
    shaded::Float64 = 0.0
end


"""
    StomatalStresses

Stresses that would reduce stomatal conductivity (0-1).

# Notes
- All parameters are unitless, between 0-1.

# Fields
- `water::Float64`: reduction due to water stress
- `vpd::Float64`: reduction due to vapor pressure deficit
- `par::Float64`: reduction due to low PAR solar irradiance
"""
@with_kw mutable struct StomatalStresses
    water::Float64 = 0.0
    vpd::Float64 = 0.0
    par::Float64 = 0.0
end


"""
    AvailEnergy

Available radiation to the soil-plant-atmosphere system (W/m2).

# Notes
- All parameters are in W/m2.

# Fields
- `total::Float64`: radiation available to both crop and soil (total)
- `crop::Float64`: radiation available to the crop
- `soil::Float64`: radiation available to the soil
- `net::Float64`: net radiation
- `g::Float64`: soil heat flux
"""
@with_kw mutable struct AvailEnergy
    total::Float64 = 0.0
    crop::Float64 = 0.0
    soil::Float64 = 0.0
    net::Float64 = 0.0
    g::Float64 = 0.0
end


"""
     Resistances

 Resistances to heat fluxes (s/m).

# Notes
- All parameters are in s/m.

# Fields
- `rsa::Float64`: aerodynamic resistance between soil and mcf (mean canopy flow)
- `raa::Float64`: aerodynamic resistance between mcf and reference level
- `rca::Float64`: boundary layer resistance
- `rst::Float64`: leaf stomatal resistance
- `rcs::Float64`: canopy resistance
- `rss::Float64`: soil resistance
"""
@with_kw mutable struct Resistances
    rsa::Float64 = 0.0
    raa::Float64 = 0.0
    rca::Float64 = 0.0
    rst::Float64 = 0.0
    rcs::Float64 = 0.0
    rss::Float64 = 0.0
end


"""
    HeatFluxes

Heat flux components.

# Notes
- Units vary depending on use (e.g., can be W/m2, mm/day, or MJ/m2/day).

# Fields
- `total::Float64`: total heat flux (crop + soil)
- `crop::Float64`: heat flux from the crop
- `soil::Float64`: heat flux from the soil
"""
@with_kw struct HeatFluxes
    total::Float64 = 0.0
    crop::Float64 = 0.0
    soil::Float64 = 0.0
end


"""
    OPHour

Instantaneous (hourly) model output results.

# Fields
- `date::Date`: date
- `doy::Int`: day of year
- `solarhour::Float64`: solar hour (hours)
- `air_temp::Float64`: air temperature (°C)
- `dew_temp::Float64`: dew temperature (°C)
- `svp::Float64`: saturated air vapor pressure (mbar)
- `vp::Float64`: air vap. pressure (mbar)
- `rh::Float64`: relative hunidity (%)
- `vpd::Float64`: vap. pressure deficit (VPD) (mbar)
- `vpd0::Float64`: VPD at mean canopy flow height (mbar)
- `svp_slope::Float64`: saturated air vapor pressure slope (mbar/°C)
- `solarinc::Float64`: solar inclination (radians)
- `solarhgt::Float64`: solar height (radians)
- `solarazi::Float64`: solar azimuth (radians)
- `solarcon::Float64`: solar constanrt (W/m2)
- `etrad::Float64`: ET solar radiation (W/m2)
- `totrad::Float64`: total solar radiation (W/m2)
- `drrad::Float64`: direct solar radiation (W/m2)
- `dfrad::Float64`: diffuse solar radiation (W/m2)
- `netrad::Float64`: net solar radiation (W/m2)
- `wind::Float64`: wind speed (m/s)

- `kdr::Float64`: canopy extinction coefficient for direct radiation (0-1)
- `kdf::Float64`: canopy extinction coefficient for diffuse radiation (0-1)
- `clump::Float64`: canopy clump factor (0-1)
- `gap::Float64`: canopy gap fraction, viewed from zenith (fraction, 0-1)
- `lsl::Float64`: sunlit leaf area index (m2 leaf / m2 area)
- `lsh::Float64`: shaded leaf area index (m2 leaf / m2 area)
- `par::PAR`: PAR components within canopies (µmol photons/m2 leaf/s)

- `co2_ambient::Float64`: ambient CO2 level (µmol CO2/mol air)
- `co2_internal::Float64`: intercellular CO2 level (µmol CO2/mol air)
- `assimcoef::AssimCoef`: photosynthesis coefficients (units depend on parameter)
- `leafassim::LeafAssim`: leaf CO2 assimilation components (µmol CO2/m2 leaf/s)
- `canopyassim::Float64`: canopy CO2 assimilation (µmol CO2/m2 ground/s)

- `A::AvailEnergy`: available energy to the system (W/m2)
- `ustar::Float64`: friction velocity (m/s)
- `utreehgt::Float64`: wind speed at tree height (m/s)
- `canopy_temp::Float64`: canopy temperature (°C)
- `stress::StomatalStresses`: reduction to stomatal conductivity due to stresses (0-1)
- `R::Resistances`: resistances to fluxes (s/m)
- `et::HeatFluxes`: latent heat flux components (W/m2)
- `h::HeatFluxes`: sensible heat flux components (W/m2)
"""
@with_kw mutable struct OPHour
    date::Date = Date(1970, 1, 1)
    doy::Int = 0
    solarhour::Float64 = 0.0
    air_temp::Float64 = 0.0
    dew_temp::Float64 = 0.0
    svp::Float64 = 0.0
    vp::Float64 = 0.0
    rh::Float64 = 0.0
    vpd::Float64 = 0.0
    vpd0::Float64 = 0.0
    svp_slope::Float64 = 0.0
    solarinc::Float64 = 0.0
    solarhgt::Float64 = 0.0
    solarazi::Float64 = 0.0
    solarcon::Float64 = 0.0
    etrad::Float64 = 0.0
    totrad::Float64 = 0.0
    drrad::Float64 = 0.0
    dfrad::Float64 = 0.0
    netrad::Float64 = 0.0
    wind::Float64 = 0.0

    kdr::Float64 = 0.0
    kdf::Float64 = 0.0
    clump::Float64 = 0.0
    gap::Float64 = 0.0
    lsl::Float64 = 0.0
    lsh::Float64 = 0.0
    par::PAR = PAR()

    co2_ambient::Float64 = 0.0
    co2_internal::Float64 = 0.0
    assimcoef::AssimCoef = AssimCoef()
    leafassim::LeafAssim = LeafAssim()
    canopyassim::Float64 = 0.0

    A::AvailEnergy = AvailEnergy()
    ustar::Float64 = 0.0
    utreehgt::Float64 = 0.0
    canopy_temp::Float64 = 0.0
    stresses::StomatalStresses = StomatalStresses()
    R::Resistances = Resistances()
    et::HeatFluxes = HeatFluxes()
    h::HeatFluxes = HeatFluxes()
end
