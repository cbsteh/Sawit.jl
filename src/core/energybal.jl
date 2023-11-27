"""
    vonK()

von Karman constant

# Returns
- `Float64`: von Karman constant (unitless)
"""
vonK() = 0.41


"""
    vol_heat_capacity()

Volumentric heat capacity (J/m3/K)

# Returns
- `Float64`: volumentric heat capacity (J/m3/K)
"""
vol_heat_capacity() = 1221.09


"""
    soil_roughlen()

Soil roughness length (m).

# Returns
- `Float64`: soil roughness length (m)
"""
soil_roughlen() = 0.0112


"""
    reference_height(treehgt)

Reference height (m)

# Arguments:
- `treehgt`: total tree height (m)

# Returns
- `Float64`: reference height (m)
"""
reference_height(treehgt) = treehgt + 0.5


"""
    effective_lai(lai)

Effective leaf area index (m2 leaf/m2 ground)

# Arguments:
- `lai`: leaf area index (m2 leaf/m2 ground)

# Returns
- `Float64`: effective leaf area index (m2 leaf/m2 ground)
"""
function effective_lai(lai)
    (lai < 2) ? lai :
    (lai > 4) ? 0.5 * lai :
    2.0
end


"""
    wind_profile(efflai, treehgt)

Vertical wind speed profile paramaters.

# Arguments
- `efflai`: effective leaf area index (m2 leaf/m2 ground)
- `treehgt`: total tree height (m)

# Returns
- `Tuple{Float64, Float64, Float64, Float64}`:
        zero plane displacement (m), crop roughness length (m), and
        wind and eddy extinction coefficients (unitless)
"""
function wind_profile(efflai, treehgt)
    # wind speed extinction coefficient (Massman, 1987):
    cd = 0.2    # mean drag coefficient
    cdlai = cd * efflai  # effective LAI used
    y = max(0.1, -0.146 + 0.319 / sqrt(cdlai))   # fitted curve for alpha=2
    kwind = min(5, 1 / y)   # no more than 5 (usual values: 2 to 5)

    n = 2 * kwind
    ustar2u_h = 0.32 - (0.264 * exp(-15.1 * cdlai))
    dhd = min(0.95, max(0.3, 1 - (exp(-n) / n * (exp(n) - 1))))
    d = dhd * treehgt     # zero plane displacement height, d (m)
    dhz0 = (1 - dhd) * exp(-vonK() / ustar2u_h)
    z0 = dhz0 * treehgt   # crop roughness length, z0 (m)

    d, z0, kwind, kwind      # assume eddy = wind extinction coefficient
end


"""
    stomatal_stresses(vpd, totrad, stress_t)

Reduction in stomatal conductance due to various stresses (unitless).

# Arguments:
- `vpd`: vapor pressure deficit (mbar)
- `totrad`: total solar irradiance (W/m2)
- `stress_t`: crop water stress (fraction, 0-1)

# Returns
- `StomatalStresses`: stresses that will reduce stomatal conductance (unitless).
"""
function stomatal_stresses(vpd, totrad, stress_t)
    gst_vpd(vp) = -0.01111 * log(vp) + 0.044172
    gst_par(pr) = 0.014614 * (1 - exp(-0.008740 * pr))

    gstmin = gst_vpd(53)     # VPD for min. conductance
    gstmax = gst_vpd(10)     # VPD for max. conductance
    gst = min(gstmax, max(gstmin, gst_vpd(max(10, vpd))))
    svpd = gst / gstmax

    gstmin = gst_par(0.1)  # PAR for min. conductance (not 0, avoid divide-by-zero later)
    gstmax = gst_par(330)  # PAR for max. conductance
    totpar = totrad * 0.5  # total PAR in W/m2
    gst = min(gstmax, max(gstmin, gst_par(totpar)))
    spar = gst / gstmax

    StomatalStresses(
        water = stress_t,
        vpd = svpd,
        par = spar
    )
end


"""
    available_energy(dg, netrad, kdr, clump, lai)

Available energy (W/m2) to crop and soil (net radiation partitioning).

# Arguments:
- `dg::Float64`: first (top) soil layer geometric particle size distribution (µm)
- `netrad`: net irradiance (W/m2)
- `kdr`: canopy extinction coefficient for direct radiation (unitless)
- `clump`: canopy clump (cluster) factor (0-1) (unitless)
- `lai`: leaf area index (m2 leaf/m2 ground)

# Returns
- `AvailEnergy`: available energy to crop and soil (W/m2)
"""
function available_energy(dg, netrad, kdr, clump, lai)
    # Rn:G clayey soils (small dg) > Rn:G sandy soils (large dg)
    ts = 0.3 + 0.274 / sqrt(dg)   # Rn:G for bare soils (0.32-0.35 in most cases)
    tc = (ts - 0.3) + 0.05        # Rn:G under full canopy (0.07-0.1 in most cases)
    alphaL = sqrt(0.5) * lai
    gap = exp(-kdr * clump * alphaL) * exp(-0.5 * alphaL)

    crop = (1 - gap) * (1 - tc) * netrad  # Rn portion to crop
    soil = gap * (1 - ts) * netrad  # Rn portion to soil
    total = crop + soil  # total Rn to both crop and soil
    g = (tc + gap * (ts - tc)) * netrad  # Rn portion as soil heat flux

    AvailEnergy(
        total = total,
        crop = crop,
        soil = soil,
        net = netrad,
        g = g
    )
end


"""
    windspd_at_refhgt(refhgt, wind, methgt)

Wind speed at reference height (m/s).

# Arguments:
- `refhgt`: reference height (m)
- `wind speed`: instantaneous wind speed (m/s)
- `methgt`: height of weather station, where wind speed measurement (m)

# Returns
- `Float64`: wind speed at reference height (m/s)
"""
function windspd_at_refhgt(refhgt, wind, methgt)
    zs0 = soil_roughlen()    # soil roughness length (m)
    # extrapolation by log law to determine wind speed at the reference height
    return wind * log(refhgt / zs0) / log(methgt / zs0)
end


"""
    friction_velocity(d, z0, refhgt, wind_refhgt)

Friction velocity (m/s).

# Arguments:
- `d`: zero plane displacement (m)
- `z0`: crop roughness length (m)
- `refhgt`: reference height (m)
- `wind_refhgt`: wind speed at reference height (m)

# Returns
- `Float64`: friction velocity (m/s)
"""
function friction_velocity(d, z0, refhgt, wind_refhgt)
    vonK() * wind_refhgt / log((refhgt - d) / z0)
end


"""
    windspd_at_treehgt(d, z0, ustar, treehgt)

Wind speed at tree height (m/s).

# Arguments:
- `d`: zero plane displacement (m)
- `z0`: crop roughness length (m)
- `ustar`: friction velocity (m/s)
- `treehgt`: total tree height (m)

# Notes
-  Wind speed measured at weather station height may not be the same as the
   reference height (thus, wind speed will have to be extrapolated to the
   reference height).

# Returns
- `Float64`: wind speed at tree height (m/s)
"""
function windspd_at_treehgt(d, z0, ustar, treehgt)
    (ustar / vonK()) * log((treehgt - d) / z0)   # extrapolation by log law
end


"""
    res_rss(firstlayer::SoilLayer)

Soil surface resistance (s/m).

# Arguments:
- `firstlayer::SoilLayer`: first (top) soil layer

# Returns
- `Float64`: soil surface resistance (s/m)
"""
function res_rss(firstlayer::SoilLayer)
    a = firstlayer.porosity
    tau = sqrt(a + 3.79 * (1 - a))  # tortuosity
    dmv = 2.47e-5     # vapor diffusion coefficient, m2/s
    rssmax = tau * firstlayer.thick / (a * dmv)
    rssmax * exp(-firstlayer.vwc / (firstlayer.psd * firstlayer.sat))
end


"""
    res_rsa(d, z0, ustar, keddy, treehgt)

Aerodynamic resistance between soil and mcf (mean canopy flow) (s/m).

# Arguments:
- `d`: zero plane displacement (m)
- `z0`: crop roughness length (m)
- `ustar`: friction velocity (m/s)
- `keddy`: eddy extinction coefficient (unitless)
- `treehgt`: total tree height (m)

# Returns
- `Float64`: aerodynamic resistance between soil and mean canopy flow (s/m)
"""
function res_rsa(d, z0, ustar, keddy, treehgt)
    a = exp(keddy) / (keddy * vonK() * ustar)
    b = exp(-keddy * soil_roughlen() / treehgt)
    c = exp(-keddy * (z0 + d) / treehgt)
    a * (b - c)
end


"""
    res_raa(d, z0, ustar, keddy, refhgt, treehgt)

Aerodynamic resistance between mcf (mean canopy flow) and reference height (s/m).

# Arguments:
- `d`: zero plane displacement (m)
- `z0`: crop roughness length (m)
- `ustar`: friction velocity (m/s)
- `keddy`: eddy extinction coefficient (unitless)
- `refhgt`: reference height (m)
- `treehgt`: total tree height (m)

# Returns
- `Float64`: aerodynamic resistance between mean canopy flow and reference height (s/m)
"""
function res_raa(d, z0, ustar, keddy, refhgt, treehgt)
    a = vonK() * ustar
    b = log((refhgt - d) / (treehgt - d)) / a
    c = 1 - (z0 + d) / treehgt
    d = (exp(keddy * c) - 1) / (keddy * a)
    b + d
end


"""
    res_rca(leafwidth, efflai, kwind, wind_treehgt)

Boundary layer resistance (s/m).

# Arguments:
- `leafwidth`: mean length of leaflets (pinnae) (m)
- `efflai`: effective leaf area index (m2 leaf/m2 ground)
- `kwind`: wind speed extinction coefficient (unitless)
- `wind_treehgt`: wind speed at tree height (m/s)

# Returns
- `Float64`: boundary layer resistance (s/m)
"""
function res_rca(leafwidth, efflai, kwind, wind_treehgt)
    a = (1 - exp(-0.5 * kwind)) * sqrt(wind_treehgt / leafwidth)
    kwind / (0.01 * efflai * a)
end


"""
    res_rcs_st(efflai, stress::StomatalStresses)

Stomatal and canopy resistances (s/m).

# Arguments:
- `efflai`: effective leaf area index (m2 leaf/m2 ground)
- `stress::StomatalStresses`: reduction to stomatal conductance due to stresses (0-1)

# Returns
- `Tuple{Float64, Float64}`: stomatal and canopy resistances (s/m)
"""
function res_rcs_st(efflai, stress::StomatalStresses)
    gstmax = 0.0125  # max. stomatal conductance (m/s) (equivalent to 500 mmol/m2/s)
    gstmin = 0.0002  # min. stomatal conductance (m/s) for C3 plants (10 mmol/m2/s)
    gst = max(gstmin, gstmax * stress.water * stress.vpd * stress.par)
    gcs = gst * efflai
    return 1 / gst, 1 / gcs
end


"""
    resistances(ophr::OPHour, mo::Output)

Flux resistances (s/m).

# Arguments
- `ophr::OPHour`: model output hourly (instantaneous) results
- `mo::Output`: model output daily results

# Returns
- `Resistances`: flux resistancess (s/m)
"""
function resistances(ophr::OPHour, mo::Output)
    @unpack d, z0, keddy, kwind, refhgt, leafwidth, efflai, treehgt, layers = mo.op
    @unpack ustar, utreehgt, stresses = ophr

    rst, rcs = res_rcs_st(efflai, stresses)
    Resistances(
        rsa = res_rsa(d, z0, ustar, keddy, treehgt),
        raa = res_raa(d, z0, ustar, keddy, refhgt, treehgt),
        rca = res_rca(leafwidth, efflai, kwind, utreehgt),
        rst = rst,
        rcs = rcs,
        rss = res_rss(layers[1])
    )
end


"""
    solve_fluxes(R::Resistances, A::AvailEnergy, Δ, vpd)

Latent and sensible heat fluxes (W/m2), based on the electrical network analogy.

# Arguments
- `R::Resistances`: flux resistances (s/m)
- `A::AvailEnergy`: available energy (W/m2) to crop and soil
- `Δ`: saturated air vapor pressure slope (mbar/°C)
- `vpd`: air vapor pressure deficit at weather station height (mbar)

# Notes
- Vapor pressure deficit (VPD) at the mean canopy flow is also returned.

# Returns
- `Tuple{HeatFluxes, HeatFluxes, Float64}`:
     latent and sensible heat fluxes (W/m2) and VPD at mean canopy flow (mbar)
"""
function solve_fluxes(R::Resistances, A::AvailEnergy, Δ, vpd)
    @unpack raa, rca, rsa, rcs, rss = R
    @unpack total, crop, soil = A

    psycho = 0.658  # psychometric constant mbar/K
    pcp = vol_heat_capacity()

    ra = (Δ + psycho) * raa
    rc = (Δ + psycho) * rca + psycho * rcs
    rs = (Δ + psycho) * rsa + psycho * rss
    cc = 1 / (1 + rc * ra / (rs * (rc + ra)))
    cs = 1 / (1 + rs * ra / (rc * (rs + ra)))
    pmc = Δ * total + (pcp * vpd - Δ * rca * soil) / (raa + rca)
    pmc /= Δ + psycho * (1 + rcs / (raa + rca))
    pms = Δ * total + (pcp * vpd - Δ * rsa * crop) / (raa + rsa)
    pms /= Δ + psycho * (1 + rss / (raa + rsa))
    et = cc * pmc + cs * pms
    vpd0 = vpd + (raa / pcp) * (Δ * total - (Δ + psycho) * et)
    etc = Δ * crop + pcp * vpd0 / rca
    etc /= Δ + psycho * (rcs + rca) / rca
    ets = Δ * soil + pcp * vpd0 / rsa
    ets /= Δ + psycho * (rss + rsa) / rsa
    hc = psycho * crop * (rcs + rca) - pcp * vpd0
    hc /= Δ * rca + psycho * (rcs + rca)
    hs = psycho * soil * (rss + rsa) - pcp * vpd0
    hs /= Δ * rsa + psycho * (rss + rsa)
    h = hc + hs

    latent = HeatFluxes(total=et, crop=etc, soil=ets)
    sensible = HeatFluxes(total=h, crop=hc, soil=hs)
    latent, sensible, vpd0
end


"""
   canopy_temp(h::HeatFluxes, rca, raa, ta)

Canopy (foliage) temperature (°C).

# Arguments
- `h::HeatFluxes`: crop and soil sensible heat fluxes
- `rca`: boundary layer resistance (s/m)
- `raa`: aerodynamic resistance between mean canopy flow and reference height (s/m)
- `ta`: air temperature (°C)

# Returns
- `Float64`: canopy temperature (°C)
"""
function canopy_temp(h::HeatFluxes, rca, raa, ta)
    delta = h.crop * rca + (h.soil + h.crop) * raa
    delta / vol_heat_capacity() + ta
end


"""
    update_energybal!(ophr::OPHour, mo::Output, mi::Input)

Instantaneous energy balance (W/m2).

# Arguments
- `ophr::OPHour`: model output hourly (instantaneous) results
- `mo::Output`: model output daily results
- `mi::Input`: model input

# Returns
- nothing
"""
function update_energybal!(ophr::OPHour, mo::Output, mi::Input)
    @unpack svp_slope, vpd, air_temp, totrad, netrad, kdr, clump, wind = ophr
    @unpack d, z0, refhgt, treehgt, lai, stress_t, layers = mo.op

    methgt = mi.methgt
    dg = layers[1].dg

    # solve the energy balance:
    ophr.stresses = stomatal_stresses(vpd, totrad, stress_t)
    ophr.A = A = available_energy(dg, netrad, kdr, clump, lai)
    wind_refhgt = windspd_at_refhgt(refhgt, wind, methgt)
    ophr.ustar = friction_velocity(d, z0, refhgt, wind_refhgt)
    ophr.utreehgt = windspd_at_treehgt(d, z0, ophr.ustar, treehgt)
    ophr.R = R = resistances(ophr, mo)
    ophr.et, ophr.h, ophr.vpd0 = solve_fluxes(R, A, svp_slope, vpd)
    ophr.canopy_temp = canopy_temp(ophr.h, R.rca, R.raa, air_temp)
end
