"""
   parscatter()

PAR scattering coefficient (unitless).

# Returns
- `Float64`: PAR scattering coefficient (unitless)
"""
parscatter() = 0.8  # PAR scattering coefficient (unitless)


"""
   parabsorb()

PAR absorption coefficient (unitless).

# Returns
- `Float64`: PAR scattering coefficient (unitless)
"""
parabsorb() = 0.84  # PAR absorption coefficient (unitless)


"""
   parsoil()

PAR reflection off soil surface (unitless).

# Returns
- `Float64`: PAR reflection off soil surface (unitless)
"""
parsoil() = 0.15  # PAR reflection off soil surface (unitless)


"""
    quantum_yield()

Quantum efficiency/yield (µmol CO2/µmol photons).

# Returns
- `Float64`: quantum efficiency/yield (µmol CO2/µmol photons)
"""
quantum_yield() = 0.051  # quantum efficiency/yield (µmol CO2/µmol photons)


"""
    o2_ambient()

Ambient O2 level (µmol O2/mol air).

# Returns
- `Float64`: ambient O2 level (µmol O2/mol air)
"""
o2_ambient() = 210000.0   # ambient O2 level (µmol O2/mol air)


"""
    gap_fraction(lai)

Canopy gap fraction, viewed from zenith (0 = no gap/openings, 1 = full opening).

# Arguments
- `lai`: leaf area index (m2 leaf/m2 ground)

# Returns
- `Float64`: canopy gap fraction (0-1)
"""
gap_fraction(lai) = min(1, max(0, 1 - lai / 3))


"""
    vcmax25(treeage)

Maximum Rubisco capacity at 25 deg. C (µmol CO2/m2 leaf/s).

# Arguments
- `treeage`: tree age (days)

# Returns
- `Float64`: maximum Rubisco capacity at 25 deg. C (µmol CO2/m2 leaf/s)
"""
vcmax25(treeage) = exp(4.7971 - 0.9095 * exp(-treeage / 365))


"""
    canopy_extinction(solarinc, lai)

Canopy extinction coefficients (unitless) for direct and diffuse solar irradiance.

# Arguments
- `solarinc`: solar inclination (radians)
- `lai`: leaf area index (m2 leaf/m2 ground)

# Returns
- `Tuple{Float64, Float64}`:
        canopy extinction coefficients (unitless) for direct and diffuse solar irradiance
"""
function canopy_extinction(solarinc, lai)
    kdr = min(10, 0.5 / cos(solarinc))  # cap k for hours before sunrise and after sunset
    kdf = exp(0.25 - 0.52 * sqrt(lai))
    kdr, kdf
end


"""
    canopy_clump(kdr, gap, solarinc, lai)

Canopy clump (cluster) factor (0-1) (unitless).

# Arguments
- `kdr`: canopy extinction coefficient (unitless) for direct solar irradiance
- `gap`: canopy gap fraction (0-1)
- `solarinc`: solar inclination (radians)
- `lai`: leaf area index (m2 leaf/m2 ground)

# Returns
- `Float64`: canopy clump (cluster) factor (0-1) (unitless)
"""
function canopy_clump(kdr, gap, solarinc, lai)
    kdrlai = kdr * lai
    w0 = -log(gap + (1 - gap) * exp(-kdrlai / (1 - gap))) / kdrlai
    w0 + 6.6557 * (1 - w0) * exp(-exp(-solarinc + 2.2103))
end


"""
    reflection_coef(kdr, kdf, clump, lai)

Reflection coefficients (unitless) for direct and diffuse PAR.

# Arguments
- `kdr`: canopy extinction coefficient (unitless) for direct solar irradiance
- `kdf`: canopy extinction coefficient (unitless) for diffuse solar irradiance
- `clump`: canopy clump (cluster) factor (0-1) (unitless)
- `lai`: leaf area index (m2 leaf/m2 ground)

# Returns
- `Tuple{Float64, Float64}`: reflection coefficients (unitless) for direct and diffuse PAR
"""
function reflection_coef(kdr, kdf, clump, lai)
    a = sqrt(parscatter()) * lai
    ps = parsoil()
    pdr = max(0.04, ps * exp(-2 * kdr * clump * a))
    pdf = max(0.04, ps * exp(-2 * kdf * a))
    pdr, pdf
end


"""
    lai_components(kdr, clump, lai)

Distinguish total LAI (leaf area index) into sunlit and shaded LAI (m2 leaf/m2 ground).

# Arguments
- `kdr`: canopy extinction coefficient (unitless) for direct solar irradiance
- `clump`: canopy clump (cluster) factor (0-1) (unitless)
- `lai`: leaf area index (m2 leaf/m2 ground)

# Returns
- `Tuple{Float64, Float64}`:
        sunlit and shaded LAI (m2 leaf/m2 ground)
"""
function lai_components(kdr, clump, lai)
    a = kdr * clump
    lsl = (1 - exp(-a * lai)) / a
    lsh = lai - lsl
    lsl, lsh
end


"""
    par_components(drrad, dfrad, kdr, kdf, clump, lai)

PAR components within canopies (µmol photons/m2 leaf/s).

# Arguments
- `drrad`: direct solar irradiance outside canopies (W/m2)
- `dfrad`: diffuse solar irradiance outside canopies (W/m2)
- `kdr`: canopy extinction coefficient (unitless) for direct solar irradiance
- `kdf`: canopy extinction coefficient (unitless) for diffuse solar irradiance
- `clump`: canopy clump (cluster) factor (0-1) (unitless)
- `lai`: leaf area index (m2 leaf/m2 ground)

# Returns
- `PAR`: PAR components within canopies (µmol photons/m2 leaf/s)
"""
function par_components(drrad, dfrad, kdr, kdf, clump, lai)
    # outside PAR (µmol photons/m2 ground/s):
    # 50% of solar radiation = PAR, and 1 W/m2 = 4.55 µmol photons/m2/s
    outdr = drrad * 0.5 * 4.55
    outdf = dfrad * 0.5 * 4.55
    a = kdr * clump * lai
    b = sqrt(parscatter())

    # within canopies PAR (mol photons/m2 (ground or leaf)/s):
    # µmol photons/m2 ground/s
    pdr, pdf = reflection_coef(kdr, kdf, clump, lai)
    indrscatter = (1 - pdr) * outdr * exp(-a * b)
    indr = (1 - pdr) * outdr * exp(-a)  # µmol photons/m2 ground/s
    inscatter = 0.5 * (indrscatter - indr)  # µmol photons/m2 leaf/s
    a = kdf * b * lai
    indf = ((1 - pdf) * outdf * (1 - exp(-a))) / a  # µmol photons/m2 leaf/s

    # absorbed PAR by sunlit and shaded leaves (µmol photons/m2 leaf/s):
    pa = parabsorb()
    abssunlit = pa * (kdr * clump * outdr + indf + inscatter)
    absshaded = pa * (indf + inscatter)

    PAR(
        outdr = outdr,
        outdf = outdf,
        indrscatter = indrscatter,
        indr = indr,
        inscatter = inscatter,
        indf = indf,
        abssunlit = abssunlit,
        absshaded = absshaded
    )
end


"""
   assim_coefs!(ophr::OPHour, treeage)

Set the temperature-dependent CO2 assimilation parameters/coefficients.

# Arguments
- `ophr::OPHour`: instantaneous/hourly properties
- `treeage`: tree age (days)

# Returns
- nothing
"""
function assim_coefs!(ophr::OPHour, treeage)
    # lists for Kc, Ko, specificity, then Vcmax
    # (each parameter's 25 deg C value then its Q10 value)
    tf = ophr.canopy_temp
    vcmax = vcmax25(treeage)
    params = [(270, 2.786), (165000, 1.355), (2388.6, 0.700), (vcmax, 2.5319)]
    adjp = map(p -> p[1] * p[2]^((tf - 25) / 10), params)

    # additional temperature correction just for Vcmax
    adjp[end] /= 1 + exp(0.2006 * (tf - 45.7521))
    append!(adjp, 0.5 * o2_ambient() / adjp[3])  # add CO2 compensation point

    ophr.assimcoef = AssimCoef(
        mmco2 = adjp[1],
        mmo2 = adjp[2],
        specificity = adjp[3],
        vcmax = adjp[4],
        co2_pt = adjp[5]
    )
end


"""
    co2_internal!(ophr::OPHour)

Internal (intercellular) CO2 concentration in plant (µmol CO2/mol air).

# Arguments
- `ophr::OPHour`: instantaneous/hourly properties

# Returns
- nothing
"""
function co2_internal!(ophr::OPHour)
    @unpack co2_ambient, assimcoef, vp, canopy_temp = ophr
    co2_pt = assimcoef.co2_pt

    # vapor pressure deficit in leaf: > 53 mbar, after which full stomatal closure
    vpdleaf = min(53.0, svp_fn(canopy_temp) - vp)
    a, b = 0.0615, 0.0213  # empirical coefficients for intercellular CO2 (Ci) vs VPD
    ophr.co2_internal = co2_ambient * (1 - (1 - co2_pt / co2_ambient) * (a + b * vpdleaf))
end


"""
    plantrad_regime!(ophr::OPHour, lai)

Plant radiation regime and portion of sunlit and shaded leaf area index.

# Arguments
- `ophr::OPHour`: instantaneous/hourly properties
- `lai`: leaf area index (m2 leaf/m2 ground)

# Returns
- nothing
"""
function plantrad_regime!(ophr::OPHour, lai)
    @unpack solarinc, drrad, dfrad = ophr

    gap = gap_fraction(lai)
    kdr, kdf = canopy_extinction(solarinc, lai)
    clump = canopy_clump(kdr, gap, solarinc, lai)
    lsl, lsh = lai_components(kdr, clump, lai)
    par = par_components(drrad, dfrad, kdr, kdf, clump, lai)
    @pack! ophr = kdr, kdf, clump, gap, lsl, lsh, par
end


"""
    leaf_assim!(ophr::OPHour)

CO2 assimilation per leaf area basis (µmol CO2/m2 leaf/s).

# Arguments
- `ophr::OPHour`: instantaneous/hourly properties

# Returns
- nothing
"""
function leaf_assim!(ophr::OPHour)
    @unpack par, assimcoef, co2_internal = ophr

    # no photosynthesis if internal CO2 is less than CO2 compensation point
    co2_diff = max(0.0, co2_internal - assimcoef.co2_pt)
    vc = assimcoef.vcmax * co2_diff
    n = 1 + o2_ambient() / assimcoef.mmo2
    vc /= assimcoef.mmco2 * n + co2_internal
    a = co2_diff / (co2_internal + 2 * assimcoef.co2_pt)
    a *= quantum_yield() * parabsorb()
    vqsl, vqsh = par.abssunlit * a, par.absshaded * a
    vs = assimcoef.vcmax * 0.5

    fm = findmin([vc, vqsl + vqsh, vs])     # find the most limiting factor
    if fm[2] == 2
        sunlit = vqsl
        shaded = vqsh
    else
        sunlit = shaded = fm[1]
    end

    ophr.leafassim = LeafAssim(
        vc = vc,
        vqsl = vqsl,
        vqsh = vqsh,
        vs = vs,
        sunlit = sunlit,
        shaded = shaded
    )
end


"""
    canopy_assim!(ophr::OPHour)

Scale up CO2 leaf assimilation to canopy level (µmol CO2/m2 ground/s).

# Arguments
- `ophr::OPHour`: instantaneous/hourly properties

# Returns
- nothing
"""
function canopy_assim!(ophr::OPHour)
    @unpack lsl, lsh, leafassim = ophr
    ophr.canopyassim = (leafassim.sunlit * lsl) + (leafassim.shaded * lsh)
end


"""
    soil_plant_atmos!(mo::Output, mi::Input)

Instantaneous canopy CO2 assimilation and energy balance of the
soil-plant-atmosphere system (SPAS).

# Arguments
- `mo::Output`: model output daily results
- `mi::Input`: model input

# Notes
- Example:

  spas = soil_plant_atmos!(mo, mi)  # set the daily properties
  spas(10)   # determine the SPAS for the local solar hour 10

  This method could have been written as:

  soil_plant_atmos!(10, mo, mi)

  but this would mean the daily properties would be set for
  each time this method is called for different hours on the
  same day. This is inefficient as the daily properties do
  not change within the same day.

  This method is used in integration over a certain period
  within the same day.

# Returns
- This method returns a function `solve!` which returns a
  `Tuple{OPHour, Float64, Float64, Float64, Float64, Float64, Float64, Float64}` for:
        instantaneous/hourly properties,
        canopy CO2 assimilation (umol CO2/m2 leaf/s),
        latent heat fluxes (total, crop, and soil) (W/m2), and
        sensible heat fluxes (total, crop, and soil) (W/m2)
"""
function soil_plant_atmos!(mo::Output, mi::Input)
    op = mo.op
    @unpack treehgt, treeage, plantdens, lai, date = op

    co2_a = mi.ccfn.co2(date)  # ambient CO2 level (µmol CO2/mol air)
    op.refhgt = reference_height(treehgt)
    op.efflai = effective_lai(lai)
    op.d, op.z0, op.kwind, op.keddy = wind_profile(op.efflai, treehgt)

    hour_wthr = est_hour_wthr(mi.lat, mo)   # hourly weather function object

    function solve!(th)   # th: local solar hour (hours)
        # 1. calculate the meteorological properties for local solar hour, th:
        ophr = hour_wthr(th)

        # 2. calculate the plant radiation regime:
        plantrad_regime!(ophr, lai)

        # 3. calculate the heat fluxes in the energy balance:
        update_energybal!(ophr, mo, mi)

        # 4. canopy temperature is now known, so finally calculate the photosynthesis:
        @unpack canopy_temp, vp = ophr
        assim_coefs!(ophr, treeage)
        ophr.co2_ambient = co2_a    # constant for whole day
        co2_internal!(ophr)
        leaf_assim!(ophr)
        canopy_assim!(ophr)

        (
            ophr = ophr,
            canopyassim = ophr.canopyassim,
            et_total = ophr.et.total,
            et_crop = ophr.et.crop,
            et_soil = ophr.et.soil,
            h_total = ophr.h.total,
            h_crop = ophr.h.crop,
            h_soil = ophr.h.soil,
        )
    end
end


"""
    update_photosyn!(mo::Output, mi::Input, date::Date)

Daily photosynthesis and heat fluxes (latent and sensible) for the current date.

# Arguments
- `mo::Output`: model output daily results
- `mi::Input`: model input
- `date::Date`: date

# Returns
- nothing
"""
function update_photosyn!(mo::Output, mi::Input)
    spas = soil_plant_atmos!(mo, mi)
    fields = [:canopyassim, :et_total, :et_crop, :et_soil, :h_total, :h_crop, :h_soil]
    itg = integrate(9, 0, 24)
    res = itg(spas; fields=fields)

    op = mo.op
    res[1] *= 1.08 / op.plantdens  # µmol CO2/m2 leaf/day to kg CH2O/palm/day
    res[2:4] .*= 3600 / 2454000    # latent fluxes in mm/day
    res[5:7] .*= 3600 / 10^6       # sensible fluxes in MJ/m2/day

    op.assimilates = res[1]
    op.et = HeatFluxes(total=res[2], crop=res[3], soil=res[4])
    op.h = HeatFluxes(total=res[5], crop=res[6], soil=res[7])
end
