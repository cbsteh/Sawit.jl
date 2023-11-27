"""
   to_dg(c, s)

Convert soil texture into its geometric particle size distribution (µm).

# Arguments
- `c`: clay (%)
- `s`: sand (%)

# Returns
- `Float64`: geometric particle size distribution (µm)
"""
to_dg(c, s) = exp(-0.0196 * c + 0.023 * (100 - s - c) + 0.0576 * s)


"""
   est_airentry(dg)

Air entry suction (bubbling presssure) (kPa).

# Arguments
- `dg`:  geometric particle size distribution (µm)

# Returns
- `Float64`: air entry suction (bubbling presssure) (kPa)
"""
est_airentry(dg) = 3.9 - 0.61 * log(dg)


"""
   est_psd(dg)

 Pore size distribution (slope of the log. suction-water curve).

# Arguments
- `dg`:  geometric particle size distribution (µm)

# Returns
- `Float64`: pore size distribution (slope of the log. suction-water curve)
"""
est_psd(dg) = 1 / (8.25 - 1.26 * log(dg))


"""
   est_critical(sat, pwp)

 Critical soil water content, below which crop water stress occurs (m3/m3).

# Arguments
- `sat`: saturation (m3/m3)
- `pwp`: permanent wilting point (m3/m3)

# Keywords
- `p`: fraction of available water that can be depleted before crop water stress occurs
       (default: 0.5 for C3 plants in general)

# Returns
- `Float64`: critical soil water content (m3/m3)
"""
est_critical(sat, pwp, p=0.5) = pwp + p * (sat - pwp)


"""
   to_wc(vwc, thick)

 Convert volumetric soil water content m3/m3 to unit mm.

# Arguments
- `vwc`: volumetric soil water content (m3/m3)
- `thick`: soil layer thickness (m)

# Returns
- `Float64`: soil water content (mm)
"""
to_wc(vwc, thick) = vwc * thick * 1000


"""
   to_vwc(wc, thick)

 Convert soil water content from unit mm to volumetric unit m3/m3.

# Arguments
- `wc`: soil water content (mm)
- `thick`: soil layer thickness (m)

# Returns
- `Float64`: volumetric soil water content (m3/m3)
"""
to_vwc(wc, thick) = wc / (thick * 1000)


"""
    est_swc(c, s)

Soil water characteristics: saturation, field capacity, and permanent wiling point (m3/m3).

# Arguments
- `c`: clay (%)
- `s`: sand (%)

# Returns
- `Tuple{Float64, Float64, Float64}`:
        saturation, field capacity and permanent wiling point (m3/m3)
"""
function est_swc(c, s)
    s2 = s^2
    coef_a = exp(-4.396 - 0.0715 * c - 0.000488 * s2 - 0.00004285 * s2 * c)
    coef_b = -3.14 - 0.00222 * c^2 - 0.00003484 * s2 * c
    cmin = 1 - 0.007 * s    # min. clay content for min. SAT
    sat = 0.332 - 0.0007251 * s + 0.1276 * log10(max(cmin, c))   # SAT
    b = 1 / coef_b
    fc = (0.33333 / coef_a) ^ b    # FC
    pwp = (15 / coef_a) ^ b   # PWP
    sat, fc, pwp
end


"""
    ksat(sat, psd, airentry)

Saturated hydraulic conductivity (m/day).

# Arguments
- `sat`: saturation (m3/m3)
- `psd`: pore size distribution
- `airentry`: air entry suction (bubbling pressure) (kPa)

# Returns
- `Float64`: saturated hydraulic conductivity (m/day)
"""
ksat(sat, psd, airentry) = 864 * 0.07 * (sat * (1 - (airentry / 33) ^ psd)) ^ 4


"""
    heads_k(layer::SoilLayer)

Matric and gravity heads (m), and unsaturated hydraulic conductivity (m/day).

# Arguments
- `vwc`: volumetric soil water content (m3/m3)
- `sat`: saturation (m3/m3)
- `fc`: field capacity (m3/m3)
- `psd`: pore size distribution
- `airentry`: air entry suction (bubbling pressure) (kPa)
- `ksat`: saturated hydraulic conductivity (m/day)
- `depth`: soil layer depth (m)

# Returns
- `Tuple{Float64, Float64, Float64}`:
        matric and gravity heads (m), and unsaturated hydraulic conductivity (m/day)
"""
function heads_k(vwc, sat, fc, psd, airentry, ksat, depth)
    b = 1 / psd

    # matric suction (Saxton et al., 2006), convert from kPa to m by dividing by 10
    if vwc >= fc
        hm = (33 - (33 - airentry) * (vwc - fc) / (sat - fc)) / 10
    else
        a = exp(3.496508 + b * log(fc))
        hm = (a * max(0.05, vwc) ^ (-b)) / 10
    end

    matric = hm         # matric head (m)
    gravity = depth     # gravity head (m)

    # unsaturated hydraulic conductivity (m/day)
    ae = airentry / 10  # air entry (convert to m)
    k = ksat * ((matric > ae) ? (ae / hm) ^ (2 + 3 / b) : 1)

    matric, gravity, k
end


"""
    heads_k(layer::SoilLayer)

Matric and gravity heads (m), and unsaturated hydraulic conductivity (m/day).

# Arguments
- `layer::SoilLayer`: soil layer

# Notes
- Less verbose call to `heads_k(vwc, sat, fc, psd, airentry, ksat, depth)`,
  but ensure required properties (i.e., vwc, sat, fc, psd, airentry, ksat,
  and depth) in the `layer` object have been set first

# Returns
- `Tuple{Float64, Float64, Float64}`:
        matric and gravity heads (m), and unsaturated hydraulic conductivity (m/day)
"""
function heads_k(layer::SoilLayer)
    vwc = layer.vwc
    sat = layer.sat
    fc = layer.fc
    psd = layer.psd
    airentry = layer.airentry
    ksat = layer.ksat
    depth = layer.depth
    heads_k(vwc, sat, fc, psd, airentry, ksat, depth)
end


"""
    new_layer(thick, c, s, vwc, prevlayer::SoilLayer=SoilLayer())

Create and initialize a soil layer.

# Arguments
- `thick`: thickness of soil layer (m)
- `c`: clay in soil layer (%)
- `s`: sand in soil layer (%)
- `vwc`: volumetric soil water content of soil layer (m3/m3)

# Keywords
- `prevlayer::SoilLayer`: previous soil layer (default: `SoilLayer()`)

# Notes
- to create the first soil layer, pass a blank `SoilLayer` object
  for `prevlayer`, else `prevlayer` argument refers to the
  previous `SoilLayer` object

# Returns
- `SoilLayer`: a new soil layer
"""
function new_layer(thick, c, s, vwc, prevlayer::SoilLayer=SoilLayer())
    dg = to_dg(c, s)
    sat, fc, pwp = est_swc(c, s)
    psd = est_psd(dg)
    porosity = sat
    airentry = est_airentry(dg)

    # check for special code for given VWC (1=SAT, 2=FC, 3=PWP):
    if vwc < 0
        xs = [1.0, 2.0, 3.0]
        ys = [sat, fc, pwp]
        fn = LinearInterpolation(xs, ys, extrapolation_bc=Line())
        vwc = max(pwp, min(sat, fn(-vwc)))
    end

    wc = to_wc(vwc, thick)
    critical = est_critical(sat, pwp)
    ks = ksat(sat, psd, airentry)
    accthick = thick + prevlayer.accthick
    d = 0.5 * (prevlayer.thick + thick)
    depth = prevlayer.depth + d
    matric, gravity, k = heads_k(vwc, sat, fc, psd, airentry, ks, depth)
    tothead = matric + gravity

    SoilLayer(
        thick = thick,
        clay = c,
        sand = s,
        dg = dg,
        accthick = accthick,
        depth = depth,
        sat = sat,
        fc = fc,
        pwp = pwp,
        psd = psd,
        porosity = porosity,
        airentry = airentry,
        critical = critical,
        ksat = ks,
        vwc = vwc,
        wc = wc,
        k = k,
        matric = matric,
        gravity = gravity,
        tothead = tothead
    )
end


"""
    roots_water!(roots::RootsLayer, layers::AbstractVector{SoilLayer})

Set the soil water content (both in m3/m3 and mm) in the root zone.

# Arguments
- `roots::RootsLayer`: root zone soil layer
- `layers::AbstractVector{SoilLayer}`: all soil layers (with and without roots)

# Notes
- Update the argument `roots::RootsLayer` to hold the root zone water content.

# Returns
- nothing
"""
function roots_water!(roots::RootsLayer, layers::AbstractVector{SoilLayer})
    nlayers = roots.last_layer
    rootdepth = roots.depth
    wgt_wc = Vector{Float64}(undef, nlayers)

    for i ∈ 1:nlayers
        layer = layers[i]
        n = max(0, layer.accthick - rootdepth)
        w = max(0, layer.thick - n)
        wgt_wc[i] = layer.vwc * w
    end

    roots.vwc = sum(wgt_wc) / rootdepth
    roots.wc = to_wc(roots.vwc, rootdepth)
end


"""
    new_roots_layer(rootdepth, layers::AbstractVector{SoilLayer})

Create and initialize a soil layer that holds the entire root zone.

# Arguments
- `rootdepth`: rooting depth (m)
- `layers::AbstractVector{SoilLayer}`: all soil layers (with and without roots)

# Returns
- `RootsLayer`: soil layer that represents the entire root zone
"""
function new_roots_layer(rootdepth, layers::AbstractVector{SoilLayer})
    # weighted averages for soil layers encompassing the root zone:
    nlayers = length(layers)
    wgt_clay = Vector{Float64}(undef, nlayers)
    wgt_sand = Vector{Float64}(undef, nlayers)
    wgt_wcsat = Vector{Float64}(undef, nlayers)
    wgt_wcfc = Vector{Float64}(undef, nlayers)
    wgt_wcpwp = Vector{Float64}(undef, nlayers)
    wgt_w = Vector{Float64}(undef, nlayers)
    wgt_wc = Vector{Float64}(undef, nlayers)
    for i ∈ eachindex(layers)
        layer = layers[i]
        n = max(0, layer.accthick - rootdepth)
        w = max(0, layer.thick - n)   # weight
        wgt_clay[i] = layer.clay * w
        wgt_sand[i] = layer.sand * w
        wgt_wcsat[i] = layer.sat * w
        wgt_wcfc[i] = layer.fc * w
        wgt_wcpwp[i] = layer.pwp * w
        wgt_wc[i] = layer.vwc * w
        wgt_w[i] = w
    end

    # find the last soil layer that holds the roots:
    i = findfirst(isapprox(0), wgt_w)
    last_layer = isnothing(i) ? nlayers : (i - 1)

    # root zone soil texture (%):
    sumw = sum(wgt_w)
    c = sum(wgt_clay) / sumw
    s = sum(wgt_sand) / sumw

    # root zone water content (m3/m3):
    sat = sum(wgt_wcsat) / rootdepth
    fc = sum(wgt_wcfc) / rootdepth
    pwp = sum(wgt_wcpwp) / rootdepth
    vwc = sum(wgt_wc) / rootdepth
    wc = to_wc(vwc, rootdepth)
    critical = est_critical(sat, pwp)

    RootsLayer(
        depth = rootdepth,
        clay = c,
        sand = s,
        vwc = vwc,
        wc = wc,
        sat = sat,
        fc = fc,
        pwp = pwp,
        critical = critical,
        last_layer = last_layer
    )
end
