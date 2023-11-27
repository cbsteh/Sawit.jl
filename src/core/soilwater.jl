"""
    Fluxes

Soil water fluxes (m/day).

# Fields
- `t::Float64`: plant water uptake (m/day)
- `e::Float64`: soil evaporation (m/day)
- `influx::Float64`: water influx (m/day)
- `outflux::Float64`: water outflux (m/day)
- `netflux::Float64`: water netflux (m/day)
"""
@with_kw mutable struct Fluxes
    t::Float64 = 0.0
    e::Float64 = 0.0
    influx::Float64 = 0.0
    outflux::Float64 = 0.0
    netflux::Float64 = 0.0
end


"""
    copy!(flx::Fluxes, layer::SoilLayer)

Copy all fields from `Fluxes` to their respective counterparts in `SoilLayer`.
Internal use.

# Arguments
- `flx::Fluxes`: water fluxes
- `layer::SoilLayer`: soil layer

# Returns
- nothing
"""
function copy!(flx::Fluxes, layer::SoilLayer)
    layer.t = flx.t
    layer.e = flx.e
    layer.influx = flx.influx
    layer.outflux = flx.outflux
    layer.netflux = flx.netflux
end


"""
    Transient

Store intermediary flux calculations. Internal use.

# Fields
- `prvpsi::Float64`: intermediary variable to determine plant water uptake
- `pf::Vector{Fluxes}`: intermediary fluxes
- `cf::Vector{Fluxes}`: cumulative fluxes
"""
mutable struct Transient
    prvpsi::Float64
    pf::Vector{Fluxes}
    cf::Vector{Fluxes}

    function Transient(sz::Int)
        pf = [Fluxes() for _ ∈ 1:sz]
        cf = [Fluxes() for _ ∈ 1:sz]
        new(0.0, pf, cf)
    end
end


"""
    Base.:+(x::Fluxes, y::Fluxes)

Addition of two `Fluxes`: x + y

# Arguments
- `x::Fluxes`: lefthand `Fluxes`
- `y::Fluxes`: righthand `Fluxes`

# Notes
- Mostly for internal use.

# Returns
- `Fluxes`: sum of two `Fluxes`
"""
function Base.:+(x::Fluxes, y::Fluxes)
    Fluxes(
        t = x.t + y.t,
        e = x.e + y.e,
        influx = x.influx + y.influx,
        outflux = x.outflux + y.outflux,
        netflux = x.netflux + y.netflux,
    )
end


"""
    Base.:/(x::Fluxes, y)

Division of a `Fluxes` object, x, by a value, y: x / y.

# Arguments
- `x::Fluxes`: object to be divided
- `y`: value

# Notes
- Mostly for internal use.

# Returns
- `Fluxes`: resultant division of a `Fluxes` by `y`
"""
function Base.:/(x::Fluxes, y)
    Fluxes(
        t = x.t / y,
        e = x.e / y,
        influx = x.influx / y,
        outflux = x.outflux / y,
        netflux = x.netflux / y,
    )
end


"""
    net_rainfall(lai, rain)

Net rainfall (below canopies) (mm/day).

# Arguments
- `lai`: leaf area index (m2/m2)
- `rain`: daily gross rain (above canopies) (mm/day)

# Returns
- `Float64`: Net rainfall (mm/day)
"""
net_rainfall(lai, rain) = min(0.7295, 1 - 0.0541 * lai) * rain


"""
    new_rooting_depth(roots::RootsLayer, stress_t, maxdepth)

Update the rooting depth (m).

# Arguments
- `roots::RootsLayer`: root zone layer
- `stress_t`: crop water stress (fraction, 0-1)
- `maxdepth`: maximum depth the rooting zone cannot exceed (m)
- `maxrate`: maximum rate of increase in rooting depth (m)

# Returns
- `Float64`: new, increased rooting depth (m)
"""
function new_rooting_depth(roots::RootsLayer, stress_t, maxdepth, maxrate)
    newdepth = roots.depth + maxrate * stress_t   # root growth is limited by water stress
    min(maxdepth, newdepth)   # do not exceed max. depth
end


"""
    et_reduction(roots::RootsLayer, firstlayer::SoilLayer)

Fractional reduction in evaporation and transpiration (0-1).

# Arguments
- `roots::RootsLayer`: root zone layer
- `firstlayer::SoilLayer`: first (top) soil layer

# Notes
- 1 = no reduction, 0 = max. reduction

# Returns
- `Tuple{Float64, Float64}`:
        fractional reduction in potential evaporation and transpiration (0-1)
"""
function et_reduction(roots::RootsLayer, firstlayer::SoilLayer)
    rde = 1 / (1 + (3.6073 * (firstlayer.vwc / firstlayer.sat)) ^ -9.3172)

    vwc = roots.vwc
    cr = roots.critical
    pwp = roots.pwp
    rdt = (vwc >= cr) ? 1.0 : # no stress because water content >= critical point
          (pwp < vwc < cr) ? (vwc - pwp) / (cr - pwp) :  # linear increase in stress
          0.0   # max. stress because water content <= PWP

    rdt, rde
    # 1.0, 1.0        # CBSTEH: No stresses to simulate no water deficits
end


"""
    actual_et(stress_t, stress_e, pet_t, pet_e)

Actual soil evaporation and traspiration (mm/day).

# Arguments
- `stress_t`: reduction to plant transpiration (fraction, 0-1)
- `stress_e`: reduction to soil evaporation (fraction, 0-1)
- `pet_t`: potential transpiration (mm/day)
- `pet_t`: potential soil evaporation (mm/day)

# Returns
- `Tuple{Float64, Float64}`:
        actual soil evaporation and transpiration (mm/day)
"""
function actual_et(stress_t, stress_e, pet_t, pet_e)
    aet_t = stress_t * pet_t
    aet_e = stress_e * pet_e
    aet_t, aet_e
end


"""
    influx_from_deeplayer(lastlayer::SoilLayer, has_watertable=false)

Influx of water from groundwater or deeper soil layers (m/day).

# Arguments
- `lastlayer::SoilLayer`: last soil layer
- `has_watertable`: `true` if groundwater is present, else `false` (default false)

# Returns
- `Float64`: influx of water from groundwater or deeper soil layers (m/day)
"""
function influx_from_deeplayer(lastlayer::SoilLayer, has_watertable=false)
    @unpack sat, fc, psd, airentry, ksat, depth = lastlayer
    vwc = has_watertable ? sat : fc

    deep_hm, deep_hg, deep_k = heads_k(vwc, sat, fc, psd, airentry, ksat, depth)
    deep_hg = lastlayer.accthick
    deep_tothead = deep_hm + deep_hg

    last_k = lastlayer.k
    last_tothead = lastlayer.matric + lastlayer.gravity

    n = log(deep_k) - log(last_k)
    kmean = isapprox(n, 0) ? deep_k : (deep_k - last_k) / n
    thick = lastlayer.thick * 0.5
    kmean * (deep_tothead - last_tothead) / thick
end


"""
    _influx!(flx::Transient, mo::Output, layer_i::Int)

Water entry (influx) into a given soil layer.

# Arguments
- `flx::Transient`: store intermediary fluxes
- `mo::Output`: model output results
- `layer_i::Int`: soil layer number

# Notes
- Internal use. Updates `flx::Transient`.

# Returns
- nothing
"""
function _influx!(flx::Transient, mo::Output, layer_i::Int)
    op = mo.op
    cur = op.layers[layer_i]
    prv = (layer_i > 1) ? op.layers[layer_i - 1] : nothing
    not_topsoil = !isnothing(prv)

    ei = not_topsoil ? 0.0 : op.aet_soil / 1000    # actual E (m/day)
    cj = min(1, cur.accthick / op.roots.depth)
    curpsi = 1.8 * cj - 0.8 * cj^2
    ti = op.aet_crop * (curpsi - flx.prvpsi) / 1000   # actual T (m/day)
    flx.prvpsi = curpsi

    curpf = flx.pf[layer_i]
    curpf.t = ti
    curpf.e = ei

    if not_topsoil
        # use Darcy's law for second soil layer onwards
        n = log(cur.k) - log(prv.k)
        k = isapprox(n, 0) ? cur.k : (cur.k - prv.k) / n
        delta = cur.tothead - prv.tothead
        curpf.influx = k * delta / (cur.depth - prv.depth)
    else
        # first layer influx is simply the net rainfall
        netrain = op.netrain / 1000       # net rainfall (convert to m)
        curpf.influx = netrain
    end
end


"""
    _outflux!(flx::Transient, mo::Output, mi::Input, layer_i::Int)

Water exit (outflux) out of a given soil layer.

# Arguments
- `flx::Transient`: store intermediary fluxes
- `mo::Output`: model output results
- `mi::Input`: model input
- `layer_i::Int`: soil layer number

# Notes
- Internal use.
- Update the arguments `flx::Transient` and `out::OutputWater`.

# Returns
- nothing
"""
function _outflux!(flx::Transient, mo::Output, mi::Input, layer_i::Int)
    op = mo.op
    cur = op.layers[layer_i]
    curpf = flx.pf[layer_i]
    nlayers = length(op.layers)
    nxt = (layer_i < nlayers) ? op.layers[layer_i + 1] : nothing
    not_lastlayer = !isnothing(nxt)

    wc = cur.vwc * cur.thick  # current water content (m)
    if not_lastlayer
        outflux = flx.pf[layer_i + 1].influx   # outflux is the next soil layer's influx
    else
        outflux = influx_from_deeplayer(op.layers[end], mi.has_watertable)
    end

    # ensure a soil layer cannot be too dry (>0.005 m3/m3)
    influx = curpf.influx  # water into current layer
    nextwc = influx + wc - outflux
    drylmt = cur.thick * 0.005
    if nextwc < drylmt
        outflux = influx + wc - drylmt
    end

    if not_lastlayer
        flx.pf[layer_i + 1].influx = outflux
    end

    # determine the net fluxes and water content for current layer:
    curpf.outflux = outflux
    curpf.netflux = curpf.influx - curpf.outflux - curpf.e - curpf.t

    # update at every sub-interval step:
    wc += curpf.netflux / mi.nintervals
    cur.vwc = max(0.005, wc / cur.thick)  # vol. water content (m3/m3)
    cur.wc = cur.vwc * cur.thick * 1000  # water content (mm)

    cur.matric, cur.gravity, cur.k = heads_k(cur)
    cur.tothead = cur.matric + cur.gravity

    # sum the water fluxes in every sub-interval step to determine final fluxes
    #    at the end of a daily time step
    flx.cf[layer_i] = flx.cf[layer_i] + (curpf / mi.nintervals)
end


"""
    _fluxes!(flx::Transient, mo::Output, mi::Input)

Iteratively determine the various water fluxes for every soil layer.

# Arguments
- `flx::Transient`: store intermediary fluxes
- `mo::Output`: model output results
- `mi::Input`: model input

# Notes
- Internal use. Update the arguments `flx::Transient` and `mo::OutputWater`

# Returns
- nothing
"""
function _fluxes!(flx::Transient, mo::Output, mi::Input)
    op = mo.op
    @unpack roots, layers, et = op

    op.stress_t, op.stress_e = et_reduction(roots, layers[1])
    op.aet_crop, op.aet_soil = actual_et(op.stress_t, op.stress_e, et.crop, et.soil)

    # 1. calculate every influx:
    flx.prvpsi = 0.0
    foreach(layer_i -> _influx!(flx, mo, layer_i), eachindex(layers))
    # 2. calculate every outflux, netflux, then water content
    foreach(layer_i -> _outflux!(flx, mo, mi, layer_i), eachindex(layers))
    # 3. update the root zone water content
    roots_water!(roots, layers)
end


"""
    update_soilwater!(mo::Output, mi::Input)

Update all soil layer properties and root zone layer.

# Arguments
- `mo::Output`: model output results
- `mi::Input`: model input

# Notes
- Update the argument `out::OutputWater`.

# Returns
- nothing
"""
function update_soilwater!(mo::Output, mi::Input)
    op = mo.op
    @unpack lai, rain, layers, roots, stress_t = op

    op.netrain = net_rainfall(lai, rain)

    maxdepth = layers[end].accthick     # max. rooting depth
    maxrate = mi.maxrate_rootdepth      # max. rate of rooting depth increase
    rootdepth = new_rooting_depth(op.roots, stress_t, maxdepth, maxrate)
    if rootdepth > roots.depth
        # new `RootsLayer` object returned, so set `op.roots` to this new
        # object and, most importantly, reinitialize `roots`, else it
        # still points to the old `RootsLayer` object
        op.roots = roots = new_roots_layer(rootdepth, layers)
    end

    # solve the water balance:
    flx = Transient(length(layers))
    foreach(i -> _fluxes!(flx, mo, mi), 1:mi.nintervals)
    copy!.(flx.cf, layers)

    # deep percolation:
    last = layers[roots.last_layer]
    d = last.accthick - roots.depth
    m = d / last.thick
    lastnetflux = last.netflux
    lastoutflux = last.outflux
    op.dp = (lastoutflux + m * lastnetflux) * 1000  # mm/day
end


"""
    init_soilwater!(mo::Output, mi::Input)

Initialize the soil water parameters.

# Arguments
- `mo::Output`: model output soil results
- `mi::Input`: model input soil parameters

# Returns
- nothing
"""
function init_soilwater!(mo::Output, mi::Input)
    layers = Vector{SoilLayer}(undef, mi.nlayers)
    for i ∈ eachindex(layers)
        prevlayer = (i > 1) ? layers[i-1] : SoilLayer()
        @unpack thick, clay, sand, vwc = mi.layers[i]
        layers[i] = new_layer(thick, clay, sand, vwc, prevlayer)
    end

    op = mo.op
    op.layers = layers
    op.roots = new_roots_layer(mi.rootdepth, op.layers)
    op.stress_t, op.stress_e = et_reduction(op.roots, op.layers[1])
end
