"""
    clone(src::OP)

Clone by deep copying the given `OP` object source.

# Arguments
- `src::OP`: source to be cloned

# Notes
- Faster object deep copy than using Julia's `deepcopy`.

# Returns
- `OP`: cloned `OP` object
"""
function clone(src::OP)
    sp = src.parts
    pinnae = sp.pinnae
    trunk = sp.trunk
    rachis = sp.rachis
    roots = sp.roots
    maleflo = sp.maleflo
    femflo = sp.femaflo
    bunches = sp.bunches

    parts = Parts(
        pinnae = Properties(
            maint = pinnae.maint,
            frac = pinnae.frac,
            growth = pinnae.growth,
            death = pinnae.death,
            weight = pinnae.weight
        ),
        rachis = Properties(
            maint = rachis.maint,
            frac = rachis.frac,
            growth = rachis.growth,
            death = rachis.death,
            weight = rachis.weight
        ),
        trunk = Properties(
            maint = trunk.maint,
            frac = trunk.frac,
            growth = trunk.growth,
            death = trunk.death,
            weight = trunk.weight
        ),
        roots = Properties(
            maint = roots.maint,
            frac = roots.frac,
            growth = roots.growth,
            death = roots.death,
            weight = roots.weight
        ),
        maleflo = Properties(
            maint = maleflo.maint,
            frac = maleflo.frac,
            growth = maleflo.growth,
            death = maleflo.death,
            weight = maleflo.weight
        ),
        femaflo = Properties(
            maint = femflo.maint,
            frac = femflo.frac,
            growth = femflo.growth,
            death = femflo.death,
            weight = femflo.weight
        ),
        bunches = Properties(
            maint = bunches.maint,
            frac = bunches.frac,
            growth = bunches.growth,
            death = bunches.death,
            weight = bunches.weight
        )
    )

    rlayer = src.roots
    roots = RootsLayer(
        depth = rlayer.depth,
        clay = rlayer.clay,
        sand = rlayer.sand,
        vwc = rlayer.vwc,
        wc = rlayer.wc,
        sat = rlayer.sat,
        fc = rlayer.fc,
        pwp = rlayer.pwp,
        critical = rlayer.critical,
        last_layer = rlayer.last_layer
    )

    et = HeatFluxes(
        total = src.et.total,
        crop = src.et.crop,
        soil = src.et.soil
    )

    h = HeatFluxes(
        total = src.h.total,
        crop = src.h.crop,
        soil = src.h.soil
    )

    layers = map(l -> SoilLayer(
        thick = l.thick,
        clay = l.clay,
        sand = l.sand,
        dg = l.dg,
        accthick = l.accthick,
        depth = l.depth,
        sat = l.sat,
        fc = l.fc,
        pwp = l.pwp,
        psd = l.psd,
        porosity = l.porosity,
        airentry = l.airentry,
        critical = l.critical,
        ksat = l.ksat,
        vwc = l.vwc,
        wc = l.wc,
        k = l.k,
        matric = l.matric,
        gravity = l.gravity,
        tothead = l.tothead,
        t = l.t,
        e = l.e,
        influx = l.influx,
        outflux = l.outflux,
        netflux = l.netflux
        ),
        src.layers
    )

    OP(
        yap = src.yap,
        date = src.date,
        doy = src.doy,
        solardec = src.solardec,
        sunrise = src.sunrise,
        sunset = src.sunset,
        etrad = src.etrad,
        totrad = src.totrad,
        drrad = src.drrad,
        dfrad = src.dfrad,
        tmin = src.tmin,
        tmax = src.tmax,
        wind = src.wind,
        rain = src.rain,
        treeage = src.treeage,
        plantdens = src.plantdens,
        parts = parts,
        treehgt = src.treehgt,
        trunkhgt = src.trunkhgt,
        assim4maint = src.assim4maint,
        assim4growth = src.assim4growth,
        assim4gen = src.assim4gen,
        vdmreq = src.vdmreq,
        tdmwgt = src.tdmwgt,
        vdmwgt = src.vdmwgt,
        flowgt = src.flowgt,
        frdwgt = src.frdwgt,
        vdmgro = src.vdmgro,
        tdmgro = src.tdmgro,
        yield = src.yield,
        sla = src.sla,
        lai = src.lai,
        laimax = src.laimax,
        vdmmax = src.vdmmax,
        layers = layers,
        roots = roots,
        stress_t = src.stress_t,
        stress_e = src.stress_e,
        netrain = src.netrain,
        aet_soil = src.aet_soil,
        aet_crop = src.aet_crop,
        dp = src.dp,
        refhgt = src.refhgt,
        d = src.d,
        z0 = src.z0,
        kwind = src.kwind,
        keddy = src.keddy,
        leaflength = src.leaflength,
        leafwidth = src.leafwidth,
        efflai = src.efflai,
        assimilates = src.assimilates,
        et = et,
        h = h
    )
end
