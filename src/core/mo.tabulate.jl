"""
    tabulate(res::AbstractVector{OP})

Convert daily model results, stored as an array of `OP`, into a `DataFrame`.

# Arguments
- `res::AbstractVector{OP}`: daily model results

# Returns
- `DataFrame`: daily model results stored as a `DataFrame`
"""
function tabulate(res::AbstractVector{OP})
    gp = getproperty
    gpfn(r) = (ar -> gp.(r, ar))    # function object

    gp_res = gpfn(res)   # iterator for each day model results
    gp_parts = gp_res(:parts)   # iterator for every tree part
    gp_pinnae = gp.(gp_parts, :pinnae) |> gpfn   # iterator for pinnae properties
    gp_rachis = gp.(gp_parts, :rachis) |> gpfn   # iterator for rachis
    gp_trunk = gp.(gp_parts, :trunk) |> gpfn   # iterator for trunk
    gp_roots = gp.(gp_parts, :roots) |> gpfn   # iterator for roots
    gp_maleflo = gp.(gp_parts, :maleflo) |> gpfn   # iterator for male flowers
    gp_femaflo = gp.(gp_parts, :femaflo) |> gpfn   # iterator for female flower
    gp_bunches = gp.(gp_parts, :bunches) |> gpfn   # iterator for bunches
    gp_rootzone = gp_res(:roots) |> gpfn   # iterator for root zone soil layer
    gp_et = gp_res(:et) |> gpfn   # iterator for latent heat heat flux
    gp_h = gp_res(:h) |> gpfn   # iterator for sensible heat flux

    df = DataFrame(
        yap = gp_res(:yap),
        date = gp_res(:date),
        doy = gp_res(:doy),
        solardec = gp_res(:solardec),
        sunrise = gp_res(:sunrise),
        sunset = gp_res(:sunset),
        etrad = gp_res(:etrad),
        totrad = gp_res(:totrad),
        drrad = gp_res(:drrad),
        dfrad = gp_res(:dfrad),
        tmin = gp_res(:tmin),
        tmax = gp_res(:tmax),
        wind = gp_res(:wind),
        rain = gp_res(:rain),
        treeage = gp_res(:treeage),
        plantdens = gp_res(:plantdens),
        treehgt = gp_res(:treehgt),
        trunkhgt = gp_res(:trunkhgt),
        assim4maint = gp_res(:assim4maint),
        assim4growth = gp_res(:assim4growth),
        assim4gen = gp_res(:assim4gen),
        vdmreq = gp_res(:vdmreq),
        tdmwgt = gp_res(:tdmwgt),
        vdmwgt = gp_res(:vdmwgt),
        flowgt = gp_res(:flowgt),
        frdwgt = gp_res(:frdwgt),
        vdmgro = gp_res(:vdmgro),
        tdmgro = gp_res(:tdmgro),
        yield = gp_res(:yield),
        sla = gp_res(:sla),
        lai = gp_res(:lai),
        laimax = gp_res(:laimax),
        vdmmax = gp_res(:vdmmax),

        stress_t = gp_res(:stress_t),
        stress_e = gp_res(:stress_e),
        netrain = gp_res(:netrain),
        aet_soil = gp_res(:aet_soil),
        aet_crop = gp_res(:aet_crop),
        dp = gp_res(:dp),

        refhgt = gp_res(:refhgt),
        d = gp_res(:d),
        z0 = gp_res(:z0),
        kwind = gp_res(:kwind),
        keddy = gp_res(:keddy),
        leaflength = gp_res(:leaflength),
        leafwidth = gp_res(:leafwidth),
        efflai = gp_res(:efflai),
        assimilates = gp_res(:assimilates),

        roots_depth = gp_rootzone(:depth),
        roots_clay = gp_rootzone(:clay),
        roots_sand = gp_rootzone(:sand),
        roots_vwc = gp_rootzone(:vwc),
        roots_wc = gp_rootzone(:wc),
        roots_sat = gp_rootzone(:sat),
        roots_fc = gp_rootzone(:fc),
        roots_pwp = gp_rootzone(:pwp),
        roots_critical = gp_rootzone(:critical),
        roots_last_layer = gp_rootzone(:last_layer),

        et_total = gp_et(:total),
        et_crop = gp_et(:crop),
        et_soil = gp_et(:soil),
        h_total = gp_h(:total),
        h_crop = gp_h(:crop),
        h_soil = gp_h(:soil),

        parts_pinnae_maint = gp_pinnae(:maint),
        parts_pinnae_frac = gp_pinnae(:frac),
        parts_pinnae_growth = gp_pinnae(:growth),
        parts_pinnae_death = gp_pinnae(:death),
        parts_pinnae_weight = gp_pinnae(:weight),

        parts_rachis_maint = gp_rachis(:maint),
        parts_rachis_frac = gp_rachis(:frac),
        parts_rachis_growth = gp_rachis(:growth),
        parts_rachis_death = gp_rachis(:death),
        parts_rachis_weight = gp_rachis(:weight),

        parts_trunk_maint = gp_trunk(:maint),
        parts_trunk_frac = gp_trunk(:frac),
        parts_trunk_growth = gp_trunk(:growth),
        parts_trunk_death = gp_trunk(:death),
        parts_trunk_weight = gp_trunk(:weight),

        parts_roots_maint = gp_roots(:maint),
        parts_roots_frac = gp_roots(:frac),
        parts_roots_growth = gp_roots(:growth),
        parts_roots_death = gp_roots(:death),
        parts_roots_weight = gp_roots(:weight),

        parts_maleflo_maint = gp_maleflo(:maint),
        parts_maleflo_frac = gp_maleflo(:frac),
        parts_maleflo_growth = gp_maleflo(:growth),
        parts_maleflo_death = gp_maleflo(:death),
        parts_maleflo_weight = gp_maleflo(:weight),

        parts_femaflo_maint = gp_femaflo(:maint),
        parts_femaflo_frac = gp_femaflo(:frac),
        parts_femaflo_growth = gp_femaflo(:growth),
        parts_femaflo_death = gp_femaflo(:death),
        parts_femaflo_weight = gp_femaflo(:weight),

        parts_bunches_maint = gp_bunches(:maint),
        parts_bunches_frac = gp_bunches(:frac),
        parts_bunches_growth = gp_bunches(:growth),
        parts_bunches_death = gp_bunches(:death),
        parts_bunches_weight = gp_bunches(:weight),

    )

    fields = fieldnames(SoilLayer)
    nlayers = length(res[1].layers)
    layers = map(i -> getindex.(gp_res(:layers), i), 1:nlayers)
    for i ∈ eachindex(layers)
        pfx = "layers$(i)_"
        lyr = layers[i]
        for f ∈ fields
            df[!, pfx * String(f)] = gp.(lyr, f)
        end
    end
    df
end


"""
    tabulate(res::AbstractVector{OPHour})

Convert hourly model results, stored as an array of `OPHour`, into a `DataFrame`.

# Arguments
- `res::AbstractVector{OPHour}`: hourly model results

# Returns
- `DataFrame`: hourly model results stored as a `DataFrame`
"""
function tabulate(res::AbstractVector{OPHour})
    gp = getproperty
    gpfn(r) = (ar -> gp.(r, ar))    # function object

    gp_res = gpfn(res)   # iterator for each hour model results
    gp_par = gp_res(:par) |> gpfn
    gp_assimcoef = gp_res(:assimcoef) |> gpfn
    gp_leafassim = gp_res(:leafassim) |> gpfn
    gp_A = gp_res(:A) |> gpfn
    gp_stresses = gp_res(:stresses) |> gpfn
    gp_R = gp_res(:R) |> gpfn
    gp_et = gp_res(:et) |> gpfn   # iterator for latent heat heat flux
    gp_h = gp_res(:h) |> gpfn     # iterator for sensible heat flux

    DataFrame(
        date = gp_res(:date),
        doy = gp_res(:doy),
        solarhour = gp_res(:solarhour),
        air_temp = gp_res(:air_temp),
        dew_temp = gp_res(:dew_temp),
        svp = gp_res(:svp),
        vp = gp_res(:vp),
        rh = gp_res(:rh),
        vpd = gp_res(:vpd),
        vpd0 = gp_res(:vpd0),
        svp_slope = gp_res(:svp_slope),
        solarinc = gp_res(:solarinc),
        solarhgt = gp_res(:solarhgt),
        solarazi = gp_res(:solarazi),
        solarcon = gp_res(:solarcon),
        etrad = gp_res(:etrad),
        totrad = gp_res(:totrad),
        drrad = gp_res(:drrad),
        dfrad = gp_res(:dfrad),
        netrad = gp_res(:netrad),
        wind = gp_res(:wind),

        kdr = gp_res(:kdr),
        kdf = gp_res(:kdf),
        clump = gp_res(:clump),
        gap = gp_res(:gap),
        lsl = gp_res(:lsl),
        lsh = gp_res(:lsh),

        par_outdr = gp_par(:outdr),
        par_outdf = gp_par(:outdf),
        par_indrscatter = gp_par(:indrscatter),
        par_indr = gp_par(:indr),
        par_inscatter = gp_par(:inscatter),
        par_indf = gp_par(:indf),
        par_abssunlit = gp_par(:abssunlit),
        par_absshaded = gp_par(:absshaded),

        co2_ambient = gp_res(:co2_ambient),
        co2_internal = gp_res(:co2_internal),

        assimcoef_mmco2 = gp_assimcoef(:mmco2),
        assimcoef_mmo2 = gp_assimcoef(:mmo2),
        assimcoef_specificity = gp_assimcoef(:specificity),
        assimcoef_vcmax = gp_assimcoef(:vcmax),
        assimcoef_co2_pt = gp_assimcoef(:co2_pt),

        leafassim_vc = gp_leafassim(:vc),
        leafassim_vqsl = gp_leafassim(:vqsl),
        leafassim_vqsh = gp_leafassim(:vqsh),
        leafassim_vs = gp_leafassim(:vs),
        leafassim_sunlit = gp_leafassim(:sunlit),
        leafassim_shaded = gp_leafassim(:shaded),

        canopyassim = gp_res(:canopyassim),

        A_total = gp_A(:total),
        A_crop = gp_A(:crop),
        A_soil = gp_A(:soil),
        A_net = gp_A(:net),
        A_g = gp_A(:g),

        ustar = gp_res(:ustar),
        utreehgt = gp_res(:utreehgt),
        canopy_temp = gp_res(:canopy_temp),

        stresses_water = gp_stresses(:water),
        stresses_vpd = gp_stresses(:vpd),
        stresses_par = gp_stresses(:par),

        R_rsa = gp_R(:rsa),
        R_raa = gp_R(:raa),
        R_rca = gp_R(:rca),
        R_rst = gp_R(:rst),
        R_rcs = gp_R(:rcs),
        R_rss = gp_R(:rss),

        et_total = gp_et(:total),
        et_crop = gp_et(:crop),
        et_soil = gp_et(:soil),
        h_total = gp_h(:total),
        h_crop = gp_h(:crop),
        h_soil = gp_h(:soil),
    )
end
