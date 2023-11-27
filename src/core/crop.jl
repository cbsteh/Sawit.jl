"""
    dm_weights(parts::Parts)

Sum the tree parts weights for TDM, VDM, flowers, and fronds (kg DM/palm).

# Arguments
- `parts::Parts`: all tree parts

# Returns
- `Tuple{Float64, Float64, Float64, Float64}`:
        total, vegetative, flower, and fronds weights (kg DM/palm)
"""
function dm_weights(parts::Parts)
    flowgt = parts.maleflo.weight + parts.femaflo.weight  # flowers
    frdwgt = parts.pinnae.weight + parts.rachis.weight    # fronds
    vdmwgt = sum(vegparts(parts, :weight))
    tdmwgt = vdmwgt + flowgt + parts.bunches.weight
    tdmwgt, vdmwgt, flowgt, frdwgt
end


"""
    dm_growth(parts::Parts, yield)

Sum the growth rates for the aboveground tree parts (kg DM/palm/day).

# Arguments
- `parts::MO,.Parts`: all tree parts
- `yield`: FFB yield (kg DM/palm/day)

# Returns
- `Tuple{Float64, Float64}`: aboveground total and vegetative growth rates (kg DM/palm/day)
"""
function dm_growth(parts::Parts, yield)
    vdm = parts.pinnae.growth + parts.rachis.growth + parts.trunk.growth
    tdm = vdm + yield
    tdm, vdm
end


"""
    canopy_height(treeage::Int)

Canopy height (m).

# Arguments
- `treeage::Int`: tree age (days)

# Returns
- `Float64`: canopy height (m)
"""
canopy_height(treeage::Int) = (0.1382 * treeage + 150.91) / 100


"""
    trunk_height(treeage::Int, pd::Int)

Trunk height based on tree age and planting density (m).

# Arguments
- `treeage::Int`: tree age (days)
- `pd::Int`: planting density (palms/ha)

# Returns
- `Float64`: trunk height (m)
"""
function trunk_height(treeage::Int, pd::Int)
    a, b, c = 4.7862, -4741.6896, -196.656
    exp(a + b / pd^2 + c / sqrt(treeage))
end


"""
    tree_height(treeage::Int, pd::Int, cr)

Tree (trunk + canopy) height and trunk height (m).

# Arguments
- `treeage::Int`: tree age (days)
- `pd::Int`: planting density (palms/ha)
- `sand`: sand content in the root zone (%)

# Notes
- This method determines trunk height differently than
  `trunk_height(treeage::Int, pd::Int)`. Here, trunk height is
   additionally a function of the root zone's sand content.

# Returns
- `Tuple{Float64, Float64}`: total tree height and trunk height (m)
"""
function tree_height(treeage::Int, pd::Int, sand)
    # linear increase for trunk height growth rate, but max at > 11 yrs:
    a, b, c = 0.274424613, 0.000306785, 0.002212069
    dhmax = a + b * pd + c * sand   # max. height growth rate
    yrs = treeage / 365     # tree age in years
    dh = (yrs < 11) ? 0.1 * dhmax * (yrs - 1) : dhmax
    trunkhgt = dh * yrs
    treehgt = trunkhgt + canopy_height(treeage)
    treehgt, trunkhgt
end


"""
    leaf_dimension(treeage)

Mean length and width of leaflets (pinnae) (m).

# Arguments:
- `treeage`: tree age (days)

# Returns
- `Tuple{Float64, Float64}`: mean length and width of leaflets (pinnae) (m)
"""
function leaf_dimension(treeage)
    treeage /= 365  # tree age in years
    leaflen = 0.2191 * log(treeage) + 0.475   # leaflet mean length (m)
    leafwid = 0.0152 * log(treeage) + 0.0165  # leaflet mean width (m)
    leaflen, leafwid
end


"""
    sla_lai(treeage, pd, pinnae_wgt, sla_table::Afgen)

Specific leaf area (m2 leaf/kg DM leaf) and leaf area index (m2 leaf/m2 ground).

# Arguments
- `treeage`: tree age (days)
- `pd`: planting density (palms/ha)
- `pinnae_wgt`: pinnae weight (kg DM/palm)
- `sl_table::Afgen`: table of tree age vs specific leaf area (m2 leaf/kg DM leaf)

# Returns
- `Tuple{Float64, Float64}`: SLA (m2 leaf/kg DM leaf) and LAI (m2 leaf/m2 ground).
"""
function sla_lai(treeage, pd, pinnae_wgt, sla_table::Afgen)
    sla = sla_table(treeage)
    lai = pinnae_wgt * sla * pd / 10000
    sla, lai
end


"""
    lai_maximum(pd::Int)

Maximum leaf area index (LAI) (m2 leaf/m2 ground).

# Arguments
- `pd::Int`: planting density (palms/ha)

# Returns
- `Float64`: maximum LAI (m2 leaf/m2 ground)
"""
lai_maximum(pd::Int) = 1.4119 + 0.0267 * pd


"""
    vdm_maximum(pd::Int)

Maximum annual production rate of vegetative dry matter (VDM) for a given
planting density (kg DM/palm/year).

# Arguments
- `pd::Int`: planting density (palms/ha)

# Returns
- `Float64`: maximum VDM annual production rate (kg DM/palm/year)
"""
vdm_maximum(pd::Int) = 160.06 - 0.2416 * pd


"""
    maintenance_respiration(mo::Output, mi::Input)

Maintenance respiration for tree parts (kg CH2O/palm/day).

# Arguments
- `mo::Output`: model output crop growth results
- `mi::Input`: model input crop parameters

# Returns
- `Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}`:
        maintenance respiration (kg CH2O/palm/day) for every tree part (pinnae,
        rachis, trunk, roots, male flowers, female flowers, and bunches). Total
        maintenance is also returned as the eigth (last) element.
"""
function maintenance_respiration(mo::Output, mi::Input)
    op = mo.op
    tmean = 0.5 * (op.tmin + op.tmax)
    treeage = op.treeage

    # temperature correction for maintenance:
    temp_corr(val25, q10=2) = val25 * q10^((tmean - 25) / 25)

    function maintcoef(n_table, m_table, xc=2.0)
        # tree part maintenance coefficient (kg CH2O/kg DM)
        # xc - correction factor to include other nutrients (not P, K, Ca, and Mg)
        val25 = 0.04 * 6.25 * n_table(treeage) + 0.08 * xc * m_table(treeage)
        temp_corr(val25)
    end

    parts = op.parts

    # pinnae:
    part = parts.pinnae
    daylen = op.sunset - op.sunrise
    mc_pinnae = maintcoef(mi.pinnae_n_table, mi.pinnae_m_table)
    m_pinnae = part.weight * mc_pinnae * (24 - daylen) / 24

    # rachis:
    part = parts.rachis
    mc_rachis = maintcoef(mi.rachis_n_table, mi.rachis_m_table)
    m_rachis = part.weight * mc_rachis

    # trunk:
    part = parts.trunk
    mc_trunk = maintcoef(mi.trunk_n_table, mi.trunk_m_table)
    upper_wgt = min(45, part.weight)
    lower_wgt = part.weight - upper_wgt
    m_trunk = upper_wgt * mc_trunk + lower_wgt * mc_trunk * 0.06

    # roots:
    part = parts.roots
    mc_roots = maintcoef(mi.roots_n_table, mi.roots_m_table)
    m_roots = part.weight * mc_roots

    # male flowers:
    part = parts.maleflo
    m_maleflo = part.weight * mc_rachis

    # female flowers (immature bunches):
    part = parts.femaflo
    m_femaflo = part.weight * mc_rachis

    # mature bunches:
    part = parts.bunches
    m_bunches = part.weight * temp_corr(0.0027)

    # total maintenance:
    if 15 < tmean < 46
        tdmwgt = sum(allparts(parts, :weight))
        m_metabolic = temp_corr(0.16 * op.assimilates / tdmwgt)
        m_total = m_pinnae + m_rachis + m_trunk + m_roots +
                  m_maleflo + m_femaflo + m_bunches + m_metabolic
    else
        # all assimilates diverted for maintenance due to unfavorable temperature
        m_total = op.assimilates
    end

    m_pinnae, m_rachis, m_trunk, m_roots, m_maleflo, m_femaflo, m_bunches, m_total
end


"""
    vdm_requirement(pd::Int, lai, Lmax, Vmax)

Vegetative dry matter daily demand/requirement (kg DM/palm/day).

# Arguments
- `pd::Int`: planting density (palms/ha)
- `lai`: leaf area index (m2 leaf/m2 ground)
- `Lmax`: max. leaf area index (m2 leaf/m2 ground)
- `Vmax`: max. annual VDM (kg/palm/yr)

# Returns
- `Float64`: daily vegetative dry matter requirement (kg DM/palm/day)
"""
function vdm_requirement(pd::Int, lai, Lmax, Vmax)
    function vdmmin(lmin)
        ca = 0.935
        inv_ca = 1 / ca
        a = ca / Vmax
        b = 0.1 * (inv_ca - 1) * (pd / 100) ^ inv_ca
        1 / (a + b / lmin^1.5)
    end

    Lmin = 0.1
    Vmin = vdmmin(Lmin)

    if lai < Lmax
        eLmax = Lmax^1.5
        eLmin = Lmin^1.5
        a = (Vmin * eLmax - Vmax * eLmin) / (Vmax * Vmin * (eLmax - eLmin))
        b = (eLmax * eLmin * (1 - Vmin / Vmax)) / (Vmin * (eLmax - eLmin))
        V = 1 / (a + b / lai^1.5)
    else
        n1 = Vmax * Vmin * (Lmax - Lmin)
        n2 = Vmin * Lmax - Vmax * Lmin
        a = n1 / n2
        n1 = (Vmax * Lmin - Vmin * Lmax) ^ 2
        n2 = Lmax * Lmin * (Lmin- Lmax) * Vmax * Vmin * (Vmax - Vmin)
        b = -n1 / n2
        V = a * (1 - 1 / (1 + a * b * lai))
    end
    V / 365  # convert annual VDM to daily VDM
end


"""
    cvf(parts::Parts)

Convert weight in glucose (CH2O) to that in dry matter (DM) (kg DM/kg CH2O).

# Arguments
- `parts::Parts`: tree parts

# Returns
- `Float64`: converted weight in CH2O to kg DM/kg CH2O
"""
function cvf(parts::Parts)
    leaves = parts.pinnae.frac + parts.rachis.frac
    0.7 * leaves + 0.66 * parts.trunk.frac + 0.65 * parts.roots.frac
end


"""
    veg_growth(parts::Parts, assim4growth)

Growth rates for the vegetative tree parts (kg DM/palm/day).

# Arguments
- `parts::Parts`: tree parts
- `assim4growth`: assimilates for growth (kg CH2O/palm/day)

# Returns
- `Tuple{Float64, Float64, Float64, Float64}`:
        growth rates for pinnae, rachis, trunk, and roots (kg DM/palm/day)
"""
function veg_growth(parts::Parts, assim4growth)
    availvdm = assim4growth * cvf(parts)  # convert weight from per CH2O to DM basis
    availvdm .* vegparts(parts, :frac)
end


"""
    veg_death(parts::Parts, treeage::Int)

Death rates for the vegetative tree parts (kg DM/palm/day).

# Arguments
- `parts::Parts`: tree parts
- `treeage::Int`: tree age (days)
- `stress_t`: crop water stress (fraction, 0-1)
- `sand`: amount of sand in the root zone (%)

# Returns
- `Tuple{Float64, Float64, Float64, Float64}`:
        death rates for pinnae, rachis, trunk, and roots (kg DM/palm/day)
"""
function veg_death(parts::Parts, treeage::Int, stress_t, sand)
    # pinnae and rachis death ∝ water stress level and sand content in root zone:
    amp = max(0, min(0.5, 0.6 - 0.01 * sand))
    leafdeath = 0.0013 * (1 + amp * (1 - stress_t))
    lower, upper = 600, 2400
    Δ = (treeage - lower) / (upper - lower)
    dleaves = leafdeath * (Δ * (lower < treeage <= upper) + (treeage > upper))
    dpinnae = dleaves * parts.pinnae.weight
    drachis = dleaves * parts.rachis.weight

    # roots death:
    lower = 1200
    rootsdeath = 0.25
    droots = (treeage > lower) ? (rootsdeath / 365) : 0.0
    droots *= parts.roots.weight

    dpinnae, drachis, 0.0, droots  # no death for trunk
end


"""
    update_veg!(mo::Output, mi::Input)

Increment the various tree part weights (kg DM/palm), and update LAI
(m2 leaf/m2 ground) and respiration rates (kg CH2O/palm/day).

# Arguments
- `mo::Output`: model output crop growth results
- `mi::Input`: model input crop parameters

# Returns
- `nothing
"""
function update_veg!(mo::Output, mi::Input)
    op = mo.op
    @unpack parts, treeage, plantdens, assimilates, lai, laimax, vdmmax = op

    c = cvf(parts)

    # 1. maintenance respiration:
    m = maintenance_respiration(mo, mi)    # kg CH2O/palm/day
    allparts!(parts, :maint, c .* m[1:7])  # store as kg DM/palm/day
    # total maintenance (note: in kg CH2O/palm/day)
    assim4maint = min(assimilates, m[end])

    # 2. growth respiration:
    # max. avail. assimilates for growth (kg CH2O/palm/day)
    maxassim = assimilates - assim4maint
    # VDM growth demand (kg DM/palm/day)
    vdmreq = vdm_requirement(plantdens, lai, laimax, vdmmax)
    # assimilates for vegetative growth (kg CH2O/palm/day)
    assim4growth = min(vdmreq / c, maxassim)
    # growth rates for the individual plant parts (kg DM/palm/day)
    vegparts!(parts, :growth, veg_growth(parts, assim4growth))
    # death rates for the individual plant parts (kg DM/palm/day)
    vegparts!(parts, :death, veg_death(parts, treeage, op.stress_t, op.roots.sand))

    # assimilates for generative growth (kg CH2O/palm/day)
    assim4gen = maxassim - assim4growth

    # 3. update weights (kg DM/palm):
    growth = vegparts(parts, :growth)
    death = vegparts(parts, :death)
    netgrowth = growth .- death
    Δ = vegparts(parts, :weight) .+ netgrowth
    vegparts!(parts, :weight, Δ)

    # 4. pinnae weight has changed, so update the SLA and LAI
    op.sla, op.lai = sla_lai(treeage, plantdens, parts.pinnae.weight, mi.sla_table)

    # 5. update the available assimilates
    op.assim4maint = assim4maint
    op.assim4growth = assim4growth
    op.assim4gen = assim4gen
    op.vdmreq = vdmreq
end


"""
    new_flower_sex(femaprob)

Sex of new flower (0 = male, 1 = female).

# Arguments
- `femaprob`: probability of getting a female flower (0-1)

# Returns
- `Int`: either 0 (male) or 1 (female)
"""
new_flower_sex(femaprob) = (rand() <= femaprob) ? 1 : 0


"""
    gen_growth(parts::Parts, boxcar::BoxCar, newflosex::Int, assim4gen)

Growth rates for generative organs (kg DM/palm/day).

# Arguments
- `parts::Parts`: tree parts
- `boxcar::BoxCar`: flower and bunch development cycles
- `newflosex::Int`: sex of new flower (0 = male or 1 = female)
- 'assim4gen`: assimilates for generative growth (kg CH2O/palm/day)

# Returns
- `Tuple{Float64, Float64, Float64}`:
        growth rates (kg DM/palm/day) for male flowers, female flowers, and bunches
"""
function gen_growth(parts::Parts, boxcar::BoxCar, newflosex::Int, assim4gen)
    # count non-zero weights in boxcar trains:
    # +1 if first gender is male
    n1 = count(>(0), boxcar.maleflo[2:end]) + (1 - newflosex)
    f1 = parts.maleflo.frac * n1 / length(boxcar.maleflo)

    # +1 if first gender is female
    n2 = count(>(0), boxcar.femaflo[2:end]) + newflosex
    f2 = parts.femaflo.frac * n2 / length(boxcar.femaflo)

    n3 = count(>(0), boxcar.bunches)  # may have no bunches
    f3 = parts.bunches.frac * n3 / length(boxcar.bunches)

    ftotal = f1 + f2 + f3  # never 0; there's at least one female or male flower
    f1 /= ftotal
    f2 /= ftotal
    f3 /= ftotal

    cvf2 = 0.7 * f1 + 0.7 * f2 + 0.44 * f3
    g1 = (n1 > 0) ? f1 * cvf2 * assim4gen / n1 : 0.0
    g2 = (n2 > 0) ? f2 * cvf2 * assim4gen / n2 : 0.0
    g3 = (n3 > 0) ? f3 * cvf2 * assim4gen / n3 : 0.0

    g1, g2, g3  # male flowers, female flowers, bunches
end


"""
    update_gen!(mo::Output, mi::Input)

Increment the weights for all generative organs (kg DM/palm). Returns FFB
yield (kg DM/palm/day).

# Arguments
- `mo::Output`: model output crop growth results
- `mi::Input`: model input crop parameters

# Returns
- `nothing
"""
function update_gen!(mo::Output, mi::Input)
    function increment!(box, wgt)
        # increment weights in the individual boxcars
        @inbounds foreach(i -> if (box[i] > 0) box[i] += wgt end, eachindex(box))
    end

    op = mo.op
    bc = mo.bc
    parts = op.parts

    # determine if water stress will abort flower (male and female):
    if rand() > op.stress_t
        bc.maleflo[90] = 0.0  # male and female flowers aborted at node 90
        bc.femaflo[90] = 0.0
    end

    # growth rates (kg DM/palm/day):
    prob = mi.femaleprob  # probability of getting female flowers
    newflosex = new_flower_sex(prob)    # 0 (male) or 1 (female)
    parts.maleflo.growth,
    parts.femaflo.growth,
    parts.bunches.growth = gen_growth(parts, bc, newflosex, op.assim4gen)

    # increment the weights (kg DM/palm) in the generative organ boxcars:
    increment!(bc.maleflo, parts.maleflo.growth)
    increment!(bc.femaflo, parts.femaflo.growth)
    increment!(bc.bunches, parts.bunches.growth)

    # get yield (kg DM/palm/day) which is the last value in the bunches boxcar train:
    op.yield = bc.bunches[end]

    # in-place shift to the right for weights in the boxcar:
    pushfirst!(bc.maleflo, pop!(bc.maleflo))
    pushfirst!(bc.femaflo, pop!(bc.femaflo))
    pushfirst!(bc.bunches, pop!(bc.bunches))

    # update the endpoints (latest flower can be either male or female):
    bc.maleflo[1] = parts.maleflo.growth * (1 - newflosex)
    bc.bunches[1] = bc.femaflo[1]
    bc.femaflo[1] = parts.femaflo.growth * newflosex

    # update the total weights (kg DM/palm):
    parts.maleflo.weight = sum(bc.maleflo)
    parts.femaflo.weight = sum(bc.femaflo)
    parts.bunches.weight = sum(bc.bunches)
end


"""
    update_crop!(mo::Output, mi::Input)

Crop growth and yield for the current date.

# Arguments
- `mo::Output`: model output crop growth results
- `mi::Input`: model input crop parameters

# Notes
- Update the arguments `out::OutputCrop` and `boxcar::BoxCar`.

# Returns
- nothing
"""
function update_crop!(mo::Output, mi::Input)
    op = mo.op
    treeage = mi.treeage + Dates.value(op.date - mi.date)

    # reduce planting density, after thinning, if any:
    thinpd = mi.thinplantdens
    cond = (thinpd > 0) && (treeage >= mi.thinage)
    pd = mi.plantdens * (1 - cond) + thinpd * cond

    op.treeage = treeage
    op.yap = div(op.treeage, 365)
    op.plantdens = pd

    sand = op.roots.sand
    op.treehgt, op.trunkhgt = tree_height(treeage, pd, sand)
    op.leaflength, op.leafwidth = leaf_dimension(treeage)
    update_veg!(mo, mi)
    update_gen!(mo, mi)
    op.tdmwgt, op.vdmwgt, op.flowgt, op.frdwgt = dm_weights(op.parts)
    op.tdmgro, op.vdmgro = dm_growth(op.parts, op.yield)
end


"""
    init_crop!(mo::Output, mi::Input)

Initialize the crop parameters.

# Arguments
- `mo::Output`: model output crop growth results
- `mi::Input`: model input crop parameters

# Returns
- `BoxCar`: male and female flowers and bunches development cycle
"""
function init_crop!(mo::Output, mi::Input)
    # initial tree part weights:
    op = mo.op
    parts = op.parts
    parts.pinnae.weight = pin_wgt = mi.pinnae_wgt
    parts.rachis.weight = mi.rachis_wgt
    parts.trunk.weight = mi.trunk_wgt
    parts.roots.weight = mi.roots_wgt
    parts.maleflo.weight = mi.maleflo_wgt
    parts.femaflo.weight = mi.femaflo_wgt
    parts.bunches.weight = mi.bunches_wgt

    # all DM partitioning are constant with tree age:
    parts.pinnae.frac = 0.23
    parts.rachis.frac = 0.5
    parts.trunk.frac = 0.17
    parts.roots.frac = 1 - parts.pinnae.frac - parts.rachis.frac - parts.trunk.frac
    parts.maleflo.frac = 0.159
    parts.femaflo.frac = 0.159
    parts.bunches.frac = 0.682

    # set all the initial crop parameters:
    op.plantdens = pd = mi.plantdens      # not necessarily constant; may be thinned
    op.treeage = treeage = mi.treeage
    op.yap = div(treeage, 365)
    op.trunkhgt = mi.trunkhgt
    op.treehgt = canopy_height(treeage) + op.trunkhgt
    op.leaflength, op.leafwidth = leaf_dimension(treeage)
    op.sla, op.lai = sla_lai(treeage, pd, pin_wgt, mi.sla_table)
    op.tdmwgt, op.vdmwgt, op.flowgt, op.frdwgt = dm_weights(parts)
    op.laimax = lai_maximum(pd)
    op.vdmmax = vdm_maximum(pd)
end
