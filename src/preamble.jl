using CairoMakie
using DataFrames

function available_sites(include_main::Bool=true, include_kalumpong::Bool=false)
    sites = []

    if include_main
        merlimau = [122, 135, 148, 164, 181, 199, 220, 243, 268, 296]
        rengam = [50, 55, 61, 67, 74, 82, 90, 100, 111, 123, 136, 166, 184]
        sg_buloh = [ 64,  70,  77,  85,  94, 103, 114, 125, 136,
                    138, 152, 167, 184, 203, 223, 246, 271, 299]
        sabrang = [136, 170]
        seri_intan = [180]
        kerayong = [136, 160, 185]

        append!(sites, [
            ("merlimau", merlimau),
            ("rengam", rengam[6:end]),
            ("sg-buloh", sg_buloh[3:end]),
            ("sabrang", sabrang),
            ("seri-intan", seri_intan),
            ("kerayong", kerayong)
        ])
    end

    if include_kalumpong
        kalumpong_1997 = [124]
        kalumpong_1998 = [126]
        kalumpong_2000A = [138]
        kalumpong_2001A1 = [154]
        kalumpong_2001A2 = [150]
        kalumpong_2002A1 = [146]
        kalumpong_2002A2 = [158]
        kalumpong_2002B1 = [150]
        kalumpong_2003A1 = [150]
        kalumpong_2006A2 = [159]

        append!(sites, [
            ("kalumpong-1997", kalumpong_1997),
            ("kalumpong-1998", kalumpong_1998),
            ("kalumpong-2000A", kalumpong_2000A),
            ("kalumpong-2001A1", kalumpong_2001A1),
            ("kalumpong-2001A2", kalumpong_2001A2),
            ("kalumpong-2002A1", kalumpong_2002A1),
            ("kalumpong-2002A2", kalumpong_2002A2),
            ("kalumpong-2002B1", kalumpong_2002B1),
            ("kalumpong-2003A1", kalumpong_2003A1),
            ("kalumpong-2006A2", kalumpong_2006A2),
        ])
    end

    sites
end


function site(loc::AbstractString; include_main::Bool=true, include_kalumpong::Bool=false)
    locs = available_sites(include_main, include_kalumpong)
    pds = locs[loc .== first.(locs)][1][2]
    dir = "data/$loc/"
    (; pds, dir)
end


function multirun(jsonfname::AbstractString, site::AbstractString, pds;
                  gof::AbstractVector{<:Function}=[GoF.NMAE, GoF.NMBE, GoF.KGE],
                  plot_valid::Bool=false, saveplot_valid::Bool=false,
                  plot_combo::Bool=false, saveplot_combo::Bool=false,
                  plot_index::Bool=false, saveplot_index::Bool=false)

    dir = dirname(jsonfname)

    scenarios = map(v -> Dict(
        :plantdens => v,
        :outputfname => dir * "/$(site)-$(v).csv",
        ),
        pds
    )

    @tictoc begin
        println("Running modeling scenarios..."); flush(stdout)
        res = start(jsonfname, scenarios)
    end

    @tictoc begin
        println("Validating model results..."); flush(stdout)
        valid = map(eachindex(pds)) do i
            validate("obs/$(site)-obs-$(pds[i]).csv",
                     res[i].annl, res[i].mi;
                     gof=gof,
                     plot=plot_valid,
                    #  title="$(site), PD = $(pds[i])",
                     title="",
                     saveplot=saveplot_valid,
                     plotfname="$(dir)/$(site)-$(pds[i]).png")
        end
    end

    round_df!([valid[i].gofdf for i in eachindex(valid)], 3)

    println("\n")
    for i ∈ eachindex(valid)
        println("**** PD = $(pds[i]) ****")
        show(valid[i].gofdf, allrows=true, allcols = true)
        println("\n")
    end

    if plot_combo
        nm = names(valid[1].pairdf)
        par_sim = nm[occursin.(r"_sim$", nm)]
        par_obs = replace.(par_sim, r"_sim$" => "")

        obs = DataFrame()
        sim = DataFrame()
        for i ∈ eachindex(par_obs)
            po, ps = par_obs[i], par_sim[i]
            o = map(v -> getproperty(v.pairdf, po), valid)
            s = map(v -> getproperty(v.pairdf, ps), valid)
            setproperty!(obs, po, vcat(o...))
            setproperty!(sim, po, vcat(s...))
        end

        plot_pairs(obs, sim;
                   gof=gof,
                   title="",
                #    title="$(site) (all PDs)",
                   saveplot=saveplot_combo,
                   plotfname="$(dir)/$(site)-all.png")
    end

    if plot_index
        println("Plotting GoF indexes..."); flush(stdout)
        nparams = size(valid[1].gofdf)[1]   # no. of parameters
        pds = string.(map(i -> scenarios[i][:plantdens], eachindex(scenarios)))

        foreach(eachindex(valid)) do i
            valid[i].gofdf.pd = repeat([pds[i]], nparams)
        end

        gofdf = vcat(map(i -> valid[i].gofdf, eachindex(valid))...)
        gdf = groupby(gofdf, :param)

        indexes = [n for n ∈ names(gofdf) if ((n != "param") && (n != "pd"))]

        p = valid[1].gofdf.param
        x = Vector{Vector{Float64}}(undef, nparams)
        y = Vector{Vector{Float64}}(undef, nparams)
        z = Vector{Vector{Float64}}(undef, nparams)

        for (i, g) ∈ enumerate(gdf)
            x[i] = g[:, indexes[1]]
            y[i] = g[:, indexes[2]]
            z[i] = g[:, indexes[3]]
        end

        plot_gof(p, pds, x, y, z;
                 zmin=0, zmax=1, xopt=0, yopt=0,
                 xlabel=indexes[1], ylabel=indexes[2], zlabel=indexes[3],
                 flabel=1.7,
                 fdist=0.4,
                 link_axes=true,
                 mesh_on=false,
                 saveplot=saveplot_index, plotfname="$(dir)/$(site)-gof.png",
                 ncolors=10)
    end

    (; scenarios, res, valid)
end
