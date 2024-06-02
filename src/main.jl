using Sawit

# include("preamble.jl")


# loc = "sg-buloh"
# pds, dir = site(loc; include_main=false, include_kalumpong=true)
# jsonfname = dir * "input.json"
# println("Site = $loc")
# println("Planting densities = $pds")
# println("Input JSON file = $(pwd())/$jsonfname")

# gof = [GoF.NMAE, GoF.NMBE, GoF.KGE]
# res = multirun(jsonfname, loc, pds;
#                gof=gof,
#                plot_valid=true, saveplot_valid=true,
#                plot_combo=true, saveplot_combo=true,
#                plot_index=true, saveplot_index=true)

# @tictoc begin
#     res = runsimhour!(jsonfname;
#                       saveoutput=true,
#                       outfname="hour.csv",
#                       plot=true,
#                       saveplot=false)
# end


jsonfname = "data/merlimau/input.json"
@tictoc begin
    res = start(jsonfname)
    plot_annual(res.annl, res.mi; saveplot=false)
    plot_daily(res.daly, res.mi; saveplot=false)
end

# pairdf = DataFrame(pd=Int[], yap=Float64[], yield=Float64[], yield_sim=Float64[])
# cnt = size(res.valid)[1]
# for i ∈ 1:cnt
#     pd = res.res[i].mi.plantdens
#     df = select(res.valid[i].pairdf, [:yap, :yield, :yield_sim])
#     dropmissing!(df)
#     df.pd .= pd
#     append!(pairdf, df)
# end

# CSV.write("pd.csv", pairdf)
# println(pairdf)
# println("done.")
;
