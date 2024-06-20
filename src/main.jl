using Sawit


jsonfname = "data/merlimau/input.json"

@tictoc begin
    res = start(jsonfname)
    plot_annual(res.annl, res.mi; saveplot=false)
    plot_daily(res.daly, res.mi; saveplot=false)
end

# @tictoc begin
#     res = runsimhour!(jsonfname;
#                       saveoutput=true,
#                       outfname="hour.csv",
#                       plot=true,
#                       saveplot=false)
# end

;
