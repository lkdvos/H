using DrWatson
@quickactivate :Heisenberg

file = jldopen(datadir("exp_pro","su3_300_gap.jld2"))
keys(file)
# file["sval_7.742636826811278e-7"]["E"]
##
close(file)