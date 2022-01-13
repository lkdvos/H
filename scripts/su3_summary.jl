using DrWatson
@quickactivate :Heisenberg

file = jldopen(datadir("exp_pro","su3_300_gap.jld2"))


for (E,sigma,sval) in zip(file["E"], file["sigma"], file["svals"])
    println("@sval $(sval): E = $(E) ± $(sigma)")
end

close(file)