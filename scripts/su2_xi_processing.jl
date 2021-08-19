using DrWatson
@quickactivate :Heisenberg

filename = datadir("sims", "xi", "spin=1_spt=1.jld2")

file = jldopen(filename)
println(file["ξs"])
println(file["σs"]
close(file)