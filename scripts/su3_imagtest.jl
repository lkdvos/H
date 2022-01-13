using DrWatson
@quickactivate :Heisenberg
using SUNRepresentations

file = jldopen(datadir("sims", "[3 0 0]", "gs", "spt0_sval-4.0.jld2"))

gs = file["gs"]
H = su3_heis_ham(3,0,2,3)

sector = SUNIrrep((2,1,0))
E, B = excitations(H, QuasiparticleAnsatz(;toler=1e-4), 2π/3, gs; sector=sector)
E
