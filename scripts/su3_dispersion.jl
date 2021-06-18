using DrWatson, JLD2, MPSKit, TensorKit, Plots
@quickactivate :Heisenberg

H = su3_heis_ham(3,0,2,3);

##

file = jldopen(datadir("sims", "su3300", "matlabstates.jld2"), "r")
keys(file)
gs = file["sval_4.6415888336127825e-5"]["gs"]
envs = file["sval_4.6415888336127825e-5"]["envs"]
close(file)

##

sectors = (su3_irrep(0,0), su3_irrep(1,1), su3_irrep(3,0), su3_irrep(0,3), su3_irrep(2,2));
momenta = range(0, π, length=31)
Es = Dict()
Bs = Dict()

for sector in sectors
    (E, B) = excitations(H, QuasiparticleAnsatz(toler=1e-4), momenta, gs, envs, sector=sector)
    Es[sector] = E
    Bs[sector] = B
end
