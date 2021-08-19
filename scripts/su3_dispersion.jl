using DrWatson, JLD2, MPSKit, TensorKit, Plots, SUNRepresentations
@quickactivate :Heisenberg

H = su3_heis_ham(3,0,2,3);

##

file = jldopen(datadir("sims", "su3300", "matlabstates.jld2"), "r")
keys(file)
gs = file["sval_0.0001291549665014884"]["gs"]
#envs = file["sval_4.6415888336127825e-5"]["envs"]
close(file)

##

sectors = (su3_irrep(0,0), su3_irrep(1,1), su3_irrep(3,0), su3_irrep(0,3), );
momenta = range(-π, π, length=3)
Es = Dict()
Bs = Dict()

for sector in sectors
    (E, B) = excitations(H, QuasiparticleAnsatz(toler=1e-4), momenta, gs, sector=sector)
    Es[sector] = E
    Bs[sector] = B
end

## Save data
