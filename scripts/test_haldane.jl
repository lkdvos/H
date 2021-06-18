
include("../src/header.jl")
plotlyjs()

##

spin = 1
H = su2_xxx_ham(spin=spin)
expand_alg = Dynamical(
    alg_groundstate = Vumps(tol_galerkin = 1e-6, maxiter=10, verbose=true),
    alg_changebonds = TwoSite(scut=1e-8,tol_eigenval=1e-6),
    maxiter = 25
)

(gs, envs, _) = find_groundstate(
    initial_state(spin, 1), H, expand_alg
)

p = entanglementplot(gs, show=true)

(gs2, envs2, _) = find_groundstate(gs, H, GradientGrassmann(maxiter=100,tol=1e-10))


println("ξ=$(correlation_length(gs2; sector=SU2Irrep(1)))")

spectrum=transfer_spectrum(gs2,sector=SU2Irrep(1))
println("ξ=$(-1/log(abs(first(spectrum))))")

##
(E, B) = excitations(H, QuasiparticleAnsatz(), 1.0π, gs2, envs2, sector=SU2Irrep(1));
println("$(real(first(E))) ± $(variance(B[1], H))")
