using DrWatson
@quickactivate :Heisenberg
using MPSKit, MPSKitModels, JLD2, TensorKit

allparams = Dict(
    "spin" => [1,], 
    "spt" => [1], 
    # "charge" => SU2Irrep(1),
    "sval" => collect(-5:-0.5:-9),
)

function gs_simulations(params::Dict)
    H = su2_xxx_ham(spin=params["spin"])
    
    expand_alg = Dynamical(
        alg_groundstate = Vumps(tol_galerkin=1e-5, maxiter=5,verbose=false),
        alg_changebonds = TwoSite(scut=10^params["sval"], tol_eigenval=1e-6),
        maxiter = 25
    )
    converge_alg = GradientGrassmann(maxiter=100,tol=1e-10)
    
    (gs, envs, delta) = find_groundstate(
        initial_state(params["spin"], params["spt"]),
        H,
        expand_alg & converge_alg
    )
    
    E = expectation_value(gs, H)
    σ = variance(gs, H, envs)
    
    println("@sval $(params["sval"]): E=$E±$σ")
    
    # (ϵ,δ,θ) = marek_gap(gs, sector=params["charge"])
    # println("@sval $(-log(params["sval"])): ξ=$(1/ϵ)")
    
    # filename = "../data/sims/su2_xi/spin_$(params["spin"])_spt_$(params["spt"])_sval_$(params["spin"]).jld2"
    
    data = Dict(
        "gs" => gs,
        "envs" => envs,
        "delta" => delta
    )
    return data
    # safesave(datadir("sims", "gs",savename(params,"jld2")), data)
    
end


## Processing the gs_simulations

extrapol = map(dict_list(allparams)) do c
    data, s = produce_or_load(datadir("sims","gs"), c, gs_simulations; suffix="jld2", tag=false)
    gs = data["gs"]
    envs = data["envs"]
    
    (ϵ,δ,θ) = marek_gap(gs, sector=SU2Irrep(1))
    
    return ϵ, δ
end

ϵs = getindex.(extrapol, 1)
δs = getindex.(extrapol, 2)
ξ,σ,p = correlation_fit(ϵs, δs, doplot=true)
display(p)
println("ξ = $ξ ± $σ")