include("../src/header.jl")
plotlyjs()

## Prepare input parameters
allparams = Dict(
    "spin" => [1,], 
    "spt" => [0, 1], 
    "chi" => [100, 200, 300],
)

## Define simulations
function spt_simulations(params::Dict)
    H = su2_xxx_ham(spin=params["spin"])
    
    expand_alg = Dynamical(
        alg_groundstate = Vumps(tol_galerkin = 1e-4,verbose=false), 
        alg_changebonds = SvdCut(trscheme=truncdim(params["chi"])) &
            OptimalExpand(trscheme = truncdim(div(params["chi"], 5)))
    )
    
    _finalize(iter, state, H, envs) = changebonds(state, H, SvdCut(trscheme=truncdim(params["chi"])), envs)
    converge_alg = Vumps(maxiter=200,verbose=false, finalize=_finalize)
    
    (gs, envs, _) = find_groundstate(
        initial_state(params["spin"], params["spt"]),
        H,
        expand_alg & converge_alg
    )
    
    E = expectation_value(gs, H)
    σ = variance(gs, H, envs)
    return gs, envs, E, σ
end

filename(params::Dict) = "gs_spt_$(params["spt"])_chi_$(params["chi"]).jld2" 

## Run simulations
Es = []
σs = []
for paramset in expand_params(allparams)
    gs, envs, E, σ = spt_simulations(paramset)
    println(filename(paramset))
    println("E = $E ± $σ")
    push!(Es, real(E))
    push!(σs, σ)
    p = entanglementplot(gs, show=true)
end
