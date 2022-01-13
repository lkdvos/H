using DrWatson
@quickactivate :Heisenberg

using Heisenberg

allparams = Dict(
    :spt => [0, 1],
    :sval => -4.0,
)

function energy_variance_sim(params)
    @unpack spt, sval = params
    scut = 10.0^sval
    
    H = su3_heis_ham(3, 0, 2, 3)
    expand_alg = Dynamical(
        alg_groundstate = Vumps(tol_galerkin=1e-4, verbose=false),
        alg_changebonds = VumpsSvdCut(tol_galerkin=1e-4, tol_eigenval=1e-6, tol_gauge=1e-6, trscheme=truncbelow(scut)),
        maxiter = 50
    )
    converge_alg = GradientGrassmann(tol=1e-10, maxiter=200)
    
    init = find_init(scut)
    
    (gs, envs, delta) = find_groundstate(init, H, expand_alg & converge_alg)
    E = first(real.(expectation_value(gs, envs)))
    sigma = variance(gs, H, envs)
    
    return @strdict gs envs delta E sigma
end

function find_init(scut)
    file = jldopen(datadir("sims", "su3300", "matlabstates.jld2"))
    svals_str = keys(file)
    
    svals_fl = broadcast(x->parse(Float64, x[6:end]), svals_str)
    
    sval_fl = maximum(filter(x-> x < scut, svals_fl))
    sval_str = "sval_" * string(sval_fl)
    init = file[sval_str]["gs"]
    close(file)
    return init
end


path = datadir("sims", "[3 0 0]", "gs")
ps = []
p2s = []

for params in dict_list(allparams)
    data, s = produce_or_load(path, params, energy_variance_sim; loadfile = true, tag=false, force=false)
    @unpack gs, envs, delta, E, sigma = data
    println("Found thingy @ spt $(params[:spt])")
    push!(ps, entanglementplot(gs; sector_formatter=label, title="SPT $(params[:spt])"))
    push!(p2s, transferplot(gs, gs; sectors=[SUNIrrep((0,0,0)), SUNIrrep((3,0,0)), SUNIrrep((3,3,0)), SUNIrrep((2,1,0)), SUNIrrep((4,2,0))], sector_formatter=label, title="SPT $(params[:spt])"))
end

p = plot(ps...; link = :y)
savefig(p, plotsdir("[3 0 0]", "gs", savename("sptplot", params, "png")))
savefig(p, plotsdir("[3 0 0]", "gs", savename("sptplot", params, "svg")))

p2 = plot(p2s...)
savefig(p2, plotsdir("[3 0 0]", "gs", savename("sptplot2", params, "png")))
savefig(p2, plotsdir("[3 0 0]", "gs", savename("sptplot2", params, "svg")))