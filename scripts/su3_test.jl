using DrWatson
@quickactivate :Heisenberg
plotlyjs()
using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(8)

## Prepare input parameters.
allparams = Dict(
    :spt => [0, 1, 2],
    :sval => [-3.0, ],
)
path = datadir("sims", "[3 0 0]", "gs")

## Save groundstates
for params in dict_list(allparams)
    data, s = produce_or_load(path, params, su3_gs_simulations; loadfile=true, tag=false, force=false)
    @unpack gs, envs, delta, E, sigma = data
    println("GS for spt $(params[:spt]) at sval $(params[:sval]): E = $E ± $sigma")
    p = entanglementplot(gs, title="SPT=$(params[:spt]), Cut=$(params[:sval])", sector_formatter=label)
    savefig(p, plotsdir("[3 0 0]", "gs", savename("entplot", params, "svg")))
end


## Compute correlation length
allparams = Dict(
    :spt => 0,
    :svals => [[-3.0, -3.5, -4.0, -4.5, -5.0, -5.5, -6.0], ],
    :charges => [[(0,0,0), (2,1,0), (4,2,0), (3,0,0), (3,3,0)], ]
)

for params in dict_list(allparams)
    data, s = produce_or_load(datadir("sims", "[3 0 0]", "xi"), params, su3_xi_simulations; tag=false, force=false)
    
    println("Correlation length for spt $(params[:spt]):")
    for (key, val) in data["ξs"]
        println("ξ[$key] = $val ± $(data["σs"][key])")
    end
    
    _, _, p = correlation_fit(data["ϵs"], data["δs"]; doplot=true, label=label)
    savefig(p, plotsdir("[3 0 0]", "xi", savename("marek", params, "svg")))
end

## Dispersion relation
allparams = Dict( 
    :spt => 0,
    :sval => -4,
    :charges => [[(2,1,0), (4,2,0), (3,0,0), (3,3,0)], ],
    :momenta => range(0, π, length=121)
)

for params in dict_list(allparams)
    data, s = produce_or_load(datadir("sims", "[3 0 0]", "gap"), params, su3_gap_simulations; tag=false, force=false)
    @unpack Es, Bs = data
    
    p = plot(
        ;
        title="Dispersion relation", 
        xguide="momentum",
        yguide="Δ"
    )
    for sector in keys(Es)
        plot!(p, params[:momenta], Es[sector], label=label(sector))
    end
    savefig(p, plotsdir("[3 0 0]", "gap", savename(params, "png")))
    savefig(p, plotsdir("[3 0 0]", "gap", savename(params, "svg")))
end
