using DrWatson
@quickactivate :Heisenberg
plotlyjs()


spins = [1, 2, 3]
## Prepare input parameters.
allparams = Dict(
    :spin => spins,
    :spt => [0, 1],
    :sval => [-5.0, ],
)
path = datadir("sims", "gs")

## Save groundstates
for params in dict_list(allparams)
    data, s = produce_or_load(path, params, su2_gs_simulations; loadfile=true, tag=false, force=false)
    @unpack gs, envs, delta, E, sigma = data
    println("GS for spin $(params[:spin]), spt $(params[:spt]) at sval $(params[:sval]): E = $E ± $sigma")
    p = entanglementplot(gs, title="SPT=$(params[:spt]), Cut=$(params[:sval])")
    savefig(p, plotsdir("gs", savename("entplot", params, "svg")))
end


## Compute correlation length
allparams = Dict(
    :spin => spins, 
    :spt => [@onlyif(iseven(:spin), 0), @onlyif(isodd(:spin), 1)],
    :svals => [[-4.0, -5.0, -6.0, -7.0, -8.0], ],
    :charges => [[0, 1, 2], ]
)

for params in dict_list(allparams)
    data, s = produce_or_load(datadir("sims", "xi"), params, su2_xi_simulations; tag=false, force=false)
    
    println("Correlation length for spin $(params[:spin]), spt $(params[:spt]):")
    for (key, val) in data["ξs"]
        println("ξ[$key] = $val ± $(data["σs"][key])")
    end
    
    _, _, p = correlation_fit(data["ϵs"], data["δs"]; doplot=true)
    savefig(p, plotsdir("xi", savename("marek", params, "svg")))
end

## Dispersion relation
allparams = Dict(
    :spin => spins, 
    :spt => [@onlyif(iseven(:spin), 0), @onlyif(isodd(:spin), 1)],
    :sval => -5.0,
    :charges => [[0, 1, 2], ],
    :momenta => range(0, π, length=50)
)

for params in dict_list(allparams)
    data, s = produce_or_load(datadir("sims", "gap"), params, su2_gap_simulations; tag=false, force=false)
    @unpack Es, Bs = data
    
    p = plot(
        ;
        title="Dispersion relation", 
        xguide="momentum",
        yguide="Δ"
    )
    for sector in keys(Es)
        plot!(p, params[:momenta], Es[sector], label=string(sector))
    end
    savefig(p, plotsdir("gap", savename(params, "png")))
    savefig(p, plotsdir("gap", savename(params, "svg")))
end
