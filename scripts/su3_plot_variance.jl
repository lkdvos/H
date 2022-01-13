using DrWatson
@quickactivate :Heisenberg
plotlyjs()

## Get all data
spt = 0

svals = []
Es = []
sigmas = []

datafiles = readdir(datadir("sims", "[3 0 0]", "gs"))
for datafile in datafiles
    prefix, params, suffix = parse_savename(datafile)
    println(params)
    if params["spt"] == spt
        file = wload(datadir("sims", "[3 0 0]", "gs", datafile))
        push!(svals, params["sval"])
        push!(Es, file["E"])
        push!(sigmas, file["sigma"])
    end
end
println(svals)
println(Es)
println(sigmas)

p = scatter(sigmas, Es, 
    title  = "Groundstate energy - variance", 
    xlabel ="σ²", 
    #xlims = (1e-6, 1),
    #ylims = (1e-10, 1e-1),
    yguide ="E₀",
    #xscale = :log10,
    #yscale = :log10,
    #xticks = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1],
    #yticks = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1],
    label = ""
)

model(V, p) = p[1] .+ p[2] .* V
p₀ = [-2.7, 1];

fit = curve_fit(model, sigmas, Es, p₀)
E₀ = first(fit.param)
σ = first(standard_errors(fit))

println("$(E₀) ± $(σ)")

plot!(p, x->model(x, fit.param), label="fit", legend=:topright)


savefig(p, plotsdir("[3 0 0]", "gs", "variance.png"))
savefig(p, plotsdir("[3 0 0]", "gs", "variance.svg"))