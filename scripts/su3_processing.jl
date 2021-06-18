using DrWatson, JLD2, MPSKit, TensorKit, Plots
@quickactivate :Heisenberg

H = su3_heis_ham(3,0,2,3)
##
file = jldopen(datadir("sims","su3300","excitations.jld2"), "r+")
##
Es = []
σs = []
svals = []

for group in keys(file)
    sval = parse(Float64, group[6:end])
    println(sval)
    if sval < 1e+4
        push!(Es, file[group]["E"])
        push!(σs, variance(first(file[group]["B"]), H))
        push!(svals, sval)
    end
end

close(file)

data = Dict(
    "E"=> Es,
    "sigma"=> σs,
    "svals" => svals
)

wsave(datadir("exp_pro", "su3_300_gap.jld2"), data)


## Correlation length
sectors = [su3_irrep(0,0), su3_irrep(1,1), su3_irrep(3,0), su3_irrep(0,3), su3_irrep(2,2)]
ϵs = Dict()
δs = Dict()
θs = Dict()

for sector in sectors
    ϵs[sector] = []
    δs[sector] = []
    θs[sector] = []
end

jldopen(datadir("sims", "su3300", "matlabstates.jld2"), "r") do file
    svals = keys(file)
    for sval in sort(svals, by=x->parse(Float64, x[6:end]), rev=true)
        println(sval)
        if parse(Float64, sval[6:end]) < 1e-4
            gs = file[sval]["gs"]
            for sector in sectors
                (ϵ,δ,θ) = marek_gap(gs, sector=sector, num_vals=10)
                push!(ϵs[sector], ϵ)
                push!(δs[sector], δ)
                push!(θs[sector], θ)
            end
        end
    end
    
end

##
data = wload(datadir("exp_pro", "epsdelta.jld2"))
ϵs = data["epsilon"]
δs = data["delta"]

##
ξ, σ, p = correlation_fit(ϵs, δs, doplot=true)


plotlyjs()
plot(p, legend=:bottomright, title="Correlation length extrapolation", xlims=(0,0.08), ylims=(0,0.1), yguide="ϵ", xguide="δ")
savefig(plotsdir("epsdelta.eps"))

data = Dict("xi"=>ξ, "sigma" => σ, "theta" => θs, "epsilon" => ϵs, "delta"=>δs)
wsave(datadir("exp_pro", "epsdelta.jld2"), data)
##

jldopen(datadir("sims", "su3300", "matlabstates.jld2"), "r") do file
    svals = keys(file)
    for sval in sort(svals, by=x->parse(Float64, x[6:end]), rev=true)
        println(sval)
        if parse(Float64, sval[6:end]) < 1e-4
            gs = file[sval]["gs"]
            E = file[sval]["E"]
            println(sval, "$(real(file[sval]["E"])) ± $(file[sval]["sigma"]) @ $(file[sval]["chi"])")
        end
    end
    
end
