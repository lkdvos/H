using DrWatson
@quickactivate :Heisenberg
plotlyjs()
using TensorKit, SUNRepresentations, Plots, MPSKit, TensorOperations

params = Dict(:spt => [0, 1], :sval => -4.0)
path = datadir("sims", "[3 0 0]", "gs")
sectors = [SUNIrrep((4, 2, 0)) SUNIrrep((3, 3, 0)) SUNIrrep((3, 0, 0)) SUNIrrep((2, 1, 0)) SUNIrrep((0, 0, 0))]

p = []


for param in dict_list(params)
data, s = produce_or_load(path, param, su3_gs_simulations; loadfile = true, tag = false, force = false)

push!(p, 
    plot(; 
        xlims = (-1, 1), 
        ylims = (-1, 1), 
        xwiden=true, 
        ywiden=true,
        framestyle=:origin 
    )
)
@unpack gs = data
spectra = map(sectors[end:end]) do sector
    
    init = TensorMap(rand, eltype(eltype(gs)), MPSKit.left_virtualspace(gs,0), ℂ[typeof(sectors[1])](sectors[1]=>1)' * MPSKit.left_virtualspace(gs,0)) 
    ding = transfer_left(init, gs.AL, gs.AL)
    println(norm(ding))
    
    # @tensor transfermatrix[-1 -2;-3 -4] = gs.AL[-1 1 -3] * conj(gs.AL[-2 1 -4])
    # println(norm(transfermatrix))
    
    spectrum = transfer_spectrum(gs; sector = sector, num_vals = 20)
    scatter!(spectrum; group=map(i->label(sector),1:length(spectrum)))
    return spectrum
end
end
plot(p...)



#     init = TensorMap(rand, eltype(eltype(gs)), MPSKit.left_virtualspace(gs,0), ℂ[typeof(sectors[1])](sectors[1]=>1)' * MPSKit.left_virtualspace(gs,0))
#     transfer_left(init, gs.AL, gs.AL)
#     transfer_spectrum(gs,num_vals=1)
#     transferplot(gs, gs)
#     push!(p, transferplot(gs, gs; sectors=sectors, sector_formatter = label, num_vals=5, thetaorigin=-pi, xticks=range(-pi, pi, length=7)))
# end
