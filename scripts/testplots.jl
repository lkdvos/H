
using MPSKit, TensorKit

using Plots
plotlyjs()


##
state = InfiniteMPS([Rep[SU₂](1 => 1)], [Rep[SU₂](1//2 => 20, 3//2 => 20, 5//2 => 10)])

entanglementplot(state,show=true)

##
χ = 10
state_nonsymm = InfiniteMPS([TensorMap(rand, ComplexF64, ℂ^χ*ℂ^3, ℂ^χ)])
entanglementplot(state_nonsymm, show=true)

##
transferplot(state, state, show=true)