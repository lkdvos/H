module Heisenberg
using Revise: push!
__precompile__(false)

using Reexport
@reexport using Revise
@reexport using TensorKit
@reexport using MPSKit
@reexport using Plots
@reexport using KrylovKit
@reexport using JLD2
using SUNRepresentations
using Parameters
using LsqFit
using MPSKitModels
using FastClosures
using DrWatson
using LinearAlgebra: diag

include("dynamical.jl")

export su2_initial_state
include("initialisation.jl")

include("groundstates.jl")
export su2_gs_simulations, su2_xi_simulations, su2_gap_simulations


export correlation_fit
include("fitting.jl")
include("utility.jl")

export su3_heis_ham
function su3_heis_ham(m::Int, n::Int, J::Number = 1.0, E₀::Number = 0.0)
    # Hamiltonian for the SU(3) Heisenberg model
    #   H = J ∑ (Sᵢ•Sⱼ) + E₀
    #   where the S are the SU(3) generators, normalised such that tr(SᵃSᵇ) = 15//2 δᵃᵇ.

    physical_space = RepresentationSpace(su3_irrep(m, n) => 1)
    adjoint_space  = RepresentationSpace(su3_irrep(1, 1) => 1)
    trivial_space  = RepresentationSpace(su3_irrep(0, 0) => 1)

    Sᵢ = TensorMap(ones, ComplexF64, trivial_space * physical_space, adjoint_space * physical_space) * sqrt(6)
    Sⱼ = TensorMap(ones, ComplexF64, adjoint_space * physical_space, trivial_space * physical_space) * sqrt(6)

    mpo_data = Array{Any, 3}(missing, 1, 3, 3)
    mpo_data[1,1,1] = 1.0
    mpo_data[1,3,3] = 1.0

    mpo_data[1,1,2] = J * Sᵢ
    mpo_data[1,2,3] = Sⱼ
    mpo_data[1,1,3] = E₀

    return MPOHamiltonian(mpo_data)

end

export su3_irrep
function su3_irrep(m, n)
    SUNRepresentations.SUNIrrep((m + n, n, 0))
end

end