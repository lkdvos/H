using Revise: push!
include("../src/header.jl")

t1 = TensorMap(rand, ComplexF64, ℂ^3*ℂ^2, ℂ^4)
t2 = TensorMap(zeros, ComplexF64, (ℂ^5) * ℂ^2, ℂ^2)

a = codomain(t1)
b, = codomain(t1)
c = (codomain(t1)...,)

V1 = codomain(t1)
V2 = codomain(t2)

@assert isdual.(V1) == isdual.(V2) 
prod(supremum.(V1, V2))

blocksectors(t1)

embed!(t2, t1)

##
t3 = TensorMap(rand, ComplexF64, SU2Space(1=>1), SU2Space(1=>1))
t4 = TensorMap(rand, ComplexF64, SU2Space(2=>1), SU2Space(2=>1))


a = space(t3,1)
GradedSpace(SU2Irrep(2)=>2)