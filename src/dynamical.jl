export Dynamical
@with_kw struct Dynamical{F,G} <: MPSKit.Algorithm where {F <: MPSKit.Algorithm, G <: MPSKit.Algorithm}
    alg_groundstate::F
    alg_changebonds::G
    maxiter::Int = 20
end

function MPSKit.find_groundstate(state, H, alg::Dynamical, envs=environments(state, H))
    iter = 0
    @info "dynamical @iteration $iter: $(virtualspace(state, 0))"
    
    while true
        iter += 1
        (state, envs, _) = find_groundstate(state, H, alg.alg_groundstate, envs)
        vspace1 = virtualspace(state, 0)
        (state, envs) = changebonds(state, H, alg.alg_changebonds, envs)
        vspace2 = virtualspace(state, 0)
        @info "dynamical @iteration $iter: $vspace2)"
        
        if vspace1 == vspace2 || iter >= alg.maxiter
            return find_groundstate(state, H, alg.alg_groundstate, envs)
        end
    end
end

export Haldane
@with_kw struct Haldane <: MPSKit.Algorithm
    scut::Float64
end

function MPSKit.changebonds(state, H, alg::Haldane, envs=environments(state,H))
    
    newspaces = []
    for i in 1:length(state)
        (_,s,_) = tsvd(state.CR[i], trunc=truncbelow(alg.scut))
        
        # Expansion of the different blocks.
        newspace = SU2Space()
        for (sector, block) in blocks(s)
            svals = diag(block)
            bond = length(svals)
            if minimum(svals) > alg.scut
                bond += div(bond, 10, RoundUp)
            end
            newspace = newspace ⊕ SU2Space(sector.j => bond)
            
            if !hassector(space(s,1), SU2Irrep(sector.j+1))
                newspace = newspace ⊕ SU2Space(sector.j+1 => 1)
            end
        end
        
        push!(newspaces, newspace)
        
    end
    
    if newspaces != space.(state.CR, (1,))
        
        newALs = similar(state.AL)
        
        for i in 1:length(state)
            newALs[i] = 1e-3 * TensorMap(
                rand, eltype(state.AL[i]), 
                newspaces[i]* space(state, i), newspaces[i]
            )
            embed!(newALs[i], state.AL[i])
        end
        
        state = InfiniteMPS(newALs)
        envs = environments(state,H)
    end
    
    return state, envs
    
end

function fillblock!(src, dst)
    rows = min(size(src, 1), size(dst, 1))
    cols = min(size(src, 2), size(dst, 2))
    dst[1:rows, 1:cols] .= src[1:rows, 1:cols]
end

function embed!(tdst::TensorMap, tsrc::TensorMap)
    numin(tsrc) == numin(tdst) || throw(SpaceMismatch("Cannot embed tensors with different number of out legs."))
    numout(tsrc) == numout(tdst) || throw(SpaceMismatch("Cannot embed tensors with different number of out legs."))
    
    V1 = codomain(tsrc)
    V2 = codomain(tdst)
    isdual.(V1) == isdual.(V2) || throw(SpaceMismatch("Cannot embed tensors whose codomain has non-matching duality."))
    
    W1 = domain(tsrc)
    W2 = domain(tdst)
    isdual.(W1) == isdual.(W2) ||throw(SpaceMismatch("Cannot embed tensors whose domain has non-matching duality"))
    
    for (c, tblock) in blocks(tdst)
        if c ∈ blocksectors(tsrc)
            fillblock!(block(tsrc, c), tblock)
        end
    end
end

export TwoSite
@with_kw struct TwoSite <: MPSKit.Algorithm
    scut::Float64
    tol_scut::Float64
    tol_eigenval::Float64
end

function MPSKit.changebonds(state::InfiniteMPS, H::Hamiltonian, alg::TwoSite, envs=environments(state,H))
    length(state) == 1 || throw(ArgumentError("Not implemented for multi-site mps"))
    state = changebonds(state, SvdCut(trscheme=truncbelow(alg.scut / alg.tol_scut)))
    @tensor AC2[-1 -2; -3 -4] := state.AC[1][-1,-2,1] * state.AR[1][1,-3,-4]
    
    (vals, vecs,_) = eigsolve(@closure(x->ac2_prime(x, 1, state, envs)), AC2, 1, :SR, tol=alg.tol_eigenval; ishermitian=false)
    AC2 = vecs[1]
    
    (_,C,_) = tsvd!(AC2)
    normalize!(C)
    
    addspaces = []
    for (c, b) in blocks(C)
        if maximum(diag(b)) > alg.scut && !hassector(virtualspace(state,1), c)
            push!(addspaces, GradedSpace(c => 1))
        end
    end
    
    (_,s,_) = tsvd(state.CR[1])
    for (c, b) in blocks(s)
        if minimum(diag(b)) > alg.scut * alg.tol_scut
            bond = min(div(length(diag(b)), 5, RoundUp), 20)
            push!(addspaces, GradedSpace(c => bond))
        end
    end
    
    
    
    if length(addspaces) > 0
        addspaces = reduce(⊕, addspaces)

        AL = 1e-3 * TensorMap(
            rand, eltype(state.AL[1]), 
            (space(state.AL[1], 1) ⊕ addspaces) * space(state.AL[1], 2), 
            space(state.AL[1], 1) ⊕ addspaces
        )
        
        embed!(AL, state.AL[1])
        state = InfiniteMPS([AL])
        envs = environments(state, H)
    end
    
    return state, envs
end