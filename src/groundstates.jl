function su2_gs_simulations(params)
    @unpack spin, spt, sval = params
    scut = 10.0^sval
    
    # Prepare input parameters for finding the groundstate.
    H = su2_xxx_ham(spin=spin)
    expand_alg = Dynamical(
        alg_groundstate = Vumps(tol_galerkin=1e-4, verbose=false),
        alg_changebonds = TwoSite(scut=scut, tol_eigenval=1e-6, tol_scut=1.5),
        maxiter = 50
    )
    converge_alg = GradientGrassmann(tol=1e-10)
    initial = su2_initial_state(spin, spt)
    
    
    (gs, envs, delta) = find_groundstate(initial, H, expand_alg & converge_alg)
    E = first(real.(expectation_value(gs, envs)))
    sigma = variance(gs, H, envs)
    
    return @strdict gs envs delta E sigma
end

function su2_xi_simulations(params)
    @unpack spin, spt, svals, charges = params
    
    ϵs = Dict()
    δs = Dict()
    θs = Dict()
    
    for charge in charges
        sector = SU2Irrep(charge)
        ϵs[sector] = []
        δs[sector] = []
        θs[sector] = []
        for sval in svals
            gs_params = Dict(:spt => spt, :spin => spin, :sval => sval)
            data, _ = produce_or_load(datadir("sims", "gs"), gs_params, su2_gs_simulations; tag=false)
            (ϵ, δ, θ) = marek_gap(data["gs"]; sector=sector)
            push!(ϵs[sector], ϵ)
            push!(δs[sector], δ)
            push!(θs[sector], θ)
        end
    end
    
    ξs, σs = correlation_fit(ϵs, δs)
    
    return @strdict ϵs δs θs ξs σs
end

function su2_gap_simulations(params)
    @unpack spin, spt, sval, charges, momenta = params
    
    gs_params = Dict(:spt => spt, :spin => spin, :sval => sval)
    data, _ = produce_or_load(datadir("sims", "gs"), gs_params, su2_gs_simulations; tag=false)
    
    H = su2_xxx_ham(spin=spin)
    
    Es = Dict()
    Bs = Dict()
    for charge in charges
        sector = SU2Irrep(charge)
        E, B = excitations(H, QuasiparticleAnsatz(), momenta, data["gs"], data["envs"], sector=sector)
        maximum(imag.(E)) < 1e-12 || @warn "Imaginairy components for the energy"
        E = real.(E)
        Es[sector] = E
        Bs[sector] = B
    end

    return @strdict Es Bs
end