function correlation_fit(ϵs, δs; doplot=false)
    # Assume a linear fit for the extrapolation.
    model(δ, p) = p[1] .+ p[2] .* δ 
    p₀ = [0.0, 1.0]
    
    fit = curve_fit(model, δs, ϵs, p₀)
    
    # Extract the data for the interception at δ = 0.
    ϵ₀ = first(fit.param)
    σ_ϵ = first(standard_errors(fit))
    
    # Compute the correlation length and propagate the error.
    ξ = 1/ϵ₀
    σ = (σ_ϵ/ϵ₀) * ξ
    
    if doplot
        p = scatter(δs, ϵs, label="ϵ(δ)")
        xlims!(0,Inf)
        plot!(p, x->model(x, fit.param), label="fit", legend=:topleft)
        return ξ, σ, p
    end
    
    return ξ, σ
end

function correlation_fit(ϵs::Dict, δs::Dict; doplot=false, label=string)
    plotlyjs()
    model(δ, p) = p[1] .+ p[2] .* δ
    p₀ = [0.0, 1.0]
    
    ξ = Dict()
    σ = Dict()
    
    if doplot
        p = plot(; xguide="δ", yguide="ϵ", title="Correlation length extrapolation")
        xlims!(0,Inf)
    end
    
    for sector in keys(ϵs)
        
        fit = curve_fit(model, δs[sector], ϵs[sector], p₀)
        
        # Extract the data for the interception at δ = 0.
        ϵ₀ = first(fit.param)
        σ_ϵ = first(standard_errors(fit))
        
        ξ[sector] = 1/ϵ₀
        σ[sector] = (σ_ϵ/ϵ₀) * ξ[sector]
        
        if doplot
            plot!(p, x->model(x, fit.param), color=:black, label="")
            scatter!(δs[sector], ϵs[sector], label=label(sector))
        end
    end
    if doplot
        return ξ, σ, p
    else
        return ξ, σ
    end
end