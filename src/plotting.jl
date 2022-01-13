using RecipesBase
using LinearAlgebra:diag
using TensorKit, MPSKit, SUNRepresentations

export label
label(sector::Sector) = string(sector)
label(sector::SU2Irrep) = string(sector.j)
label(sector::SUNIrrep) = "[$(sector.I[1]) $(sector.I[2]) $(sector.I[3])]"
label(sector::Trivial) = "Trivial"
#= 
#= const matlab_colors = [RGB(0.000, 0.447, 0.741), RGB(0.850, 0.325, 0.098),
    RGB(0.929, 0.694, 0.125), RGB(0.494, 0.184, 0.556), 
    RGB(0.466, 0.674, 0.188), RGB(0.301, 0.745, 0.933), 
    RGB(0.635, 0.078, 0.184), RGB(0.250, 0.250, 0.250)]
 =#


@userplot struct EntanglementPlot{T<:Tuple{InfiniteMPS}}
    args::T
end

@recipe function f(h::EntanglementPlot; site=0, expand_symmetry=false, sortby=maximum)

    mps = h.args[1]
    site <= length(mps) || 
        throw(ArgumentError("Not a valid site for the given mps."))
    
    
    (_, s, _) = tsvd(mps.CR[site])
    sectors = blocksectors(s)
    spectrum = []
    for sector in sectors
        partial_spectrum = diag(block(s, sector))
        
        # Duplicate entries according to the quantum dimension.
        if expand_symmetry
            partial_spectrum = repeat(partial_spectrum, dim(sector))
            sort!(partial_spectrum, rev = true)
        end
        
        push!(spectrum, diag(block(s, sector)))
    end
    
    if length(spectrum) > 1
        order = sortperm(spectrum, by=sortby, rev=true)
        spectrum = spectrum[order]
        sectors = sectors[order]
    end
    
    
    for (i, (partial_spectrum, sector)) in enumerate(zip(spectrum, sectors))
        @series begin
            seriestype := :scatter
            label := label(sector)
            n_spectrum = length(partial_spectrum)
            x = n_spectrum == 1 ? [i+1//2] : range(i+1//10,i+9//10, length=length(partial_spectrum))
            return x, partial_spectrum
        end
    end
    
    
    title --> "Entanglement Spectrum"
    legend --> false
    grid --> :xy
    widen --> true
    
    xguide --> "χ = $(dim(domain(s)))"
    xticks --> (1:length(sectors), label.(sectors))
    xtickfonthalign --> :center
    xtick_direction --> :out
    xrotation --> 45
    xlims --> (1, length(sectors)+1)
    
    ylims --> (-Inf, 1+1e-1)
    yscale --> :log10
    
    return ([])
end

@userplot struct TransferPlot{T<:Tuple{InfiniteMPS, InfiniteMPS}}
    args::T
end

@recipe function f(h::TransferPlot; sectors=[], tol=MPSKit.Defaults.tol, num_vals = 20, branchcut=0)
    
    if isempty(sectors)
        mps = h.args[1]
        sectors = [first(TensorKit.sectors(oneunit(virtualspace(mps, 1))))]
    end
    
    for sector in sectors
        spectrum = transfer_spectrum(h.args[1];below=h.args[2], tol=tol, num_vals=num_vals, sector=sector)
        
        @series begin
            yguide --> "r"
            ylims --> (-Inf,1.05)
            
            xguide --> "θ"
            xlims --> (-0.1,2pi)
            xticks --> range(0,2pi,length=7)
            xformatter --> x-> "$(rationalize(x/π, tol=0.05))π"
            xwiden --> true
            seriestype := :scatter
            markershape --> :auto
            label := label(sector)
            return mod2pi.(angle.(spectrum).+branchcut).-branchcut, abs.(spectrum)
        end
    end
    
    title --> "Transfer Spectrum"
    legend --> false
    grid --> :xy
    framestyle --> :zerolines
    
    return nothing
end =#