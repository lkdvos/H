using DrWatson
@quickactivate :Heisenberg

H = su3_heis_ham(3,0,2,3)

function converge_qp(filename_in, sval, filename_out, tols)
    H = su3_heis_ham(3,0,2,3)
    
    local B
    local E
    jldopen(filename_in, "r") do file
        B = file["sval_$sval/B"]
        E = file["sval_$sval/E"]
    end
    
    for tol in tols
        E, B = excitations(H, QuasiparticleAnsatz(toler=10.0^tol), first(B))
        σ = variance(first(B), H)

        println("QP at tol 1e$tol: E = $(real(first(E))) (±$σ).")
        
        jldopen(filename_out, "a+") do file
            file["$tol/B"] = B
            file["$tol/E"] = E
            file["$tol/sigma"] = σ
        end
    end
    
    return nothing
end

function variance_qp(filename)
    H = su3_heis_ham(3,0,2,3)
    
    jldopen(filename, "r+") do file
        
        for tol in keys(file)
            σ = variance(first(file[tol]["B"]), H)
            println("QP at tol 1e$tol: E = $(real(first(file[tol]["E"]))) (±$σ).")
            file[tol]["sigma"] = σ
        end
    end
    return nothing
end

sval = "2.1544346900318822e-6"
filename_in = datadir("sims", "su3300", "excitations.jld2")
filename_out = datadir("sims", "su3300", "exci_conv_$sval.jld2")

tols = [-5.0, -6.0, -7.0, -8.0, -9.0]

converge_qp(filename_in, sval, filename_out, tols)

#variance_qp(filename_out)
