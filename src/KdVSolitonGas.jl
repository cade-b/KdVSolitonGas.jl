module KdVSolitonGas

using OperatorApproximation, LinearAlgebra

include("CauchyInts.jl")
include("AuxiliaryFunctions.jl")
include("RHStructs.jl")
include("Precomputations.jl")

export precompute, deformed_RHP, undeformed_RHP, rhp, deformed_RHP_soliton, undeformed_RHP_soliton, rhp_soliton, deformed_RHP_solitons, undeformed_RHP_solitons, rhp_solitons, get_g, gfunction, get_h, hfunction

function (dp::deformed_RHP)(x,t)
    g = size(dp.bands,1)-1
    ntot = sum(sum(dp.nmat))
    nptot(j) = sum(sum(dp.nmat[1:j,:]))
    cc(j) = (dp.bands[j,1]+dp.bands[j,2])/2
    rr(j) = dp.circ_size*(dp.bands[j,2]-dp.bands[j,1])/2 
    
    Ωvec = dp.Ω(x,t)
    #perform the remaining computations
    Avec = h_coeffs_post(dp.Ah_pre, dp.Bh_pre, Ωvec)
    
    A₂₁ = zeros(ComplexF64,size(dp.A₁₁,1),size(dp.A₂₂,2))
    A₁₂ = zeros(ComplexF64,size(dp.A₂₂,1),size(dp.A₁₁,2))
    rhs₁₂ = zeros(ComplexF64,size(dp.A₁₁,2))
    rhs₂₁ = zeros(ComplexF64,size(dp.A₂₂,2))
    
    for j = 1:g+1
        #get values of 𝔤,𝔥 on circles
        g_vals = map(k->compute_g_post(dp.gvalsp[j][k], x, t, dp.Ω), 1:dp.nmat[j,1])
        h_vals = map(k->compute_h_post(dp.hvalsp[j][k], Avec, Ωvec), 1:dp.nmat[j,1])
        
        if j > (g+1)/2 #Σ₊
            jump_vals = -im*dp.rvals[j,3].*exp.(2x*dp.gridpts[j]-8t*dp.gridpts[j].^3-2g_vals-2h_vals)
            jump1_vals = -im*dp.rvals[j,2]*exp(-Avec[j])
            jump2_vals = -im*dp.rvals[j,1]*exp(Avec[j])
            #println(maximum(abs.(exp.(2x*dp.gridpts[j]-2g_vals[j]))))
        
            #finish top right corner
            A₂₁[nptot(j)-dp.nmat[j,2]+1:nptot(j),:] = -Diagonal(jump2_vals)*dp.Cmats[j,2]
            rhs₂₁[nptot(j)-dp.nmat[j,2]+1:nptot(j)] = jump2_vals.-1

            #finish bottom left corner
            A₁ = -Diagonal(jump_vals)*dp.Cmats[j,3]
            A₂ = -Diagonal(jump1_vals)*dp.Cmats[j,4]
            A₁₂[nptot(j-1)+1:nptot(j),:] = [A₁; A₂]
            rhs₁₂[nptot(j-1)+1:nptot(j)] = [jump_vals; jump1_vals.-1]
        
        else #Σ₋
            jump_vals = im*dp.rvals[j,3].*exp.(2h_vals+2g_vals-2x*dp.gridpts[j]+8t*dp.gridpts[j].^3)
            #println(jump_vals)
            jump1_vals = im*dp.rvals[j,1]*exp(-Avec[j])
            jump2_vals = im*dp.rvals[j,2]*exp(Avec[j])
            #println(maximum(real.(2g_vals[j]-2x*dp.gridpts[j]+8t*dp.gridpts[j].^3)))
            #global jv = jump_vals
            #global jv1 = jump1_vals
            #global jv2 = jump2_vals

            #finish top right corner
            A₁ = -Diagonal(jump_vals)*dp.Cmats[j,1] 
            A₂ = -Diagonal(jump2_vals)*dp.Cmats[j,2]
            A₂₁[nptot(j-1)+1:nptot(j),:] = [A₁; A₂]
            rhs₂₁[nptot(j-1)+1:nptot(j)] = [jump_vals; jump2_vals.-1]

            #finish bottom left corner
            A₁₂[nptot(j)-dp.nmat[j,2]+1:nptot(j),:] = -Diagonal(jump1_vals)*dp.Cmats[j,4]
            rhs₁₂[nptot(j)-dp.nmat[j,2]+1:nptot(j)] = jump1_vals.-1
        end
    end
    
    #build system
    A = [dp.A₁₁ A₂₁; A₁₂ dp.A₂₂]
    #smalleig = minimum(abs.(eigvals(A)))
    rhs = [rhs₂₁; rhs₁₂]
    coeffs = A\rhs
    #println(cond(A))
    
    s1 = zeros(ComplexF64,1,2); s2 = zeros(ComplexF64,1,2)
    for j = 1:g+1
        for k = 1:2
            #Laurent series
            s1[k] -= coeffs[nptot(j-1)+(k-1)*ntot+N₋(dp.nmat[j,1])]*rr(j)
            s2[k] -= coeffs[nptot(j-1)+(k-1)*ntot+N₋(dp.nmat[j,1])]*rr(j)*cc(j)
            s2[k] -= coeffs[nptot(j-1)+(k-1)*ntot+N₋(dp.nmat[j,1])-1]*rr(j)^2

            #Chebyshev series
            s1[k] += coeffs[nptot(j)-dp.nmat[j,2]+(k-1)*ntot+1]*im/2π
            if dp.Chebymat[j,k].kind == 3
                s2[k] += coeffs[nptot(j)-dp.nmat[j,2]+(k-1)*ntot+1]*(dp.bands[j,1]+3dp.bands[j,2])*im/8π
            elseif dp.Chebymat[j,k].kind == 4
                s2[k] += coeffs[nptot(j)-dp.nmat[j,2]+(k-1)*ntot+1]*(3dp.bands[j,1]+dp.bands[j,2])*im/8π
            else
                s2[k] += coeffs[nptot(j)-dp.nmat[j,2]+(k-1)*ntot+1]*(dp.bands[j,1]+dp.bands[j,2])*im/4π
            end
            if dp.Chebymat[j,k].kind != 1
                s2[k] += coeffs[nptot(j)-dp.nmat[j,2]+(k-1)*ntot+2]*(dp.bands[j,2]-dp.bands[j,1])*im/8π
            else
                s2[k] += √2*coeffs[nptot(j)-dp.nmat[j,2]+(k-1)*ntot+2]*(dp.bands[j,2]-dp.bands[j,1])*im/8π
            end
        end
    end
    #println(s1)
    #println(s2)
    2*(s1[1]*s1[2]+s2[1]+s2[2])
end

function (up::undeformed_RHP)(x,t)
    g = size(up.bands,1)-1
    ntot = sum(up.nvec)
    nptot(j) = sum(up.nvec[1:j])

    A₂₁ = zeros(ComplexF64,ntot,ntot); rhs₂₁=zeros(ComplexF64,ntot,1)
    A₁₂ = zeros(ComplexF64,ntot,ntot); rhs₁₂=zeros(ComplexF64,ntot,1)
    for j = 1:g+1
        if j > (g+1)/2
            jump_vals = -im*up.rvals[j].*exp.(-2x*up.gridpts[j]+8t*up.gridpts[j].^3)
            A₂ = -Diagonal(jump_vals)*up.Cmats[j,1]
            A₂₁[nptot(j-1)+1:nptot(j),:] = A₂
            rhs₂₁[nptot(j-1)+1:nptot(j)] = jump_vals
        else
            jump_vals = im*up.rvals[j].*exp.(2x*up.gridpts[j]-8t*up.gridpts[j].^3)
            A₂ = -Diagonal(jump_vals)*up.Cmats[j,2]
            A₁₂[nptot(j-1)+1:nptot(j),:] = A₂
            rhs₁₂[nptot(j-1)+1:nptot(j)] = jump_vals
        end
    end
        
    A = [up.A₁₁ A₂₁; A₁₂ up.A₂₂]
    #smalleig = minimum(abs.(eigvals(A)))
    rhs = [rhs₂₁; rhs₁₂]
    solvec = A\rhs
    
    s1 = zeros(ComplexF64,1,2); s2 = zeros(ComplexF64,1,2)
    for j = 1:g+1
        for k = 1:2
            s1[k] += solvec[nptot(j-1)+(k-1)*ntot+1]*im/2π
            if up.Chebymat[j,k].kind == 3
                s2[k] += solvec[nptot(j-1)+(k-1)*ntot+1]*(up.bands[j,1]+3up.bands[j,2])*im/8π
            elseif up.Chebymat[j,k].kind == 4
                s2[k] += solvec[nptot(j-1)+(k-1)*ntot+1]*(3up.bands[j,1]+up.bands[j,2])*im/8π
            else
                s2[k] += solvec[nptot(j-1)+(k-1)*ntot+1]*(up.bands[j,1]+up.bands[j,2])*im/4π
            end
            if up.Chebymat[j,k].kind != 1
                s2[k] += solvec[nptot(j-1)+(k-1)*ntot+2]*(up.bands[j,2]-up.bands[j,1])*im/8π
            else
                s2[k] += √2*solvec[nptot(j-1)+(k-1)*ntot+2]*(up.bands[j,2]-up.bands[j,1])*im/8π
            end
        end
    end
    2*(s1[1]*s1[2]+s2[1]+s2[2])     
end

function (rp::rhp)(x,t)
    if x>4t*rp.dp.bands[end,2]^2
        return rp.up(x,t)#rhp = undeformed_RHP(rp.bands,rp.A₁₁p,rp.A₂₂p,rp.Cmatsp,rp.Chebymatp,rp.nmat[:,2],rp.gridmat[:,2],rp.rvals[:,1])     
    else
        return rp.dp(x,t)#rhp = deformed_RHP(rp.bands,rp.A₁₁,rp.A₂₂,rp.Cmats,rp.Ω,rp.Chebymat,rp.nmat,rp.gridmat[:,1],rp.gvalsp,rp.hvalsp,rp.Ah_pre,rp.Bh_pre,rp.rvals)
    end
end

### Add distinguished soliton ####################################################################################################################

function (dp::deformed_RHP_soliton)(x,t; flip_tol = 10., pole_circ = 0.001, flip=nothing, max_deriv_terms=20, verbose=true)
    g = size(dp.bands,1)-1
    ntot = sum(sum(dp.nmat))
    nptot(j) = sum(sum(dp.nmat[1:j,:]))
    cc(j) = (dp.bands[j,1]+dp.bands[j,2])/2
    rr(j) = dp.circ_size*(dp.bands[j,2]-dp.bands[j,1])/2 
    
    Ωvec = dp.Ω(x,t)
    #perform the remaining computations
    Avec = h_coeffs_post(dp.Ah_pre, dp.Bh_pre, Ωvec)

    gk0 = compute_g_post(dp.gvalsp[end][1], x, t, dp.Ω)
    hk0 = compute_h_post(dp.hvalsp[end][1], Avec, Ωvec)
    res_cond = dp.χ*exp(-2x*dp.κ₀+8t*dp.κ₀^3+2gk0+2hk0)
    #check modulus of residue condition
    if flip === nothing
        flip = abs(res_cond)>flip_tol #check whether or not to flip problem
    end
    #println("Size of residue is: ",abs(res_cond))
    #println(flip)
    gf(j) = v -> j == 1 ? v : 1. #function to apply flip
    
    A₂₁ = zeros(ComplexF64,size(dp.A₁₁,1),size(dp.A₂₂,2))
    A₁₂ = zeros(ComplexF64,size(dp.A₂₂,1),size(dp.A₁₁,2))
    rhs₁₂ = zeros(ComplexF64,size(dp.A₁₁,2))
    rhs₂₁ = zeros(ComplexF64,size(dp.A₂₂,2))
    
    for j = 1:g+1
        #get values of 𝔤,𝔥 on circles
        g_vals = map(k->compute_g_post(dp.gvalsp[j][k], x, t, dp.Ω), 1:dp.nmat[j,1])
        h_vals = map(k->compute_h_post(dp.hvalsp[j][k], Avec, Ωvec), 1:dp.nmat[j,1])
        v_vals_circ = gf(flip).(dp.vvals[j,2])
        v_vals_int = gf(flip).(dp.vvals[j,1])

        if j > (g+1)/2 #Σ₊
            jump_vals = -im*dp.rvals[j,3].*exp.(2x*dp.gridpts[j]-8t*dp.gridpts[j].^3-2g_vals-2h_vals)./v_vals_circ
            jump1_vals = -im*dp.rvals[j,2]*exp(-Avec[j])./v_vals_int
            jump2_vals = -im*dp.rvals[j,1]*exp(Avec[j]).*v_vals_int
        
            #finish top right corner
            A₂₁[nptot(j)-dp.nmat[j,2]+1:nptot(j),:] = -Diagonal(jump2_vals)*dp.Cmats[j,2]
            rhs₂₁[nptot(j)-dp.nmat[j,2]+1:nptot(j)] = jump2_vals

            #finish bottom left corner
            A₁ = -Diagonal(jump_vals)*dp.Cmats[j,3]
            A₂ = -Diagonal(jump1_vals)*dp.Cmats[j,4]
            A₁₂[nptot(j-1)+1:nptot(j),:] = [A₁; A₂]
            rhs₁₂[nptot(j-1)+1:nptot(j)] = [jump_vals; jump1_vals]
        
        else #Σ₋
            jump_vals = im*dp.rvals[j,3].*exp.(2h_vals+2g_vals-2x*dp.gridpts[j]+8t*dp.gridpts[j].^3).*v_vals_circ
            jump1_vals = im*dp.rvals[j,1]*exp(-Avec[j])./v_vals_int
            jump2_vals = im*dp.rvals[j,2]*exp(Avec[j]).*v_vals_int

            #finish top right corner
            A₁ = -Diagonal(jump_vals)*dp.Cmats[j,1] 
            A₂ = -Diagonal(jump2_vals)*dp.Cmats[j,2]
            A₂₁[nptot(j-1)+1:nptot(j),:] = [A₁; A₂]
            rhs₂₁[nptot(j-1)+1:nptot(j)] = [jump_vals; jump2_vals]

            #finish bottom left corner
            A₁₂[nptot(j)-dp.nmat[j,2]+1:nptot(j),:] = -Diagonal(jump1_vals)*dp.Cmats[j,4]
            rhs₁₂[nptot(j)-dp.nmat[j,2]+1:nptot(j)] = jump1_vals
        end
    end
    
    #build system
    A = [dp.A₁₁ A₂₁; A₁₂ dp.A₂₂]
    rhs = [dp.rhs₁₁ rhs₂₁; rhs₁₂ dp.rhs₂₂]
    
    mat_coeffs = A\rhs
    #println("Condition number for continuous problem: ",cond(A))

    #functon to get out matrix solution
    function eval_sol(z)
        S = zeros(ComplexF64,2,2)
        for j = 1:g+1
            for k = 1:2
                #sum Chebyshev series
                polys = dp.Chebymat[j,k]
                S[:,k] += transpose(CauchyInterval(z,polys,dp.nmat[j,2]-1)*mat_coeffs[nptot(j)-dp.nmat[j,2]+(k-1)*ntot+1:nptot(j)+(k-1)*ntot,:])

                #sum Laurent series
                S[:,k] += transpose(CauchyVec(z,dp.nmat[j,1],cc(j),rr(j))*mat_coeffs[nptot(j-1)+(k-1)*ntot+1:nptot(j-1)+dp.nmat[j,1]+(k-1)*ntot,:])
            end
        end
        return S+I(2)
    end
    #zzz = 1+√2; r2(z) = 2*sqrt(3-z)*sqrt(z-2) 
    #println("Check interval jump condition: ",eval_sol(zzz+1e-14im)-eval_sol(zzz-1e-14im)*[0 -im*exp(-Avec[2])/r2(zzz); -im*r2(zzz)*exp(Avec[2]) 0.])
    #zz = 2.5+0.625im
    #println("Check circle jump condition: ",eval_sol(zz-1e-14im)-eval_sol(zz+1e-14im)*[1 -im*exp(2x*zz-8t*zz^3-2gk25-2hk25)/r2(zz); 0 1.])
    Vr1 = eval_sol(dp.κ₀)
    #println(Vr1)
    #println("Continuous solution determinant: ",det(Vr1))
    Vr2 = eval_sol(-dp.κ₀)

    #get derivative terms
    function get_der(κ;epsc=pole_circ, max_terms=max_deriv_terms)
        d11 = BasisExpansion(z->eval_sol(z)[1,1],Laurent(PeriodicMappedCircle(κ,epsc)),max_terms)
        d12 = BasisExpansion(z->eval_sol(z)[1,2],Laurent(PeriodicMappedCircle(κ,epsc)),max_terms)
        d21 = BasisExpansion(z->eval_sol(z)[2,1],Laurent(PeriodicMappedCircle(κ,epsc)),max_terms)
        d22 = BasisExpansion(z->eval_sol(z)[2,2],Laurent(PeriodicMappedCircle(κ,epsc)),max_terms)
        tol_met = norm([norm(d11.c[end-4:end]) norm(d12.c[end-4:end]); norm(d21.c[end-4:end]) norm(d22.c[end-4:end])])
        if  tol_met > 1e-11 && verbose
            @warn "Last terms in Δ Taylor expansion are not small. Tolerance met: "*string(tol_met)
        end
        f1p = [CauchyTransform()*(Derivative(1)*d11) CauchyTransform()*(Derivative(1)*d12);
        CauchyTransform()*(Derivative(1)*d21) CauchyTransform()*(Derivative(1)*d22)]
        return f1p(κ)
    end
    Vr1p = get_der(dp.κ₀)
    Vr2p = get_der(-dp.κ₀)
    #println(det(Vr2))
    #println(abs.(sum(Vr1,dims=1)-reverse(sum(Vr2,dims=1))))

    #Solve residue problem
    poles = [dp.κ₀, -dp.κ₀]
    if flip == 1
        Rp = Vr1*[0 4dp.κ₀^2/res_cond; 0 0]*inv(Vr1-Vr1p*[0 4dp.κ₀^2/res_cond; 0 0])#*inv(Vr1)*inv(I(2)-Vr1p*[0 4dp.κ₀^2/res_cond; 0 0]*inv(Vr1))
        Rm = Vr2*[0 0; -4dp.κ₀^2/res_cond 0]*inv(Vr2-Vr2p*[0 0; -4dp.κ₀^2/res_cond 0])#*inv(Vr2)*inv(I(2)-Vr2p*[0 0; -4dp.κ₀^2/res_cond 0]*inv(Vr2))
    else
        Rp = Vr1*[0 0; res_cond 0]*inv(Vr1-Vr1p*[0 0; res_cond 0])#*inv(Vr1)*inv(I(2)-Vr1p*[0 0; res_cond 0]*inv(Vr1))
        #println(norm(I(2)-Vr1p*[0 0; res_cond 0]*inv(Vr1)))
        Rm = Vr2*[0 -res_cond; 0 0]*inv(Vr2-Vr2p*[0 -res_cond; 0 0])#*inv(Vr2)*inv(I(2)-Vr2p*[0 -res_cond; 0 0]*inv(Vr2))
    end
    R = [Rp, Rm]

    #If residue is too large, flip the problem
    if norm(Rp)>10*flip_tol && flip == 0 && verbose
        @warn "Continuous solution is large. Flipping problem."
        return dp(x,t; flip_tol = flip_tol, pole_circ = pole_circ, flip=1, max_deriv_terms=max_deriv_terms)
    end
    #println(poles)
    #println(Rp)
    pole_rhp = RHP([],[],poles,R)
    rhsolver = RHSolver(pole_rhp);
    res_sol = rhsolver([1 1], 300)
    #some_val = moment(res_sol[1],0)*1im/(2pi)
    #println(some_val)

    #get behavior at infinity for residue solution
    res_a1 = [moment(res_sol[1],0) moment(res_sol[2],0)]*1im/(2pi)
    res_a2 = [moment(res_sol[1],1) moment(res_sol[2],1)]*1im/(2pi)
    #println(res_a1)
    #println(res_a2)

    #get behavior at infinity for gas solution
    cont_a1 = zeros(ComplexF64,2,2); cont_a2 = zeros(ComplexF64,2,2)
    for j = 1:g+1
        for k = 1:2
            #Laurent series
            #println(mat_coeffs[nptot(j-1)+(k-1)*ntot+N₋(dp.nmat[j,1]),:]*rr(j))
            cont_a1[:,k] -= mat_coeffs[nptot(j-1)+(k-1)*ntot+N₋(dp.nmat[j,1]),:]*rr(j)
            cont_a2[:,k] -= mat_coeffs[nptot(j-1)+(k-1)*ntot+N₋(dp.nmat[j,1]),:]*rr(j)*cc(j)
            cont_a2[:,k] -= mat_coeffs[nptot(j-1)+(k-1)*ntot+N₋(dp.nmat[j,1])-1,:]*rr(j)^2

            #Chebyshev series
            cont_a1[:,k] += mat_coeffs[nptot(j)-dp.nmat[j,2]+(k-1)*ntot+1,:]*im/2π
            if dp.Chebymat[j,k].kind == 3
                cont_a2[:,k] += mat_coeffs[nptot(j)-dp.nmat[j,2]+(k-1)*ntot+1,:]*(dp.bands[j,1]+3dp.bands[j,2])*im/8π
            elseif dp.Chebymat[j,k].kind == 4
                cont_a2[:,k] += mat_coeffs[nptot(j)-dp.nmat[j,2]+(k-1)*ntot+1,:]*(3dp.bands[j,1]+dp.bands[j,2])*im/8π
            else
                cont_a2[:,k] += mat_coeffs[nptot(j)-dp.nmat[j,2]+(k-1)*ntot+1,:]*(dp.bands[j,1]+dp.bands[j,2])*im/4π
            end
            if dp.Chebymat[j,k].kind != 1
                cont_a2[:,k] += mat_coeffs[nptot(j)-dp.nmat[j,2]+(k-1)*ntot+2,:]*(dp.bands[j,2]-dp.bands[j,1])*im/8π
            else
                cont_a2[:,k] += √2*mat_coeffs[nptot(j)-dp.nmat[j,2]+(k-1)*ntot+2,:]*(dp.bands[j,2]-dp.bands[j,1])*im/8π
            end
        end
    end
    #println(cont_a1)
    #println(cont_a2)
    #println("Check correct residue for continuous problem: ",cont_a1-(eval_sol(1e6)-I(2))*1e6)
    #println("Check correct 2nd order for continuous problem: ",cont_a2-(eval_sol(1e4)-I(2)-cont_a1*1e-4)*1e8)

    #get behavior at infinity for product
    s1 = res_a1+[1 1]*cont_a1
    s2 = res_a2+[1 1]*cont_a2+res_a1*cont_a1
    #println(s1)
    #println(s2)
    #res_recov = 2*(res_a1[1]*res_a1[2]+res_a2[1]+res_a2[2])
    #cont1 = [1 1]*cont_a1
    #cont2 = [1 1]*cont_a2
    #cont_recov = 2*(cont1[1]*cont1[2]+cont2[1]+cont2[2])

    return 2*(s1[1]*s1[2]+s2[1]+s2[2])
end

function (up::undeformed_RHP_soliton)(x,t; flip_tol = 10., pole_circ = 0.001, flip=nothing, max_deriv_terms=20, verbose=true)
    g = size(up.bands,1)-1
    ntot = sum(up.nvec)
    nptot(j) = sum(up.nvec[1:j])

    res_cond = up.χ*exp(-2x*up.κ₀+8t*up.κ₀^3)
    #println(abs(res_cond))
    if flip === nothing
        flip = abs(res_cond)>flip_tol #check whether or not to flip problem
    end
    #println(flip)
    gf(j) = v -> j == 1 ? v : 1. #function to apply flip

    A₂₁ = zeros(ComplexF64,ntot,ntot); rhs₂₁=zeros(ComplexF64,ntot,1)
    A₁₂ = zeros(ComplexF64,ntot,ntot); rhs₁₂=zeros(ComplexF64,ntot,1)
    
    for j = 1:g+1
        vvals = gf(flip).(up.vvals[j])
        #println(vvals)
        if j > (g+1)/2
            jump_vals = -im*up.rvals[j].*exp.(-2x*up.gridpts[j]+8t*up.gridpts[j].^3).*vvals
            A₂ = -Diagonal(jump_vals)*up.Cmats[j,1]
            A₂₁[nptot(j-1)+1:nptot(j),:] = A₂
            rhs₂₁[nptot(j-1)+1:nptot(j)] = jump_vals
        else
            jump_vals = im*up.rvals[j].*exp.(2x*up.gridpts[j]-8t*up.gridpts[j].^3)./vvals
            A₂ = -Diagonal(jump_vals)*up.Cmats[j,2]
            A₁₂[nptot(j-1)+1:nptot(j),:] = A₂
            rhs₁₂[nptot(j-1)+1:nptot(j)] = jump_vals
        end
    end
        
    A = [up.A₁₁ A₂₁; A₁₂ up.A₂₂]
    rhs = [zeros(ntot) rhs₂₁; rhs₁₂ zeros(ntot)]
    mat_coeffs = A\rhs
    #println(cond(A))

    #functon to get out matrix solution
    function eval_sol(z;flag=0)
        S = zeros(ComplexF64,2,2)
        for j = 1:g+1
            for k = 1:2
                #sum Chebyshev series
                polys = up.Chebymat[j,k]
                S[:,k] += transpose(CauchyInterval(z,polys,up.nvec[j]-1;flag=flag)*mat_coeffs[nptot(j-1)+(k-1)*ntot+1:nptot(j)+(k-1)*ntot,:])
            end
        end
        return S+I(2)
    end
    Vr1 = eval_sol(up.κ₀)
    #println("Continuous solution determinant: ",det(Vr1))
    #println(Vr1)
    Vr2 = eval_sol(-up.κ₀)

    #get derivative terms
    #Vr1p = (eval_sol(up.κ₀+1e-8)-eval_sol(up.κ₀-1e-8))/(2e-8)
    #Vr2p = (eval_sol(-up.κ₀+1e-8)-eval_sol(-up.κ₀-1e-8))/(2e-8)
    function get_der(κ;epsc=pole_circ, max_terms=max_deriv_terms)
        d11 = BasisExpansion(z->eval_sol(z)[1,1],Laurent(PeriodicMappedCircle(κ,epsc)),max_terms)
        d12 = BasisExpansion(z->eval_sol(z)[1,2],Laurent(PeriodicMappedCircle(κ,epsc)),max_terms)
        d21 = BasisExpansion(z->eval_sol(z)[2,1],Laurent(PeriodicMappedCircle(κ,epsc)),max_terms)
        d22 = BasisExpansion(z->eval_sol(z)[2,2],Laurent(PeriodicMappedCircle(κ,epsc)),max_terms)
        tol_met = norm([norm(d11.c[end-4:end]) norm(d12.c[end-4:end]); norm(d21.c[end-4:end]) norm(d22.c[end-4:end])])
        if  tol_met > 1e-11 && verbose
            @warn "Last terms in Δ Taylor expansion are not small. Tolerance met: "*string(tol_met)
        end
        f1p = [CauchyTransform()*(Derivative(1)*d11) CauchyTransform()*(Derivative(1)*d12);
        CauchyTransform()*(Derivative(1)*d21) CauchyTransform()*(Derivative(1)*d22)]
        return f1p(κ)
    end
    Vr1p = get_der(up.κ₀)
    Vr2p = get_der(-up.κ₀)

    #Solve residue problem
    poles = [up.κ₀, -up.κ₀]
    if flip == 1
        Rp = Vr1*[0 4up.κ₀^2 ./ res_cond; 0 0]*inv(Vr1-Vr1p*[0 4up.κ₀^2 ./ res_cond; 0 0])#*inv(Vr1)*inv(I(2)-Vr1p*[0 4up.κ₀^2 ./ res_cond; 0 0]*inv(Vr1))
        Rm = Vr2*[0 0; -4up.κ₀^2 ./res_cond 0]*inv(Vr2-Vr2p*[0 0; -4up.κ₀^2 ./res_cond 0])#*inv(Vr2)*inv(I(2)-Vr2p*[0 0; -4up.κ₀^2 ./res_cond 0]*inv(Vr2))
    else
        Rp = Vr1*[0 0; res_cond 0]*inv(Vr1-Vr1p*[0 0; res_cond 0])#*inv(Vr1)*inv(I(2)-Vr1p*[0 0; res_cond 0]*inv(Vr1))
        Rm = Vr2*[0 -res_cond; 0 0]*inv(Vr2-Vr2p*[0 -res_cond; 0 0])#*inv(Vr2)*inv(I(2)-Vr2p*[0 -res_cond; 0 0]*inv(Vr2))
    end
    #println(Rp[1,2])
    #println(Rm[2,1])
    R = [Rp, Rm]
    #If residue is too large, flip the problem
    if norm(Rp)>10*flip_tol && flip == 0 && verbose
        @warn "Continuous solution is large. Flipping problem."
        return up(x,t; flip_tol = flip_tol, pole_circ = pole_circ, flip=1, max_deriv_terms=max_deriv_terms)
    end

    pole_rhp = RHP([],[],poles,R)
    rhsolver = RHSolver(pole_rhp);
    res_sol = rhsolver([1 1], 300)
    #some_val = moment(res_sol[1],0)*1im/(2pi)
    #println(some_val)

    #get behavior at infinity for residue solution
    res_a1 = [moment(res_sol[1],0) moment(res_sol[2],0)]*1im/(2pi)
    res_a2 = [moment(res_sol[1],1) moment(res_sol[2],1)]*1im/(2pi)
    #println(res_a1)
    #println(res_a2)
    
    #get behavior at infinity for gas solution
    cont_a1 = zeros(ComplexF64,2,2); cont_a2 = zeros(ComplexF64,2,2)
    for j = 1:g+1
        for k = 1:2
            #Chebyshev series
            cont_a1[:,k] += mat_coeffs[nptot(j-1)+(k-1)*ntot+1,:]*im/2π
            if up.Chebymat[j,k].kind == 3
                cont_a2[:,k] += mat_coeffs[nptot(j-1)+(k-1)*ntot+1,:]*(up.bands[j,1]+3up.bands[j,2])*im/8π
            elseif up.Chebymat[j,k].kind == 4
                cont_a2[:,k] += mat_coeffs[nptot(j-1)+(k-1)*ntot+1,:]*(3up.bands[j,1]+up.bands[j,2])*im/8π
            else
                cont_a2[:,k] += mat_coeffs[nptot(j-1)+(k-1)*ntot+1,:]*(up.bands[j,1]+up.bands[j,2])*im/4π
            end
            if up.Chebymat[j,k].kind != 1
                cont_a2[:,k] += mat_coeffs[nptot(j-1)+(k-1)*ntot+2,:]*(up.bands[j,2]-up.bands[j,1])*im/8π
            else
                cont_a2[:,k] += √2*mat_coeffs[nptot(j-1)+(k-1)*ntot+2,:]*(up.bands[j,2]-up.bands[j,1])*im/8π
            end
        end
    end
    
    #get behavior at infinity for product
    s1 = res_a1+[1 1]*cont_a1
    s2 = res_a2+[1 1]*cont_a2+res_a1*cont_a1
    #println(s1)
    #println(s2)

    return 2*(s1[1]*s1[2]+s2[1]+s2[2])
end

function (rp::rhp_soliton)(x,t; flip_tol = 10., pole_circ = 0.001, flip=nothing, max_deriv_terms=25, verbose=true)
    if x>4t*rp.dp.bands[end,2]^2
        #println("Solving undeformed problem")
        return rp.up(x,t; flip_tol = flip_tol, pole_circ = pole_circ, flip=flip, max_deriv_terms=max_deriv_terms, verbose=verbose)    
    else
        #println("Solving deformed problem")
        return rp.dp(x,t; flip_tol = flip_tol, pole_circ = pole_circ, flip=flip, max_deriv_terms=max_deriv_terms, verbose=verbose)
    end
end

### add several solitons #############################################################################################
function (dp::deformed_RHP_solitons)(x,t; flip_tol = 10., pole_circ = 0.001, flips=nothing, max_deriv_terms=20, verbose=true)
    g = size(dp.bands,1)-1
    ntot = sum(sum(dp.nmat))
    nptot(j) = sum(sum(dp.nmat[1:j,:]))
    cc(j) = (dp.bands[j,1]+dp.bands[j,2])/2
    rr(j) = dp.circ_size*(dp.bands[j,2]-dp.bands[j,1])/2 
    
    Ωvec = dp.Ω(x,t)
    #perform the remaining computations
    Avec = h_coeffs_post(dp.Ah_pre, dp.Bh_pre, Ωvec)

    #check modulus of residue condition
    gks = compute_g_post.(dp.gvalsp[end], x, t, dp.Ω)
    #println(gk0)
    #println(gk0m)
    hks = map(k->compute_h_post(dp.hvalsp[end][k], Avec, Ωvec), 1:length(dp.κvec))
    res_conds = dp.χvec.*exp.(-2x*dp.κvec.+8t*dp.κvec.^3+2gks+2hks)
    if flips === nothing
        flips = abs.(res_conds).>flip_tol #check whether or not to flip pole
    end
    flipped_poles = dp.κvec[flips]
    
    A₂₁ = zeros(ComplexF64,size(dp.A₁₁,1),size(dp.A₂₂,2))
    A₁₂ = zeros(ComplexF64,size(dp.A₂₂,1),size(dp.A₁₁,2))
    rhs₁₂ = zeros(ComplexF64,size(dp.A₁₁,2))
    rhs₂₁ = zeros(ComplexF64,size(dp.A₂₂,2))
    
    for j = 1:g+1
        #get values of 𝔤,𝔥 on circles
        g_vals = map(k->compute_g_post(dp.gvalsp[j][k], x, t, dp.Ω), 1:dp.nmat[j,1])
        h_vals = map(k->compute_h_post(dp.hvalsp[j][k], Avec, Ωvec), 1:dp.nmat[j,1])
        v_vals_circ = map(z->v(z,flipped_poles)^2, dp.gridpts[j,1])
        v_vals_int = map(z->v(z,flipped_poles)^2, dp.gridpts[j,2])

        if j > (g+1)/2 #Σⱼ₊
            jump_vals = -im*dp.rvals[j,3].*exp.(2x*dp.gridpts[j]-8t*dp.gridpts[j].^3-2g_vals-2h_vals)./v_vals_circ
            jump1_vals = -im*dp.rvals[j,2]*exp(-Avec[j])./v_vals_int
            jump2_vals = -im*dp.rvals[j,1]*exp(Avec[j]).*v_vals_int
        
            #finish top right corner
            A₂₁[nptot(j)-dp.nmat[j,2]+1:nptot(j),:] = -Diagonal(jump2_vals)*dp.Cmats[j,2]
            rhs₂₁[nptot(j)-dp.nmat[j,2]+1:nptot(j)] = jump2_vals

            #finish bottom left corner
            A₁ = -Diagonal(jump_vals)*dp.Cmats[j,3]
            A₂ = -Diagonal(jump1_vals)*dp.Cmats[j,4]
            A₁₂[nptot(j-1)+1:nptot(j),:] = [A₁; A₂]
            rhs₁₂[nptot(j-1)+1:nptot(j)] = [jump_vals; jump1_vals]
        
        else #Σⱼ₋
            jump_vals = im*dp.rvals[j,3].*exp.(2h_vals+2g_vals-2x*dp.gridpts[j]+8t*dp.gridpts[j].^3).*v_vals_circ
            jump1_vals = im*dp.rvals[j,1]*exp(-Avec[j])./v_vals_int
            jump2_vals = im*dp.rvals[j,2]*exp(Avec[j]).*v_vals_int

            #finish top right corner
            A₁ = -Diagonal(jump_vals)*dp.Cmats[j,1] 
            A₂ = -Diagonal(jump2_vals)*dp.Cmats[j,2]
            A₂₁[nptot(j-1)+1:nptot(j),:] = [A₁; A₂]
            rhs₂₁[nptot(j-1)+1:nptot(j)] = [jump_vals; jump2_vals]

            #finish bottom left corner
            A₁₂[nptot(j)-dp.nmat[j,2]+1:nptot(j),:] = -Diagonal(jump1_vals)*dp.Cmats[j,4]
            rhs₁₂[nptot(j)-dp.nmat[j,2]+1:nptot(j)] = jump1_vals
        end
    end
    
    #build system
    A = [dp.A₁₁ A₂₁; A₁₂ dp.A₂₂]
    rhs = [dp.rhs₁₁ rhs₂₁; rhs₁₂ dp.rhs₂₂]
    
    mat_coeffs = A\rhs
    #println("Condition number for continuous problem: ",cond(A))

    #functon to get out matrix solution
    function eval_sol(z)
        S = zeros(ComplexF64,2,2)
        for j = 1:g+1
            for k = 1:2
                #sum Chebyshev series
                polys = dp.Chebymat[j,k]
                S[:,k] += transpose(CauchyInterval(z,polys,dp.nmat[j,2]-1)*mat_coeffs[nptot(j)-dp.nmat[j,2]+(k-1)*ntot+1:nptot(j)+(k-1)*ntot,:])

                #sum Laurent series
                S[:,k] += transpose(CauchyVec(z,dp.nmat[j,1],cc(j),rr(j))*mat_coeffs[nptot(j-1)+(k-1)*ntot+1:nptot(j-1)+dp.nmat[j,1]+(k-1)*ntot,:])
            end
        end
        return S+I(2)
    end
    #zzz = 1+√2; r2(z) = 2*sqrt(3-z)*sqrt(z-2) 
    #println("Check interval jump condition: ",eval_sol(zzz+1e-14im)-eval_sol(zzz-1e-14im)*[0 -im*exp(-Avec[2])/r2(zzz); -im*r2(zzz)*exp(Avec[2]) 0.])
    #zz = 2.5+0.625im
    #println("Check circle jump condition: ",eval_sol(zz-1e-14im)-eval_sol(zz+1e-14im)*[1 -im*exp(2x*zz-8t*zz^3-2gk25-2hk25)/r2(zz); 0 1.])
    poles = vcat(dp.κvec, -dp.κvec)
    Vr = eval_sol.(poles)
    #println("Norm of continuous solution at pole: ",norm(Vr[1]))

    #get derivative terms
    function get_der(κ;epsc=pole_circ, max_terms=max_deriv_terms)
        d11 = BasisExpansion(z->eval_sol(z)[1,1],Laurent(PeriodicMappedCircle(κ,epsc)),max_terms)
        d12 = BasisExpansion(z->eval_sol(z)[1,2],Laurent(PeriodicMappedCircle(κ,epsc)),max_terms)
        d21 = BasisExpansion(z->eval_sol(z)[2,1],Laurent(PeriodicMappedCircle(κ,epsc)),max_terms)
        d22 = BasisExpansion(z->eval_sol(z)[2,2],Laurent(PeriodicMappedCircle(κ,epsc)),max_terms)
        tol_met = norm([norm(d11.c[end-4:end]) norm(d12.c[end-4:end]); norm(d21.c[end-4:end]) norm(d22.c[end-4:end])])
        if  tol_met > 1e-11 && verbose
            @warn "Last terms in Δ Taylor expansion are not small. Tolerance met: "*string(tol_met)
        end
        f1p = [CauchyTransform()*(Derivative(1)*d11) CauchyTransform()*(Derivative(1)*d12);
        CauchyTransform()*(Derivative(1)*d21) CauchyTransform()*(Derivative(1)*d22)]
        return f1p(κ)
    end
    Vrp = get_der.(poles)
    #println(det(Vr2))
    #println(abs.(sum(Vr1,dims=1)-reverse(sum(Vr2,dims=1))))

    #Solve residue problem
    R = Vector{Matrix{ComplexF64}}(undef,length(poles))
    for (j,κ) in enumerate(dp.κvec)
        k = j+length(dp.κvec)
        if flips[j] == 1
            R[j] = Vr[j]*[0 1/(res_conds[j]*dv(κ,flipped_poles)^2); 0 0]*inv(Vr[j]-Vrp[j]*[0 1/(res_conds[j]*dv(κ,flipped_poles)^2); 0 0])#*inv(Vr[j])*inv(I(2)-Vrp[j]*[0 1/(res_conds[j]*dv(κ,flipped_poles)^2); 0 0]*inv(Vr[j]))
            R[k] = Vr[k]*[0 0; -1/(res_conds[j]*dv(κ,flipped_poles)^2) 0]*inv(Vr[k]-Vrp[k]*[0 0; -1/(res_conds[j]*dv(κ,flipped_poles)^2) 0])#*inv(Vr[k])*inv(I(2)-Vrp[k]*[0 0; -1/(res_conds[j]*dv(κ,flipped_poles)^2) 0]*inv(Vr[k]))
        else
            R[j] = Vr[j]*[0 0; res_conds[j]*v(κ,flipped_poles)^2 0]*inv(Vr[j]-Vrp[j]*[0 0; res_conds[j]*v(κ,flipped_poles)^2 0])#*inv(Vr[j])*inv(I(2)-Vrp[j]*[0 0; res_conds[j]*v(κ,flipped_poles)^2 0]*inv(Vr[j]))
            R[k] = Vr[k]*[0 -res_conds[j]*v(κ,flipped_poles)^2; 0 0]*inv(Vr[k]-Vrp[k]*[0 -res_conds[j]*v(κ,flipped_poles)^2; 0 0])#*inv(Vr[k])*inv(I(2)-Vrp[k]*[0 -res_conds[j]*v(κ,flipped_poles)^2; 0 0]*inv(Vr[k]))
        end
    end
    #println(poles)
    #println(R[2])
    pole_rhp = RHP([],[],poles,R)
    rhsolver = RHSolver(pole_rhp);
    res_sol = rhsolver([1 1], 300)
    #some_val = moment(res_sol[1],0)*1im/(2pi)
    #println(some_val)

    #get behavior at infinity for residue solution
    res_a1 = [moment(res_sol[1],0) moment(res_sol[2],0)]*1im/(2pi)
    res_a2 = [moment(res_sol[1],1) moment(res_sol[2],1)]*1im/(2pi)
    #println(res_a1)
    #println(res_a2)

    #get behavior at infinity for gas solution
    cont_a1 = zeros(ComplexF64,2,2); cont_a2 = zeros(ComplexF64,2,2)
    for j = 1:g+1
        for k = 1:2
            #Laurent series
            #println(mat_coeffs[nptot(j-1)+(k-1)*ntot+N₋(dp.nmat[j,1]),:]*rr(j))
            cont_a1[:,k] -= mat_coeffs[nptot(j-1)+(k-1)*ntot+N₋(dp.nmat[j,1]),:]*rr(j)
            cont_a2[:,k] -= mat_coeffs[nptot(j-1)+(k-1)*ntot+N₋(dp.nmat[j,1]),:]*rr(j)*cc(j)
            cont_a2[:,k] -= mat_coeffs[nptot(j-1)+(k-1)*ntot+N₋(dp.nmat[j,1])-1,:]*rr(j)^2

            #Chebyshev series
            cont_a1[:,k] += mat_coeffs[nptot(j)-dp.nmat[j,2]+(k-1)*ntot+1,:]*im/2π
            if dp.Chebymat[j,k].kind == 3
                cont_a2[:,k] += mat_coeffs[nptot(j)-dp.nmat[j,2]+(k-1)*ntot+1,:]*(dp.bands[j,1]+3dp.bands[j,2])*im/8π
            elseif dp.Chebymat[j,k].kind == 4
                cont_a2[:,k] += mat_coeffs[nptot(j)-dp.nmat[j,2]+(k-1)*ntot+1,:]*(3dp.bands[j,1]+dp.bands[j,2])*im/8π
            else
                cont_a2[:,k] += mat_coeffs[nptot(j)-dp.nmat[j,2]+(k-1)*ntot+1,:]*(dp.bands[j,1]+dp.bands[j,2])*im/4π
            end
            if dp.Chebymat[j,k].kind != 1
                cont_a2[:,k] += mat_coeffs[nptot(j)-dp.nmat[j,2]+(k-1)*ntot+2,:]*(dp.bands[j,2]-dp.bands[j,1])*im/8π
            else
                cont_a2[:,k] += √2*mat_coeffs[nptot(j)-dp.nmat[j,2]+(k-1)*ntot+2,:]*(dp.bands[j,2]-dp.bands[j,1])*im/8π
            end
        end
    end
    #println("Check correct residue for continuous problem: ",cont_a1-(eval_sol(1e6)-I(2))*1e6)
    #println("Check correct 2nd order for continuous problem: ",cont_a2-(eval_sol(1e4)-I(2)-cont_a1*1e-4)*1e8)
    #println(cont_a1)
    #println(cont_a2)

    #get behavior at infinity for product
    s1 = res_a1+[1 1]*cont_a1
    s2 = res_a2+[1 1]*cont_a2+res_a1*cont_a1

    return 2*(s1[1]*s1[2]+s2[1]+s2[2])
end

function (up::undeformed_RHP_solitons)(x,t; flip_tol = 10., pole_circ = 0.001, flips=nothing, max_deriv_terms=20, verbose=true)
    g = size(up.bands,1)-1
    ntot = sum(up.nvec)
    nptot(j) = sum(up.nvec[1:j])

    res_conds = up.χvec.*exp.(-2x*up.κvec.+8t*up.κvec.^3)

    if flips === nothing#check whether or not to flip problem
        flips = abs.(res_conds).>flip_tol
    end
    flipped_poles = up.κvec[flips]

    A₂₁ = zeros(ComplexF64,ntot,ntot); rhs₂₁=zeros(ComplexF64,ntot,1)
    A₁₂ = zeros(ComplexF64,ntot,ntot); rhs₁₂=zeros(ComplexF64,ntot,1)
    
    for j = 1:g+1
        vvals = map(z->v(z,flipped_poles)^2, up.gridpts[j])
        if j > (g+1)/2
            jump_vals = -im*up.rvals[j].*exp.(-2x*up.gridpts[j]+8t*up.gridpts[j].^3).*vvals
            A₂ = -Diagonal(jump_vals)*up.Cmats[j,1]
            A₂₁[nptot(j-1)+1:nptot(j),:] = A₂
            rhs₂₁[nptot(j-1)+1:nptot(j)] = jump_vals
        else
            jump_vals = im*up.rvals[j].*exp.(2x*up.gridpts[j]-8t*up.gridpts[j].^3)./vvals
            A₂ = -Diagonal(jump_vals)*up.Cmats[j,2]
            A₁₂[nptot(j-1)+1:nptot(j),:] = A₂
            rhs₁₂[nptot(j-1)+1:nptot(j)] = jump_vals
        end
    end
        
    A = [up.A₁₁ A₂₁; A₁₂ up.A₂₂]
    rhs = [zeros(ntot) rhs₂₁; rhs₁₂ zeros(ntot)]
    mat_coeffs = A\rhs
    #println(cond(A))

    #functon to get out matrix solution
    function eval_sol(z;flag=0)
        S = zeros(ComplexF64,2,2)
        for j = 1:g+1
            for k = 1:2
                #sum Chebyshev series
                polys = up.Chebymat[j,k]
                S[:,k] += transpose(CauchyInterval(z,polys,up.nvec[j]-1;flag=flag)*mat_coeffs[nptot(j-1)+(k-1)*ntot+1:nptot(j)+(k-1)*ntot,:])
            end
        end
        return S+I(2)
    end
    poles = vcat(up.κvec, -up.κvec)
    Vr = eval_sol.(poles)

    #get derivative terms
    function get_der(κ;epsc=pole_circ, max_terms=max_deriv_terms)
        d11 = BasisExpansion(z->eval_sol(z)[1,1],Laurent(PeriodicMappedCircle(κ,epsc)),max_terms)
        d12 = BasisExpansion(z->eval_sol(z)[1,2],Laurent(PeriodicMappedCircle(κ,epsc)),max_terms)
        d21 = BasisExpansion(z->eval_sol(z)[2,1],Laurent(PeriodicMappedCircle(κ,epsc)),max_terms)
        d22 = BasisExpansion(z->eval_sol(z)[2,2],Laurent(PeriodicMappedCircle(κ,epsc)),max_terms)
        tol_met = norm([norm(d11.c[end-4:end]) norm(d12.c[end-4:end]); norm(d21.c[end-4:end]) norm(d22.c[end-4:end])])
        if  tol_met > 1e-11 && verbose
            @warn "Last terms in Δ Taylor expansion are not small. Tolerance met: "*string(tol_met)
        end
        f1p = [CauchyTransform()*(Derivative(1)*d11) CauchyTransform()*(Derivative(1)*d12);
        CauchyTransform()*(Derivative(1)*d21) CauchyTransform()*(Derivative(1)*d22)]
        return f1p(κ)
    end
    Vrp = get_der.(poles)
    #println(det(Vr2))
    #println(abs.(sum(Vr1,dims=1)-reverse(sum(Vr2,dims=1))))

    #Solve residue problem
    
    R = Vector{Matrix{ComplexF64}}(undef,length(poles))
    for (j,κ) in enumerate(up.κvec)
        k = j+length(up.κvec)
        if flips[j] == 1
            R[j] = Vr[j]*[0 1/(res_conds[j]*dv(κ,flipped_poles)^2); 0 0]*inv(Vr[j]-Vrp[j]*[0 1/(res_conds[j]*dv(κ,flipped_poles)^2); 0 0])#*inv(Vr[j])*inv(I(2)-Vrp[j]*[0 1/(res_conds[j]*dv(κ,flipped_poles)^2); 0 0]*inv(Vr[j]))
            R[k] = Vr[k]*[0 0; -1/(res_conds[j]*dv(κ,flipped_poles)^2) 0]*inv(Vr[k]-Vrp[k]*[0 0; -1/(res_conds[j]*dv(κ,flipped_poles)^2) 0])#*inv(Vr[k])*inv(I(2)-Vrp[k]*[0 0; -1/(res_conds[j]*dv(κ,flipped_poles)^2) 0]*inv(Vr[k]))
        else
            R[j] = Vr[j]*[0 0; res_conds[j]*v(κ,flipped_poles)^2 0]*inv(Vr[j]-Vrp[j]*[0 0; res_conds[j]*v(κ,flipped_poles)^2 0])#*inv(Vr[j])*inv(I(2)-Vrp[j]*[0 0; res_conds[j]*v(κ,flipped_poles)^2 0]*inv(Vr[j]))
            R[k] = Vr[k]*[0 -res_conds[j]*v(κ,flipped_poles)^2; 0 0]*inv(Vr[k]-Vrp[k]*[0 -res_conds[j]*v(κ,flipped_poles)^2; 0 0])#*inv(Vr[k])*inv(I(2)-Vrp[k]*[0 -res_conds[j]*v(κ,flipped_poles)^2; 0 0]*inv(Vr[k]))
        end
    end
    pole_rhp = RHP([],[],poles,R)
    rhsolver = RHSolver(pole_rhp);
    res_sol = rhsolver([1 1], 300)
    #some_val = moment(res_sol[1],0)*1im/(2pi)
    #println(some_val)

    #get behavior at infinity for residue solution
    res_a1 = [moment(res_sol[1],0) moment(res_sol[2],0)]*1im/(2pi)
    res_a2 = [moment(res_sol[1],1) moment(res_sol[2],1)]*1im/(2pi)
    #println(res_a1)
    #println(res_a2)
    
    #get behavior at infinity for gas solution
    cont_a1 = zeros(ComplexF64,2,2); cont_a2 = zeros(ComplexF64,2,2)
    for j = 1:g+1
        for k = 1:2
            #Chebyshev series
            cont_a1[:,k] += mat_coeffs[nptot(j-1)+(k-1)*ntot+1,:]*im/2π
            if up.Chebymat[j,k].kind == 3
                cont_a2[:,k] += mat_coeffs[nptot(j-1)+(k-1)*ntot+1,:]*(up.bands[j,1]+3up.bands[j,2])*im/8π
            elseif up.Chebymat[j,k].kind == 4
                cont_a2[:,k] += mat_coeffs[nptot(j-1)+(k-1)*ntot+1,:]*(3up.bands[j,1]+up.bands[j,2])*im/8π
            else
                cont_a2[:,k] += mat_coeffs[nptot(j-1)+(k-1)*ntot+1,:]*(up.bands[j,1]+up.bands[j,2])*im/4π
            end
            if up.Chebymat[j,k].kind != 1
                cont_a2[:,k] += mat_coeffs[nptot(j-1)+(k-1)*ntot+2,:]*(up.bands[j,2]-up.bands[j,1])*im/8π
            else
                cont_a2[:,k] += √2*mat_coeffs[nptot(j-1)+(k-1)*ntot+2,:]*(up.bands[j,2]-up.bands[j,1])*im/8π
            end
        end
    end
    
    #get behavior at infinity for product
    s1 = res_a1+[1 1]*cont_a1
    s2 = res_a2+[1 1]*cont_a2+res_a1*cont_a1
    #println(s1)
    #println(s2)

    return 2*(s1[1]*s1[2]+s2[1]+s2[2])
end

function (rp::rhp_solitons)(x,t; flip_tol = 10., pole_circ = 0.001, flips=nothing, max_deriv_terms=25, verbose=true)
    if x>4t*rp.dp.bands[end,2]^2
        #println("Solving undeformed problem")
        return rp.up(x,t; flip_tol = flip_tol, pole_circ = pole_circ, flips=flips, max_deriv_terms=max_deriv_terms, verbose=verbose)    
    else
        #println("Solving deformed problem")
        return rp.dp(x,t; flip_tol = flip_tol, pole_circ = pole_circ, flips=flips, max_deriv_terms=max_deriv_terms, verbose=verbose)
    end
end

end # module KdVSolitonGas