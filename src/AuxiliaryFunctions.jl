R = pts -> ( z ->  map(y -> sqrt(z - y |> Complex),pts) |> prod)
R⁺ = pts -> ( z ->  map(y -> sqrt(z+eps()im - y |> Complex),pts) |> prod)

function p(j,x,bands) 
    gaps = hcat(bands[1:end-1,2],bands[2:end,1])
    prod(x .- sum(gaps[1:end .!= j,:],dims=2)/2)
end

function g_coeffs(bands)
    gaps = hcat(bands[1:end-1,2],bands[2:end,1])
    g = size(bands,1)-1
    A = zeros(ComplexF64,g,g)
    bx = zeros(ComplexF64,g)
    bt = zeros(ComplexF64,g)
    
    #integrate over gaps
    for k = 1:g
        out_points = vcat(bands[:,1][1:end .!= k+1],bands[:,2][1:end .!= k])
        gd = JacobiMappedInterval(gaps[k,1],gaps[k,2],-0.5,-0.5)
        sp = OperatorApproximation.Jacobi(-0.5,-0.5,gd)
        for j = 1:g
            f = BasisExpansion( z -> p(j,z,bands)/R(out_points)(z), sp)    
            if length(f.c) >= 10000
                @warn "OperatorApproximation error from constant term in linear system for 𝔤-function"
            end
            A[j,k] = -im*π*f.c[1]
        end
    end
    #println("g matrix condition number: ",cond(A))
    
    #integrate over bands
    for k = 1:g+1
        out_points = bands[1:end .!= k,:]
        gd = JacobiMappedInterval(bands[k,1],bands[k,2],-0.5,-0.5)
        sp = OperatorApproximation.Jacobi(-0.5,-0.5,gd)
        for j = 1:g
            #x integral
            fx = BasisExpansion(z -> 2z*p(j,z,bands)/R⁺(out_points)(z), sp)
            if length(fx.c) >= 10000
                @warn "OperatorApproximation error from x integral in linear system for 𝔤-function"
            end
            bx[j] += im*π*fx.c[1]
            
            #t integral
            ft = BasisExpansion(z -> 8z^3*p(j,z,bands)/R⁺(out_points)(z), sp)
            if length(ft.c) >= 10000
                @warn "OperatorApproximation error from t integral in linear system for 𝔤-function"
            end
            bt[j] -= im*π*ft.c[1]
        end
    end
    Ωx = A\bx
    Ωt = A\bt
    Ω(x,t) = x*Ωx+t*Ωt
end

function cheby_gap(bands)
    g = size(bands,1)-1
    gaps = hcat(bands[1:end-1,2],bands[2:end,1])
    cheby_coeffs = Array{Array{ComplexF64,1}}(undef,g)
    for j = 1:g
        out_points = vcat(bands[:,1][1:end .!= j+1],bands[:,2][1:end .!= j])
        gd = JacobiMappedInterval(gaps[j,1],gaps[j,2],-0.5,-0.5)
        sp = OperatorApproximation.Jacobi(-0.5,-0.5,gd)
        f = BasisExpansion( z -> 1/R(out_points)(z), sp)
        if length(f.c) >=  10000
            @warn "OperatorApproximation error from integrating over gaps"
        end
        cheby_coeffs[j] = f.c
    end
    cheby_coeffs
end

function cheby_int_2z(bands)
    g = size(bands,1)-1
    cheby_coeffs = Array{Array{ComplexF64,1}}(undef,g+1)
    for j = 1:g+1
        out_points = bands[1:end .!= j,:]
        gd = JacobiMappedInterval(bands[j,1],bands[j,2],-0.5,-0.5)
        sp = OperatorApproximation.Jacobi(-0.5,-0.5,gd)
        f = BasisExpansion( z -> 2z/R⁺(out_points)(z), sp)
        if length(f.c) >= 10000
            @warn "OperatorApproximation error from x-dependent term of 𝔤-function"
        end
        cheby_coeffs[j] = f.c
    end
    cheby_coeffs
end

function cheby_int_8z3(bands)
    g = size(bands,1)-1
    cheby_coeffs = Array{Array{ComplexF64,1}}(undef,g+1)
    for j = 1:g+1
        out_points = bands[1:end .!= j,:]
        gd = JacobiMappedInterval(bands[j,1],bands[j,2],-0.5,-0.5)
        sp = OperatorApproximation.Jacobi(-0.5,-0.5,gd)
        f = BasisExpansion( z -> 8z^3/R⁺(out_points)(z), sp)
        if length(f.c) >= 10000
            @warn "OperatorApproximation error from t-dependent term of 𝔤-function"
        end
        cheby_coeffs[j] = f.c
    end
    cheby_coeffs
end

function cheby_int(bands)
    g = size(bands,1)-1
    cheby_coeffs = Array{Array{ComplexF64,1}}(undef,g+1)
    for j = 1:g+1
        out_points = bands[1:end .!= j,:]
        gd = JacobiMappedInterval(bands[j,1],bands[j,2],-0.5,-0.5)
        sp = OperatorApproximation.Jacobi(-0.5,-0.5,gd)
        f = BasisExpansion( z -> 1/R⁺(out_points)(z), sp)
        if length(f.c) >= 10000
            @warn "OperatorApproximation error from integrating over bands"
        end
        cheby_coeffs[j] = f.c
    end
    cheby_coeffs
end

function compute_g_pre(z,bands,int_vals_2z,int_vals_8z3,gap_vals)
    g = size(bands,1)-1
    #gaps = hcat(bands[1:end-1,2],bands[2:end,1])
    gvals = zeros(ComplexF64,g+1,3)
    
    #integrate over bands
    for j = 1:g+1
        #out_points = bands[1:end .!= j,:]
        coeffsx = int_vals_2z[j]
        ncx = length(coeffsx)
        #println(nc)
        ChebyT = buildCheby(bands[j,1],bands[j,2],1)
        CauchyTx = CauchyInterval(z,ChebyT,ncx-1)
        gvals[j,1] = -im*π*(CauchyTx*coeffsx)[1]
        
        coeffst = int_vals_8z3[j]
        nct = length(coeffst)
        CauchyTt = CauchyInterval(z,ChebyT,nct-1)
        gvals[j,2] += im*π*(CauchyTt*coeffst)[1]
    end
    
    #integrate over gaps
    for j = 1:g
        #out_points = vcat(bands[:,1][1:end .!= j+1],bands[:,2][1:end .!= j])
        coeffs = gap_vals[j]
        nc = length(coeffs)
        #println(nc)
        ChebyT = buildCheby(bands[j,2],bands[j+1,1],1)
        CauchyT = CauchyInterval(z,ChebyT,nc-1)
        gvals[j,3] = -im*π*(CauchyT*coeffs)[1]  
    end
    
    return R(bands)(z)*gvals
end

function compute_g_post(gvals,x,t,Ω)
    Ωvec = Ω(x,t)
    sum(x*gvals[:,1]+t*gvals[:,2]+[Ωvec; 0].*gvals[:,3])
end

q(j,x,bands) = prod(x .- sum(bands[1:end .!= j,:],dims=2)/2)

function h_coeffs_pre(bands)
    # solve for coefficients to remove jump on gap
        gaps = hcat(bands[1:end-1,2],bands[2:end,1])
        g = size(bands,1)-1
        A = zeros(ComplexF64,g+1,g+1)
        B = zeros(ComplexF64,g+1,g)
        #integrate over bands
        for k = 1:g+1
            out_points = bands[1:end .!= k,:]
            gd = JacobiMappedInterval(bands[k,1],bands[k,2],-0.5,-0.5)
            sp = OperatorApproximation.Jacobi(-0.5,-0.5,gd)
            for j = 1:g+1
                f = BasisExpansion(z -> q(j,z,bands)/R(out_points)(z), sp)
                if length(f.c) >= 10000
                    @warn "OperatorApproximation error from integrating over bands in linear system for 𝔥-function"
                end
                A[j,k] = -im*π*f.c[1]
            end
        end
        #println("h matrix condition number:",cond(A))   
        #integrate over gaps
        for k = 1:g
            bk = zeros(ComplexF64,g+1,1)
            out_points = vcat(bands[:,1][1:end .!= k+1],bands[:,2][1:end .!= k])
            gd = JacobiMappedInterval(gaps[k,1],gaps[k,2],-0.5,-0.5)
            sp = OperatorApproximation.Jacobi(-0.5,-0.5,gd)
            for j = 1:g+1
                f = BasisExpansion(z -> q(j,z,bands)/R(out_points)(z), sp)   
                if length(f.c) >= 10000
                    @warn "OperatorApproximation error from integrating over gaps in linear system for 𝔥-function"
                end
                bk[j] = -im*π*f.c[1]
            end
            B[:,k] = bk
        end
(A,B)
end

function h_coeffs_post(A,B,Ωs)
    b = -B*log.(exp.(-Ωs))   
    Avec = A\b
    #=if maximum(imag(Avec))>1e-12
        #println("Warning: computed h coefficients nonreal. Maximum imaginary part printed")
        #println(maximum(imag(Avec)))
    end=#
    #Avec = real(Avec)
    Avec
end

function compute_h_pre(z,bands,int_vals,gap_vals)
    g = size(bands,1)-1
    h_vals = zeros(ComplexF64,2,g+1)
    for i = 1:g+1
        #out_points = bands[1:end .!= i,:]
        coeffs = int_vals[i,1]
        nc = length(coeffs)
        #println(nc)
        ChebyT = buildCheby(bands[i,1],bands[i,2],1)
        CauchyT = CauchyInterval(z,ChebyT,nc-1)
        #global svde = CauchyT
        h_vals[1,i] = (-π*im*CauchyT*coeffs)[1]
    end
    for i=1:g
        #out_points = vcat(bands[:,1][1:end .!= i+1],bands[:,2][1:end .!= i])
        coeffs = gap_vals[i,1]
        nc = length(coeffs)
        #println(nc)
        ChebyT = buildCheby(bands[i,2],bands[i+1,1],1)
        CauchyT = CauchyInterval(z,ChebyT,nc-1)
        h_vals[2,i] = (-π*im*CauchyT*coeffs)[1]        
    end
    return R(bands)(z)*h_vals
end

function compute_h_post(h_vals,Avec,Ωs)
    sum(diag(h_vals*[Avec [log.(exp.(-Ωs)); 0]]))  
end

function v(z,A) #partial Blaschke product evaluated at a positive flipped residue
    p = 1.
    for a in A
        p *= (z - a)/(z + a)
    end
    p
end

function dv(c,A) #derivative evaluated at a positive flipped residue
    p = 1.
    for a in A
        if a != c
            p *= (c - a)/(c + a)
        end
    end
    p/2c
end

#outputted functions to compute g/h
struct gfunction
    Ω::Function
    bands::Array{Float64,2}
    int_vals_2z
    int_vals_8z3
    gap_vals
end

function get_g(intervals)
    bands = [-reverse(intervals); intervals]
    Ω = g_coeffs(bands)
    int_vals_2z = cheby_int_2z(bands)
    int_vals_8z3 = cheby_int_8z3(bands)
    gap_vals = cheby_gap(bands)
    return gfunction(Ω, bands, int_vals_2z, int_vals_8z3, gap_vals)
end

function (g::gfunction)(z,x,t)
    gvals = compute_g_pre(z,g.bands,g.int_vals_2z,g.int_vals_8z3,g.gap_vals)
    return compute_g_post(gvals,x,t,g.Ω)
end

struct hfunction
    Ω::Function
    bands::Array{Float64,2}
    int_vals
    gap_vals
end

function get_h(intervals)
    bands = [-reverse(intervals); intervals]
    Ω = g_coeffs(bands)
    int_vals = cheby_int(bands)
    gap_vals = cheby_gap(bands)
    return hfunction(Ω, bands, int_vals, gap_vals)
end

function (h::hfunction)(z,x,t)
    (A,B) = h_coeffs_pre(h.bands)
    Ωvec = h.Ω(x,t)
    Avec = h_coeffs_post(A,B,Ωvec)

    h_vals = compute_h_pre(z,h.bands,h.int_vals,h.gap_vals)
    return compute_h_post(h_vals,Avec,Ωvec)
end