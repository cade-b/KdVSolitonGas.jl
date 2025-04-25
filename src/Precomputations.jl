function precompute(intervals::Array{Float64,2},h::Function,types::Vector{Int}; nmat=nothing, circ_size=1.25)
    bands = [-reverse(intervals); intervals]
    typevec = [reverse(types); types]
    if nmat === nothing
        nmat = [120*ones(size(bands,1)) 20*ones(size(bands,1))] .|> Int
    else 
        nmat = [reverse(nmat;dims=1); nmat]
    end
    
    int_vals = cheby_int(bands)
    int_vals_2z = cheby_int_2z(bands)
    int_vals_8z3 = cheby_int_8z3(bands)
    gap_vals = cheby_gap(bands)
    g = size(bands,1)-1

    function r(j)
        if j > (g+1)/2 #Σ₊
            if typevec[j] == 1 #T
                return z -> h(j-(g+1)/2)(z)/(√(bands[j,2]-z |> Complex)*√(z-bands[j,1] |> Complex))
            elseif typevec[j] == 2 #U
                return z -> h(j-(g+1)/2)(z)*(√(bands[j,2]-z |> Complex)*√(z-bands[j,1] |> Complex))
            elseif typevec[j] == 3 #V
                return z -> h(j-(g+1)/2)(z)*(√(z-bands[j,1] |> Complex)/√(bands[j,2]-z |> Complex))
            else #W
                return z -> h(j-(g+1)/2)(z)*(√(bands[j,2]-z |> Complex)/√(z-bands[j,1] |> Complex))
            end
        else
            return z -> r(g+2-j)(-z)
        end
    end

    cc(j) = (bands[j,1]+bands[j,2])/2
    rr(j) = circ_size*(bands[j,2]-bands[j,1])/2  

    gridmat = Array{Array{ComplexF64,1},2}(undef,g+1,2)
    for j = 1:g+1
        gridmat[j,1] = rr(j)*zgrid(nmat[j,1]).+cc(j).+eps()im
        gridmat[j,2] = M(bands[j,1],bands[j,2]).(Tgrid(nmat[j,2])) .|> Complex
    end

    basisperm1 = [2, 1, 3, 4]
    basisperm2 = [1, 2, 4, 3];
    Chebymat = Array{ChebyParams}(undef,g+1,2)
    for j = 1:g+1
        if j > (g+1)/2 #Σ₊
            Chebymat[j,1] = buildCheby(bands[j,1],bands[j,2],typevec[j]) #left column
            Chebymat[j,2] = buildCheby(bands[j,1],bands[j,2],typevec[j]+2*(typevec[j]%2)-1) #right column
        else #Σ₋
            Chebymat[j,1] = buildCheby(bands[j,1],bands[j,2],basisperm1[typevec[j]]) 
            Chebymat[j,2] = buildCheby(bands[j,1],bands[j,2],basisperm2[typevec[j]]) 
        end
    end

    ntot = sum(sum(nmat))
    nptot(j) = sum(sum(nmat[1:j,:]))
    
    #build top left corner
    #build the first row of the block system
    A₁₁ = zeros(ComplexF64,ntot,ntot); #rhs₁₁=zeros(ComplexF64,ntot,1)
    for j = 1:g+1
        CM₊ = zeros(ComplexF64,nmat[j,1],ntot)
        CM₋ = zeros(ComplexF64,nmat[j,1],ntot)
        CM₁₊= CauchyMat(gridmat[j,1],nmat[j,1],cc(j),rr(j);flag=1)
        CM₁₋= CauchyMat(gridmat[j,1],nmat[j,1],cc(j),rr(j);flag=-1)
        CM₂ = CauchyInterval(gridmat[j,1],Chebymat[j,1],nmat[j,2]-1)
        CM₊[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁₊
        CM₋[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁₋
        CM₊[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂
        CM₋[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,1],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,1],Chebymat[k,1],nmat[k,2]-1)
            CM₊[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]= CM₁
            CM₋[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]= CM₁
            CM₊[:,nptot(k)-nmat[k,2]+1:nptot(k)]= CM₂
            CM₋[:,nptot(k)-nmat[k,2]+1:nptot(k)]= CM₂
        end
        A₁ = CM₊-CM₋
        #build the second row of the block system
        CM₊ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₁ = CauchyMat(gridmat[j,2],nmat[j,1],cc(j),rr(j))
        CM₂₊= CauchyInterval(gridmat[j,2],Chebymat[j,1],nmat[j,2]-1;flag=1)
        CM₊[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁
        CM₊[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂₊
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,2],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,1],nmat[k,2]-1)
            CM₊[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]=CM₁
            CM₊[:,nptot(k)-nmat[k,2]+1:nptot(k)]=CM₂
        end
        #println(maximum(abs.(CM₊-CM₋)))
        A₂ = CM₊
        #global testguy=A₂
        #combine and build RHS
        A₁₁[nptot(j-1)+1:nptot(j),:] = [A₁; A₂]
        #rhs₁₁[nptot(j-1)+1:nptot(j)] = [zeros(nmat[j,1]); -ones(nmat[j,2])]
    end
    
    Cmats = Array{Matrix{ComplexF64}}(undef,g+1,4)
    for j = 1:g+1
        #build portions of the first row of the block system
        CM₋ = zeros(ComplexF64,nmat[j,1],ntot)
        CM₁₋= CauchyMat(gridmat[j,1],nmat[j,1],cc(j),rr(j);flag=-1)
        CM₂ = CauchyInterval(gridmat[j,1],Chebymat[j,2],nmat[j,2]-1)
        CM₋[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁₋
        CM₋[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,1],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,1],Chebymat[k,2],nmat[k,2]-1)
            CM₋[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]= CM₁
            CM₋[:,nptot(k)-nmat[k,2]+1:nptot(k)]= CM₂
        end
        Cmats[j,1] = CM₋
        #build the second row of the block system
        CM₋ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₁ = CauchyMat(gridmat[j,2],nmat[j,1],cc(j),rr(j))
        CM₂₋= CauchyInterval(gridmat[j,2],Chebymat[j,2],nmat[j,2]-1;flag=-1)
        CM₋[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁
        CM₋[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂₋
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,2],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,2],nmat[k,2]-1)
            CM₋[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]=CM₁
            CM₋[:,nptot(k)-nmat[k,2]+1:nptot(k)]=CM₂
        end
        Cmats[j,2] = CM₋
    end

    #build bottom left corner
    #build the first row of the block system
    for j = 1:g+1
        CM₋ = zeros(ComplexF64,nmat[j,1],ntot)
        CM₁₋= CauchyMat(gridmat[j,1],nmat[j,1],cc(j),rr(j);flag=-1)
        CM₂ = CauchyInterval(gridmat[j,1],Chebymat[j,1],nmat[j,2]-1)
        CM₋[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁₋
        CM₋[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,1],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,1],Chebymat[k,1],nmat[k,2]-1)
            CM₋[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]= CM₁
            CM₋[:,nptot(k)-nmat[k,2]+1:nptot(k)]= CM₂
        end
        Cmats[j,3] = CM₋

        #build the second row of the block system
        CM₋ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₁ = CauchyMat(gridmat[j,2],nmat[j,1],cc(j),rr(j))
        CM₂₋ = CauchyInterval(gridmat[j,2],Chebymat[j,1],nmat[j,2]-1;flag=-1)
        CM₋[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁
        CM₋[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂₋
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,2],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,1],nmat[k,2]-1)
            CM₋[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]= CM₁
            CM₋[:,nptot(k)-nmat[k,2]+1:nptot(k)]= CM₂
        end
        Cmats[j,4] = CM₋
    end
    
    #build bottom right corner
    #build the first row of the block system
    A₂₂ = zeros(ComplexF64,ntot,ntot); #rhs₂₂=zeros(ComplexF64,ntot,1)
    for j = 1:g+1
        CM₊ = zeros(ComplexF64,nmat[j,1],ntot)
        CM₋ = zeros(ComplexF64,nmat[j,1],ntot)
        CM₁₊ = CauchyMat(gridmat[j,1],nmat[j,1],cc(j),rr(j);flag=1)
        CM₁₋ = CauchyMat(gridmat[j,1],nmat[j,1],cc(j),rr(j);flag=-1)
        CM₂ = CauchyInterval(gridmat[j,1],Chebymat[j,2],nmat[j,2]-1)
        CM₊[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁₊
        CM₋[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁₋
        CM₊[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂
        CM₋[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,1],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,1],Chebymat[k,2],nmat[k,2]-1)
            CM₊[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]=CM₁
            CM₋[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]=CM₁
            CM₊[:,nptot(k)-nmat[k,2]+1:nptot(k)]=CM₂
            CM₋[:,nptot(k)-nmat[k,2]+1:nptot(k)]=CM₂
        end
        gn₁ = ones(nmat[j,1])
        A₁ = CM₊-Diagonal(gn₁)*CM₋
        #build the second row of the block system
        CM₊ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₁ = CauchyMat(gridmat[j,2],nmat[j,1],cc(j),rr(j))
        CM₂₊ = CauchyInterval(gridmat[j,2],Chebymat[j,2],nmat[j,2]-1;flag=1)
        CM₊[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]=CM₁
        CM₊[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂₊
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,2],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,2],nmat[k,2]-1)
            CM₊[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]= CM₁
            CM₊[:,nptot(k)-nmat[k,2]+1:nptot(k)]= CM₂
        end
        A₂ = CM₊
        #combine and build RHS
        A₂₂[nptot(j-1)+1:nptot(j),:] = [A₁; A₂]
        #rhs₂₂[nptot(j-1)+1:nptot(j)] = [zeros(nmat[j,1]); -ones(nmat[j,2])]
    end
    
    Ω = g_coeffs(bands)
    #store terms for efficient computation of 𝔤(z) on circles
    get_g_pre(z) = compute_g_pre(z,bands,int_vals_2z,int_vals_8z3,gap_vals)
    g_vals_pre = Array{Vector{Matrix{ComplexF64}}}(undef,g+1)
    for j = 1:g+1
        g_vals_pre[j] = map(λ -> get_g_pre(λ), gridmat[j,1])
    end

    #store values of reflection coefficients on intervals and circles
    m_r(j) = z -> imag(z)>0 ? 1/r(j)(z) : -1/r(j)(z)
    r_vals = Array{Vector{ComplexF64}}(undef,g+1,3)
    for j = 1:g+1
        r_vals[j,1] = map(z->r(j)(z), gridmat[j,2])
        r_vals[j,2] = map(z->1/r(j)(z), gridmat[j,2])
        r_vals[j,3] = map(z->m_r(j)(z), gridmat[j,1])
    end

    #store terms for efficient computation of 𝔥(z)
    (Ah_pre,Bh_pre) = h_coeffs_pre(bands)

    get_h_pre(z) = compute_h_pre(z, bands, int_vals, gap_vals)
    h_vals_pre = Array{Vector{Matrix{ComplexF64}}}(undef,g+1)
    for j = 1:g+1
        h_vals_pre[j] = map(z->get_h_pre(z), gridmat[j,1])
    end
    
    #get everything for undeformed problem
    Chebymatp = Array{ChebyParams}(undef,g+1,2)
    for j = 1:g+1
        if j > (g+1)/2 #Σ₊
            Chebymatp[j,1] = buildCheby(bands[j,1],bands[j,2],typevec[j]) #left column
            Chebymatp[j,2] = buildCheby(bands[j,1],bands[j,2],2) #right column
        else #Σ₋
            Chebymatp[j,1] = buildCheby(bands[j,1],bands[j,2],2) 
            Chebymatp[j,2] = buildCheby(bands[j,1],bands[j,2],basisperm2[typevec[j]]) 
        end
    end
    
    ntot = sum(nmat[:,2])
    nptot1(j) = sum(nmat[1:j,2])
    A₁₁p = zeros(ComplexF64,ntot,ntot); #rhs₁₁=zeros(ComplexF64,ntot,1)
    for j = 1:g+1
        CM₊ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₋ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₂₊= CauchyInterval(gridmat[j,2],Chebymat[j,1],nmat[j,2]-1;flag=1)
        CM₂₋= CauchyInterval(gridmat[j,2],Chebymat[j,1],nmat[j,2]-1;flag=-1)
        CM₊[:,nptot1(j-1)+1:nptot1(j)]= CM₂₊
        CM₋[:,nptot1(j-1)+1:nptot1(j)]= CM₂₋
        for k=(1:g+1)[1:end .!= j,:]
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,1],nmat[k,2]-1)
            CM₊[:,nptot1(k-1)+1:nptot1(k)] = CM₂
            CM₋[:,nptot1(k-1)+1:nptot1(k)] = CM₂
        end
        A₂ = CM₊-CM₋
        #combine and build block
        A₁₁p[nptot1(j-1)+1:nptot1(j),:] = A₂
    end
    
    Cmatsp = Array{Matrix{ComplexF64}}(undef,g+1,2)
    for j = 1:g+1
        CM₋ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₂₋= CauchyInterval(gridmat[j,2],Chebymat[j,2],nmat[j,2]-1;flag=-1)
        CM₋[:,nptot1(j-1)+1:nptot1(j)]= CM₂₋
        for k=(1:g+1)[1:end .!= j,:]
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,2],nmat[k,2]-1)
            CM₋[:,nptot1(k-1)+1:nptot1(k)]=CM₂
        Cmatsp[j,1] = CM₋
        end
    end
    
    for j = 1:g+1
        CM₋ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₂₋ = CauchyInterval(gridmat[j,2],Chebymat[j,1],nmat[j,2]-1;flag=-1)
        CM₋[:,nptot1(j-1)+1:nptot1(j)]= CM₂₋
        for k=(1:g+1)[1:end .!= j,:]
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,1],nmat[k,2]-1)
            CM₋[:,nptot1(k-1)+1:nptot1(k)]= CM₂
        end
        Cmatsp[j,2] = CM₋
    end

    A₂₂p = zeros(ComplexF64,ntot,ntot);
    for j = 1:g+1
        CM₊ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₋ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₂₊ = CauchyInterval(gridmat[j,2],Chebymat[j,2],nmat[j,2]-1;flag=1)
        CM₂₋ = CauchyInterval(gridmat[j,2],Chebymat[j,2],nmat[j,2]-1;flag=-1)
        CM₊[:,nptot1(j-1)+1:nptot1(j)]= CM₂₊
        CM₋[:,nptot1(j-1)+1:nptot1(j)]= CM₂₋
        for k=(1:g+1)[1:end .!= j,:]
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,2],nmat[k,2]-1)
            CM₊[:,nptot1(k-1)+1:nptot1(k)] = CM₂
            CM₋[:,nptot1(k-1)+1:nptot1(k)] = CM₂
        end
        A₂ = CM₊-CM₋
        #combine and build block
        A₂₂p[nptot1(j-1)+1:nptot1(j),:] = A₂
    end

    up = undeformed_RHP(bands,A₁₁p,A₂₂p,Cmatsp,Chebymatp,nmat[:,2],gridmat[:,2],r_vals[:,1])
    dp = deformed_RHP(bands,A₁₁,A₂₂,Cmats,circ_size,Ω,Chebymat,nmat,gridmat[:,1],g_vals_pre,h_vals_pre,Ah_pre,Bh_pre,r_vals)
    rhp(dp,up)
end

### Add distinguished soliton ####################################################################################################################
function precompute(intervals::Array{Float64,2},κ₀::Float64,χ::Float64,h::Function,types::Vector{Int}; nmat=nothing, circ_size=1.25)
    bands = [-reverse(intervals); intervals]
    typevec = [reverse(types); types]
    if nmat === nothing
        nmat = [120*ones(size(bands,1)) 20*ones(size(bands,1))] .|> Int
    else 
        nmat = [reverse(nmat;dims=1); nmat]
    end
    
    int_vals = cheby_int(bands)
    int_vals_2z = cheby_int_2z(bands)
    int_vals_8z3 = cheby_int_8z3(bands)
    gap_vals = cheby_gap(bands)
    g = size(bands,1)-1

    function r(j)
        if j > (g+1)/2 #Σ₊
            if typevec[j] == 1 #T
                return z -> h(j-(g+1)/2)(z)/(√(bands[j,2]-z |> Complex)*√(z-bands[j,1] |> Complex))
            elseif typevec[j] == 2 #U
                return z -> h(j-(g+1)/2)(z)*(√(bands[j,2]-z |> Complex)*√(z-bands[j,1] |> Complex))
            elseif typevec[j] == 3 #V
                return z -> h(j-(g+1)/2)(z)*(√(z-bands[j,1] |> Complex)/√(bands[j,2]-z |> Complex))
            else #W
                return z -> h(j-(g+1)/2)(z)*(√(bands[j,2]-z |> Complex)/√(z-bands[j,1] |> Complex))
            end
        else
            return z -> r(g+2-j)(-z)
        end
    end

    cc(j) = (bands[j,1]+bands[j,2])/2
    rr(j) = circ_size*(bands[j,2]-bands[j,1])/2

    gridmat = Array{Array{ComplexF64,1},2}(undef,g+1,2)
    for j = 1:g+1
        gridmat[j,1] = rr(j)*zgrid(nmat[j,1]).+cc(j).+eps()im
        gridmat[j,2] = M(bands[j,1],bands[j,2]).(Tgrid(nmat[j,2])) .|> Complex
    end

    basisperm1 = [2, 1, 3, 4]
    basisperm2 = [1, 2, 4, 3];
    Chebymat = Array{ChebyParams}(undef,g+1,2)
    for j = 1:g+1
        if j > (g+1)/2 #Σ₊
            Chebymat[j,1] = buildCheby(bands[j,1],bands[j,2],typevec[j]) #left column
            Chebymat[j,2] = buildCheby(bands[j,1],bands[j,2],typevec[j]+2*(typevec[j]%2)-1) #right column
        else #Σ₋
            Chebymat[j,1] = buildCheby(bands[j,1],bands[j,2],basisperm1[typevec[j]]) 
            Chebymat[j,2] = buildCheby(bands[j,1],bands[j,2],basisperm2[typevec[j]]) 
        end
    end

    ntot = sum(sum(nmat))
    nptot(j) = sum(sum(nmat[1:j,:]))
    
    #build top left corner
    #build the first row of the block system
    A₁₁ = zeros(ComplexF64,ntot,ntot); rhs₁₁=zeros(ComplexF64,ntot,1)
    for j = 1:g+1
        CM₊ = zeros(ComplexF64,nmat[j,1],ntot)
        CM₋ = zeros(ComplexF64,nmat[j,1],ntot)
        CM₁₊= CauchyMat(gridmat[j,1],nmat[j,1],cc(j),rr(j);flag=1)
        CM₁₋= CauchyMat(gridmat[j,1],nmat[j,1],cc(j),rr(j);flag=-1)
        CM₂ = CauchyInterval(gridmat[j,1],Chebymat[j,1],nmat[j,2]-1)
        CM₊[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁₊
        CM₋[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁₋
        CM₊[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂
        CM₋[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,1],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,1],Chebymat[k,1],nmat[k,2]-1)
            CM₊[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]= CM₁
            CM₋[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]= CM₁
            CM₊[:,nptot(k)-nmat[k,2]+1:nptot(k)]= CM₂
            CM₋[:,nptot(k)-nmat[k,2]+1:nptot(k)]= CM₂
        end
        A₁ = CM₊-CM₋
        #build the second row of the block system
        CM₊ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₁ = CauchyMat(gridmat[j,2],nmat[j,1],cc(j),rr(j))
        CM₂₊= CauchyInterval(gridmat[j,2],Chebymat[j,1],nmat[j,2]-1;flag=1)
        CM₊[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁
        CM₊[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂₊
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,2],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,1],nmat[k,2]-1)
            CM₊[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]=CM₁
            CM₊[:,nptot(k)-nmat[k,2]+1:nptot(k)]=CM₂
        end
        #println(maximum(abs.(CM₊-CM₋)))
        A₂ = CM₊
        #global testguy=A₂
        #combine and build RHS
        A₁₁[nptot(j-1)+1:nptot(j),:] = [A₁; A₂]
        rhs₁₁[nptot(j-1)+1:nptot(j)] = [zeros(nmat[j,1]); -ones(nmat[j,2])]
    end
    
    Cmats = Array{Matrix{ComplexF64}}(undef,g+1,4)
    for j = 1:g+1
        #build portions of the first row of the block system
        CM₋ = zeros(ComplexF64,nmat[j,1],ntot)
        CM₁₋= CauchyMat(gridmat[j,1],nmat[j,1],cc(j),rr(j);flag=-1)
        CM₂ = CauchyInterval(gridmat[j,1],Chebymat[j,2],nmat[j,2]-1)
        CM₋[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁₋
        CM₋[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,1],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,1],Chebymat[k,2],nmat[k,2]-1)
            CM₋[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]= CM₁
            CM₋[:,nptot(k)-nmat[k,2]+1:nptot(k)]= CM₂
        end
        Cmats[j,1] = CM₋
        #build the second row of the block system
        CM₋ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₁ = CauchyMat(gridmat[j,2],nmat[j,1],cc(j),rr(j))
        CM₂₋= CauchyInterval(gridmat[j,2],Chebymat[j,2],nmat[j,2]-1;flag=-1)
        CM₋[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁
        CM₋[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂₋
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,2],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,2],nmat[k,2]-1)
            CM₋[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]=CM₁
            CM₋[:,nptot(k)-nmat[k,2]+1:nptot(k)]=CM₂
        end
        Cmats[j,2] = CM₋
    end

    #build bottom left corner
    #build the first row of the block system
    for j = 1:g+1
        CM₋ = zeros(ComplexF64,nmat[j,1],ntot)
        CM₁₋= CauchyMat(gridmat[j,1],nmat[j,1],cc(j),rr(j);flag=-1)
        CM₂ = CauchyInterval(gridmat[j,1],Chebymat[j,1],nmat[j,2]-1)
        CM₋[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁₋
        CM₋[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,1],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,1],Chebymat[k,1],nmat[k,2]-1)
            CM₋[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]= CM₁
            CM₋[:,nptot(k)-nmat[k,2]+1:nptot(k)]= CM₂
        end
        Cmats[j,3] = CM₋

        #build the second row of the block system
        CM₋ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₁ = CauchyMat(gridmat[j,2],nmat[j,1],cc(j),rr(j))
        CM₂₋ = CauchyInterval(gridmat[j,2],Chebymat[j,1],nmat[j,2]-1;flag=-1)
        CM₋[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁
        CM₋[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂₋
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,2],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,1],nmat[k,2]-1)
            CM₋[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]= CM₁
            CM₋[:,nptot(k)-nmat[k,2]+1:nptot(k)]= CM₂
        end
        Cmats[j,4] = CM₋
    end
    
    #build bottom right corner
    #build the first row of the block system
    A₂₂ = zeros(ComplexF64,ntot,ntot); rhs₂₂=zeros(ComplexF64,ntot,1)
    for j = 1:g+1
        CM₊ = zeros(ComplexF64,nmat[j,1],ntot)
        CM₋ = zeros(ComplexF64,nmat[j,1],ntot)
        CM₁₊ = CauchyMat(gridmat[j,1],nmat[j,1],cc(j),rr(j);flag=1)
        CM₁₋ = CauchyMat(gridmat[j,1],nmat[j,1],cc(j),rr(j);flag=-1)
        CM₂ = CauchyInterval(gridmat[j,1],Chebymat[j,2],nmat[j,2]-1)
        CM₊[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁₊
        CM₋[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁₋
        CM₊[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂
        CM₋[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,1],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,1],Chebymat[k,2],nmat[k,2]-1)
            CM₊[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]=CM₁
            CM₋[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]=CM₁
            CM₊[:,nptot(k)-nmat[k,2]+1:nptot(k)]=CM₂
            CM₋[:,nptot(k)-nmat[k,2]+1:nptot(k)]=CM₂
        end
        gn₁ = ones(nmat[j,1])
        A₁ = CM₊-Diagonal(gn₁)*CM₋
        #build the second row of the block system
        CM₊ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₁ = CauchyMat(gridmat[j,2],nmat[j,1],cc(j),rr(j))
        CM₂₊ = CauchyInterval(gridmat[j,2],Chebymat[j,2],nmat[j,2]-1;flag=1)
        CM₊[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]=CM₁
        CM₊[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂₊
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,2],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,2],nmat[k,2]-1)
            CM₊[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]= CM₁
            CM₊[:,nptot(k)-nmat[k,2]+1:nptot(k)]= CM₂
        end
        A₂ = CM₊
        #combine and build RHS
        A₂₂[nptot(j-1)+1:nptot(j),:] = [A₁; A₂]
        rhs₂₂[nptot(j-1)+1:nptot(j)] = [zeros(nmat[j,1]); -ones(nmat[j,2])]
    end
    
    Ω = g_coeffs(bands)
    #store terms for efficient computation of 𝔤(z) on circles and at soliton
    get_g_pre(z) = compute_g_pre(z,bands,int_vals_2z,int_vals_8z3,gap_vals)
    g_vals_pre = Array{Vector{Matrix{ComplexF64}}}(undef,g+2)
    for j = 1:g+1
        g_vals_pre[j] = map(λ -> get_g_pre(λ), gridmat[j,1])
    end
    g_vals_pre[g+2] = map(λ -> get_g_pre(λ), [κ₀])

    #store values of reflection coefficients and soliton transformation on intervals and circles
    m_r(j) = z -> imag(z)>0 ? 1/r(j)(z) : -1/r(j)(z)
    v(λ) = (λ-κ₀)/(λ+κ₀)
    r_vals = Array{Vector{ComplexF64}}(undef,g+1,3)
    v_vals = Array{Vector{ComplexF64}}(undef,g+1,2)
    for j = 1:g+1
        r_vals[j,1] = map(z->r(j)(z), gridmat[j,2])
        r_vals[j,2] = map(z->1/r(j)(z), gridmat[j,2])
        r_vals[j,3] = map(z->m_r(j)(z), gridmat[j,1])

        v_vals[j,1] = map(λ->v(λ)^2, gridmat[j,2]) #store values on intervals
        v_vals[j,2] = map(λ->v(λ)^2, gridmat[j,1]) #store values on circles
    end

    #store terms for efficient computation of 𝔥(z)
    (Ah_pre,Bh_pre) = h_coeffs_pre(bands)

    get_h_pre(z) = compute_h_pre(z, bands, int_vals, gap_vals)
    h_vals_pre = Array{Vector{Matrix{ComplexF64}}}(undef,g+2)
    for j = 1:g+1
        h_vals_pre[j] = map(z->get_h_pre(z), gridmat[j,1])
    end
    h_vals_pre[g+2] = map(z->get_h_pre(z), [κ₀])
    
    #get everything for undeformed problem
    Chebymatp = Array{ChebyParams}(undef,g+1,2)
    for j = 1:g+1
        if j > (g+1)/2 #Σ₊
            Chebymatp[j,1] = buildCheby(bands[j,1],bands[j,2],typevec[j]) #left column
            Chebymatp[j,2] = buildCheby(bands[j,1],bands[j,2],2) #right column
        else #Σ₋
            Chebymatp[j,1] = buildCheby(bands[j,1],bands[j,2],2) 
            Chebymatp[j,2] = buildCheby(bands[j,1],bands[j,2],basisperm2[typevec[j]]) 
        end
    end
    
    ntot = sum(nmat[:,2])
    nptot1(j) = sum(nmat[1:j,2])
    A₁₁p = zeros(ComplexF64,ntot,ntot); #rhs₁₁=zeros(ComplexF64,ntot,1)
    for j = 1:g+1
        CM₊ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₋ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₂₊= CauchyInterval(gridmat[j,2],Chebymat[j,1],nmat[j,2]-1;flag=1)
        CM₂₋= CauchyInterval(gridmat[j,2],Chebymat[j,1],nmat[j,2]-1;flag=-1)
        CM₊[:,nptot1(j-1)+1:nptot1(j)]= CM₂₊
        CM₋[:,nptot1(j-1)+1:nptot1(j)]= CM₂₋
        for k=(1:g+1)[1:end .!= j,:]
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,1],nmat[k,2]-1)
            CM₊[:,nptot1(k-1)+1:nptot1(k)] = CM₂
            CM₋[:,nptot1(k-1)+1:nptot1(k)] = CM₂
        end
        A₂ = CM₊-CM₋
        #combine and build RHS
        A₁₁p[nptot1(j-1)+1:nptot1(j),:] = A₂
    end
    
    Cmatsp = Array{Matrix{ComplexF64}}(undef,g+1,2)
    for j = 1:g+1
        CM₋ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₂₋= CauchyInterval(gridmat[j,2],Chebymat[j,2],nmat[j,2]-1;flag=-1)
        CM₋[:,nptot1(j-1)+1:nptot1(j)]= CM₂₋
        for k=(1:g+1)[1:end .!= j,:]
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,2],nmat[k,2]-1)
            CM₋[:,nptot1(k-1)+1:nptot1(k)]=CM₂
        Cmatsp[j,1] = CM₋
        end
    end
    
    for j = 1:g+1
        CM₋ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₂₋ = CauchyInterval(gridmat[j,2],Chebymat[j,1],nmat[j,2]-1;flag=-1)
        CM₋[:,nptot1(j-1)+1:nptot1(j)]= CM₂₋
        for k=(1:g+1)[1:end .!= j,:]
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,1],nmat[k,2]-1)
            CM₋[:,nptot1(k-1)+1:nptot1(k)]= CM₂
        end
        Cmatsp[j,2] = CM₋
    end

    A₂₂p = zeros(ComplexF64,ntot,ntot);
    for j = 1:g+1
        CM₊ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₋ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₂₊ = CauchyInterval(gridmat[j,2],Chebymat[j,2],nmat[j,2]-1;flag=1)
        CM₂₋ = CauchyInterval(gridmat[j,2],Chebymat[j,2],nmat[j,2]-1;flag=-1)
        CM₊[:,nptot1(j-1)+1:nptot1(j)]= CM₂₊
        CM₋[:,nptot1(j-1)+1:nptot1(j)]= CM₂₋
        for k=(1:g+1)[1:end .!= j,:]
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,2],nmat[k,2]-1)
            CM₊[:,nptot1(k-1)+1:nptot1(k)] = CM₂
            CM₋[:,nptot1(k-1)+1:nptot1(k)] = CM₂
        end
        A₂ = CM₊-CM₋
        #combine and build RHS
        A₂₂p[nptot1(j-1)+1:nptot1(j),:] = A₂
    end

    up = undeformed_RHP_soliton(bands,A₁₁p,A₂₂p,Cmatsp,Chebymatp,nmat[:,2],gridmat[:,2],r_vals[:,1],v_vals[:,1],κ₀,χ)
    dp = deformed_RHP_soliton(bands,A₁₁,A₂₂,rhs₁₁,rhs₂₂,Cmats,circ_size,Ω,Chebymat,nmat,gridmat[:,1],g_vals_pre,h_vals_pre,Ah_pre,Bh_pre,r_vals,v_vals,κ₀,χ)
    rhp_soliton(dp,up)
end

### add several solitons #############################################################################################
function precompute(intervals::Array{Float64,2},κvec::Vector{Float64},χvec::Vector{Float64},h::Function,types::Vector{Int}; nmat=nothing, circ_size = 1.25)
    bands = [-reverse(intervals); intervals]
    typevec = [reverse(types); types]
    if nmat === nothing
        nmat = [120*ones(size(bands,1)) 20*ones(size(bands,1))] .|> Int
    else 
        nmat = [reverse(nmat;dims=1); nmat]
    end
    
    int_vals = cheby_int(bands)
    int_vals_2z = cheby_int_2z(bands)
    int_vals_8z3 = cheby_int_8z3(bands)
    gap_vals = cheby_gap(bands)
    g = size(bands,1)-1

    if length(κvec) != length(χvec)
        throw(ArgumentError("Number of poles and norming constants must agree."))
    end

    function r(j)
        if j > (g+1)/2 #Σ₊
            if typevec[j] == 1 #T
                return z -> h(j-(g+1)/2)(z)/(√(bands[j,2]-z |> Complex)*√(z-bands[j,1] |> Complex))
            elseif typevec[j] == 2 #U
                return z -> h(j-(g+1)/2)(z)*(√(bands[j,2]-z |> Complex)*√(z-bands[j,1] |> Complex))
            elseif typevec[j] == 3 #V
                return z -> h(j-(g+1)/2)(z)*(√(z-bands[j,1] |> Complex)/√(bands[j,2]-z |> Complex))
            else #W
                return z -> h(j-(g+1)/2)(z)*(√(bands[j,2]-z |> Complex)/√(z-bands[j,1] |> Complex))
            end
        else
            return z -> r(g+2-j)(-z)
        end
    end

    cc(j) = (bands[j,1]+bands[j,2])/2
    rr(j) = circ_size*(bands[j,2]-bands[j,1])/2  

    gridmat = Array{Array{ComplexF64,1},2}(undef,g+1,2)
    for j = 1:g+1
        gridmat[j,1] = rr(j)*zgrid(nmat[j,1]).+cc(j).+eps()im
        gridmat[j,2] = M(bands[j,1],bands[j,2]).(Tgrid(nmat[j,2])) .|> Complex
    end

    basisperm1 = [2, 1, 3, 4]
    basisperm2 = [1, 2, 4, 3];
    Chebymat = Array{ChebyParams}(undef,g+1,2)
    for j = 1:g+1
        if j > (g+1)/2 #Σ₊
            Chebymat[j,1] = buildCheby(bands[j,1],bands[j,2],typevec[j]) #left column
            Chebymat[j,2] = buildCheby(bands[j,1],bands[j,2],typevec[j]+2*(typevec[j]%2)-1) #right column
        else #Σ₋
            Chebymat[j,1] = buildCheby(bands[j,1],bands[j,2],basisperm1[typevec[j]]) 
            Chebymat[j,2] = buildCheby(bands[j,1],bands[j,2],basisperm2[typevec[j]]) 
        end
    end

    ntot = sum(sum(nmat))
    nptot(j) = sum(sum(nmat[1:j,:]))
    
    #build top left corner
    #build the first row of the block system
    A₁₁ = zeros(ComplexF64,ntot,ntot); rhs₁₁=zeros(ComplexF64,ntot,1)
    for j = 1:g+1
        CM₊ = zeros(ComplexF64,nmat[j,1],ntot)
        CM₋ = zeros(ComplexF64,nmat[j,1],ntot)
        CM₁₊= CauchyMat(gridmat[j,1],nmat[j,1],cc(j),rr(j);flag=1)
        CM₁₋= CauchyMat(gridmat[j,1],nmat[j,1],cc(j),rr(j);flag=-1)
        CM₂ = CauchyInterval(gridmat[j,1],Chebymat[j,1],nmat[j,2]-1)
        CM₊[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁₊
        CM₋[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁₋
        CM₊[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂
        CM₋[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,1],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,1],Chebymat[k,1],nmat[k,2]-1)
            CM₊[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]= CM₁
            CM₋[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]= CM₁
            CM₊[:,nptot(k)-nmat[k,2]+1:nptot(k)]= CM₂
            CM₋[:,nptot(k)-nmat[k,2]+1:nptot(k)]= CM₂
        end
        A₁ = CM₊-CM₋
        #build the second row of the block system
        CM₊ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₁ = CauchyMat(gridmat[j,2],nmat[j,1],cc(j),rr(j))
        CM₂₊= CauchyInterval(gridmat[j,2],Chebymat[j,1],nmat[j,2]-1;flag=1)
        CM₊[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁
        CM₊[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂₊
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,2],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,1],nmat[k,2]-1)
            CM₊[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]=CM₁
            CM₊[:,nptot(k)-nmat[k,2]+1:nptot(k)]=CM₂
        end
        #println(maximum(abs.(CM₊-CM₋)))
        A₂ = CM₊
        #global testguy=A₂
        #combine and build RHS
        A₁₁[nptot(j-1)+1:nptot(j),:] = [A₁; A₂]
        rhs₁₁[nptot(j-1)+1:nptot(j)] = [zeros(nmat[j,1]); -ones(nmat[j,2])]
    end
    
    Cmats = Array{Matrix{ComplexF64}}(undef,g+1,4)
    for j = 1:g+1
        #build portions of the first row of the block system
        CM₋ = zeros(ComplexF64,nmat[j,1],ntot)
        CM₁₋= CauchyMat(gridmat[j,1],nmat[j,1],cc(j),rr(j);flag=-1)
        CM₂ = CauchyInterval(gridmat[j,1],Chebymat[j,2],nmat[j,2]-1)
        CM₋[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁₋
        CM₋[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,1],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,1],Chebymat[k,2],nmat[k,2]-1)
            CM₋[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]= CM₁
            CM₋[:,nptot(k)-nmat[k,2]+1:nptot(k)]= CM₂
        end
        Cmats[j,1] = CM₋
        #build the second row of the block system
        CM₋ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₁ = CauchyMat(gridmat[j,2],nmat[j,1],cc(j),rr(j))
        CM₂₋= CauchyInterval(gridmat[j,2],Chebymat[j,2],nmat[j,2]-1;flag=-1)
        CM₋[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁
        CM₋[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂₋
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,2],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,2],nmat[k,2]-1)
            CM₋[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]=CM₁
            CM₋[:,nptot(k)-nmat[k,2]+1:nptot(k)]=CM₂
        end
        Cmats[j,2] = CM₋
    end

    #build bottom left corner
    #build the first row of the block system
    for j = 1:g+1
        CM₋ = zeros(ComplexF64,nmat[j,1],ntot)
        CM₁₋= CauchyMat(gridmat[j,1],nmat[j,1],cc(j),rr(j);flag=-1)
        CM₂ = CauchyInterval(gridmat[j,1],Chebymat[j,1],nmat[j,2]-1)
        CM₋[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁₋
        CM₋[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,1],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,1],Chebymat[k,1],nmat[k,2]-1)
            CM₋[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]= CM₁
            CM₋[:,nptot(k)-nmat[k,2]+1:nptot(k)]= CM₂
        end
        Cmats[j,3] = CM₋

        #build the second row of the block system
        CM₋ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₁ = CauchyMat(gridmat[j,2],nmat[j,1],cc(j),rr(j))
        CM₂₋ = CauchyInterval(gridmat[j,2],Chebymat[j,1],nmat[j,2]-1;flag=-1)
        CM₋[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁
        CM₋[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂₋
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,2],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,1],nmat[k,2]-1)
            CM₋[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]= CM₁
            CM₋[:,nptot(k)-nmat[k,2]+1:nptot(k)]= CM₂
        end
        Cmats[j,4] = CM₋
    end
    
    #build bottom right corner
    #build the first row of the block system
    A₂₂ = zeros(ComplexF64,ntot,ntot); rhs₂₂=zeros(ComplexF64,ntot,1)
    for j = 1:g+1
        CM₊ = zeros(ComplexF64,nmat[j,1],ntot)
        CM₋ = zeros(ComplexF64,nmat[j,1],ntot)
        CM₁₊ = CauchyMat(gridmat[j,1],nmat[j,1],cc(j),rr(j);flag=1)
        CM₁₋ = CauchyMat(gridmat[j,1],nmat[j,1],cc(j),rr(j);flag=-1)
        CM₂ = CauchyInterval(gridmat[j,1],Chebymat[j,2],nmat[j,2]-1)
        CM₊[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁₊
        CM₋[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]= CM₁₋
        CM₊[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂
        CM₋[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,1],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,1],Chebymat[k,2],nmat[k,2]-1)
            CM₊[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]=CM₁
            CM₋[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]=CM₁
            CM₊[:,nptot(k)-nmat[k,2]+1:nptot(k)]=CM₂
            CM₋[:,nptot(k)-nmat[k,2]+1:nptot(k)]=CM₂
        end
        gn₁ = ones(nmat[j,1])
        A₁ = CM₊-Diagonal(gn₁)*CM₋
        #build the second row of the block system
        CM₊ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₁ = CauchyMat(gridmat[j,2],nmat[j,1],cc(j),rr(j))
        CM₂₊ = CauchyInterval(gridmat[j,2],Chebymat[j,2],nmat[j,2]-1;flag=1)
        CM₊[:,nptot(j-1)+1:nptot(j)-nmat[j,2]]=CM₁
        CM₊[:,nptot(j)-nmat[j,2]+1:nptot(j)]= CM₂₊
        for k=(1:g+1)[1:end .!= j,:]
            CM₁ = CauchyMat(gridmat[j,2],nmat[k,1],cc(k),rr(k))
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,2],nmat[k,2]-1)
            CM₊[:,nptot(k-1)+1:nptot(k)-nmat[k,2]]= CM₁
            CM₊[:,nptot(k)-nmat[k,2]+1:nptot(k)]= CM₂
        end
        A₂ = CM₊
        #combine and build RHS
        A₂₂[nptot(j-1)+1:nptot(j),:] = [A₁; A₂]
        rhs₂₂[nptot(j-1)+1:nptot(j)] = [zeros(nmat[j,1]); -ones(nmat[j,2])]
    end
    
    Ω = g_coeffs(bands)
    #store terms for efficient computation of 𝔤(z) on circles and at soliton
    get_g_pre(z) = compute_g_pre(z,bands,int_vals_2z,int_vals_8z3,gap_vals)
    g_vals_pre = Array{Vector{Matrix{ComplexF64}}}(undef,g+2)
    for j = 1:g+1
        g_vals_pre[j] = map(λ -> get_g_pre(λ), gridmat[j,1])
    end
    g_vals_pre[g+2] = map(λ -> get_g_pre(λ), κvec)

    #store values of reflection coefficients and soliton transformation on intervals and circles
    m_r(j) = z -> imag(z)>0 ? 1/r(j)(z) : -1/r(j)(z)
    r_vals = Array{Vector{ComplexF64}}(undef,g+1,3)
    v_vals = Array{Vector{ComplexF64}}(undef,g+1,2)
    for j = 1:g+1
        r_vals[j,1] = map(z->r(j)(z), gridmat[j,2])
        r_vals[j,2] = map(z->1/r(j)(z), gridmat[j,2])
        r_vals[j,3] = map(z->m_r(j)(z), gridmat[j,1])
    end

    #store terms for efficient computation of 𝔥(z)
    (Ah_pre,Bh_pre) = h_coeffs_pre(bands)

    get_h_pre(z) = compute_h_pre(z, bands, int_vals, gap_vals)
    h_vals_pre = Array{Vector{Matrix{ComplexF64}}}(undef,g+2)
    for j = 1:g+1
        h_vals_pre[j] = map(z->get_h_pre(z), gridmat[j,1])
    end
    h_vals_pre[g+2] = map(z->get_h_pre(z), κvec)
    
    #get everything for undeformed problem
    Chebymatp = Array{ChebyParams}(undef,g+1,2)
    for j = 1:g+1
        if j > (g+1)/2 #Σ₊
            Chebymatp[j,1] = buildCheby(bands[j,1],bands[j,2],typevec[j]) #left column
            Chebymatp[j,2] = buildCheby(bands[j,1],bands[j,2],2) #right column
        else #Σ₋
            Chebymatp[j,1] = buildCheby(bands[j,1],bands[j,2],2) 
            Chebymatp[j,2] = buildCheby(bands[j,1],bands[j,2],basisperm2[typevec[j]]) 
        end
    end
    
    ntot = sum(nmat[:,2])
    nptot1(j) = sum(nmat[1:j,2])
    A₁₁p = zeros(ComplexF64,ntot,ntot); #rhs₁₁=zeros(ComplexF64,ntot,1)
    for j = 1:g+1
        CM₊ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₋ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₂₊= CauchyInterval(gridmat[j,2],Chebymat[j,1],nmat[j,2]-1;flag=1)
        CM₂₋= CauchyInterval(gridmat[j,2],Chebymat[j,1],nmat[j,2]-1;flag=-1)
        CM₊[:,nptot1(j-1)+1:nptot1(j)]= CM₂₊
        CM₋[:,nptot1(j-1)+1:nptot1(j)]= CM₂₋
        for k=(1:g+1)[1:end .!= j,:]
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,1],nmat[k,2]-1)
            CM₊[:,nptot1(k-1)+1:nptot1(k)] = CM₂
            CM₋[:,nptot1(k-1)+1:nptot1(k)] = CM₂
        end
        A₂ = CM₊-CM₋
        #combine and build RHS
        A₁₁p[nptot1(j-1)+1:nptot1(j),:] = A₂
    end
    
    Cmatsp = Array{Matrix{ComplexF64}}(undef,g+1,2)
    for j = 1:g+1
        CM₋ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₂₋= CauchyInterval(gridmat[j,2],Chebymat[j,2],nmat[j,2]-1;flag=-1)
        CM₋[:,nptot1(j-1)+1:nptot1(j)]= CM₂₋
        for k=(1:g+1)[1:end .!= j,:]
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,2],nmat[k,2]-1)
            CM₋[:,nptot1(k-1)+1:nptot1(k)]=CM₂
        Cmatsp[j,1] = CM₋
        end
    end
    
    for j = 1:g+1
        CM₋ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₂₋ = CauchyInterval(gridmat[j,2],Chebymat[j,1],nmat[j,2]-1;flag=-1)
        CM₋[:,nptot1(j-1)+1:nptot1(j)]= CM₂₋
        for k=(1:g+1)[1:end .!= j,:]
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,1],nmat[k,2]-1)
            CM₋[:,nptot1(k-1)+1:nptot1(k)]= CM₂
        end
        Cmatsp[j,2] = CM₋
    end

    A₂₂p = zeros(ComplexF64,ntot,ntot);
    for j = 1:g+1
        CM₊ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₋ = zeros(ComplexF64,nmat[j,2],ntot)
        CM₂₊ = CauchyInterval(gridmat[j,2],Chebymat[j,2],nmat[j,2]-1;flag=1)
        CM₂₋ = CauchyInterval(gridmat[j,2],Chebymat[j,2],nmat[j,2]-1;flag=-1)
        CM₊[:,nptot1(j-1)+1:nptot1(j)]= CM₂₊
        CM₋[:,nptot1(j-1)+1:nptot1(j)]= CM₂₋
        for k=(1:g+1)[1:end .!= j,:]
            CM₂ = CauchyInterval(gridmat[j,2],Chebymat[k,2],nmat[k,2]-1)
            CM₊[:,nptot1(k-1)+1:nptot1(k)] = CM₂
            CM₋[:,nptot1(k-1)+1:nptot1(k)] = CM₂
        end
        A₂ = CM₊-CM₋
        #combine and build RHS
        A₂₂p[nptot1(j-1)+1:nptot1(j),:] = A₂
    end

    up = undeformed_RHP_solitons(bands,A₁₁p,A₂₂p,Cmatsp,Chebymatp,nmat[:,2],gridmat[:,2],r_vals[:,1],κvec,χvec)
    dp = deformed_RHP_solitons(bands,A₁₁,A₂₂,rhs₁₁,rhs₂₂,Cmats,circ_size,Ω,Chebymat,nmat,gridmat,g_vals_pre,h_vals_pre,Ah_pre,Bh_pre,r_vals,κvec,χvec)
    rhp_solitons(dp,up)
end