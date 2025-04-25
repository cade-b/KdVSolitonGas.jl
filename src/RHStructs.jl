struct undeformed_RHP
    bands::Array{Float64,2}
    A₁₁::Array{ComplexF64}
    A₂₂::Array{ComplexF64}
    Cmats::Array{Matrix{ComplexF64}}
    Chebymat::Array{ChebyParams,2}
    nvec::Array{Int,1}
    gridpts::Array{Array{ComplexF64,1}}
    rvals::Array{Vector{ComplexF64},1}
end

struct deformed_RHP
    bands::Array{Float64,2}
    A₁₁::Array{ComplexF64}
    A₂₂::Array{ComplexF64}
    Cmats::Array{Matrix{ComplexF64}}
    circ_size::Float64
    Ω::Function
    Chebymat::Array{ChebyParams,2}
    nmat::Array{Int,2}
    gridpts::Array{Array{ComplexF64,1}}
    gvalsp::Array{Vector{Matrix{ComplexF64}}}
    hvalsp::Array{Vector{Matrix{ComplexF64}}}
    Ah_pre::Matrix{ComplexF64}
    Bh_pre::Matrix{ComplexF64}
    rvals::Matrix{Vector{ComplexF64}}
end

struct rhp
    dp::deformed_RHP
    up::undeformed_RHP
end

struct undeformed_RHP_soliton
    bands::Array{Float64,2}
    A₁₁::Array{ComplexF64}
    A₂₂::Array{ComplexF64}
    Cmats::Array{Matrix{ComplexF64}}
    Chebymat::Array{ChebyParams,2}
    nvec::Array{Int,1}
    gridpts::Array{Array{ComplexF64,1}}
    rvals::Array{Vector{ComplexF64},1}
    vvals::Array{Vector{ComplexF64},1}
    κ₀::Float64
    χ::Float64
end

struct deformed_RHP_soliton
    bands::Array{Float64,2}
    A₁₁::Array{ComplexF64}
    A₂₂::Array{ComplexF64}
    rhs₁₁::Array{ComplexF64}
    rhs₂₂::Array{ComplexF64}
    Cmats::Array{Matrix{ComplexF64}}
    circ_size::Float64
    Ω::Function
    Chebymat::Array{ChebyParams,2}
    nmat::Array{Int,2}
    gridpts::Array{Array{ComplexF64,1}}
    gvalsp::Array{Vector{Matrix{ComplexF64}}}
    hvalsp::Array{Vector{Matrix{ComplexF64}}}
    Ah_pre::Matrix{ComplexF64}
    Bh_pre::Matrix{ComplexF64}
    rvals::Matrix{Vector{ComplexF64}}
    vvals::Matrix{Vector{ComplexF64}}
    κ₀::Float64
    χ::Float64
end

struct rhp_soliton
    dp::deformed_RHP_soliton
    up::undeformed_RHP_soliton
end

struct undeformed_RHP_solitons
    bands::Array{Float64,2}
    A₁₁::Array{ComplexF64}
    A₂₂::Array{ComplexF64}
    Cmats::Array{Matrix{ComplexF64}}
    Chebymat::Array{ChebyParams,2}
    nvec::Array{Int,1}
    gridpts::Array{Array{ComplexF64,1}}
    rvals::Array{Vector{ComplexF64},1}
    κvec::Vector{Float64}
    χvec::Vector{Float64}
end

struct deformed_RHP_solitons
    bands::Array{Float64,2}
    A₁₁::Array{ComplexF64}
    A₂₂::Array{ComplexF64}
    rhs₁₁::Array{ComplexF64}
    rhs₂₂::Array{ComplexF64}
    Cmats::Array{Matrix{ComplexF64}}
    circ_size::Float64
    Ω::Function
    Chebymat::Array{ChebyParams,2}
    nmat::Array{Int,2}
    gridpts::Array{Array{ComplexF64,1},2}
    gvalsp::Array{Vector{Matrix{ComplexF64}}}
    hvalsp::Array{Vector{Matrix{ComplexF64}}}
    Ah_pre::Matrix{ComplexF64}
    Bh_pre::Matrix{ComplexF64}
    rvals::Matrix{Vector{ComplexF64}}
    κvec::Vector{Float64}
    χvec::Vector{Float64}
end

struct rhp_solitons
    dp::deformed_RHP_solitons
    up::undeformed_RHP_solitons
end