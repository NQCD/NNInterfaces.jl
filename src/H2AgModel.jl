
using PeriodicTable
using H2AgModel_jll

export H2AgModel

"""
J. Phys. Chem. Lett. 2019, 10, 4962−4967
J. Phys. Chem. C 2020, 124, 186−195
"""
struct H2AgModel{T} <: AdiabaticFrictionModel
    h2indices::Vector{Int}
    tmp_coordinates::Matrix{T}
    tmp_friction_coordinates::Matrix{T}
    tmp_friction::Matrix{T}
    tmp_energy::Matrix{T}
    function H2AgModel(h2indices=[1, 2])

        cd(splitdir(H2AgModel_jll.h2ag111_pes_path)[1]) do
            ccall((:pes_init_, h2ag111_pes), Cvoid, ())
        end

        @assert h2indices == [1, 2]
        new{Float64}(h2indices, zeros(3, 2), zeros(3, 2), zeros(6, 6), zeros(1, 1))
    end
end

function NonadiabaticModels.potential!(model::H2AgModel, V::AbstractVector, R::AbstractMatrix)
    set_coordinates!(model, R)
    ccall(((:pot0_, h2ag111_pes)), Cvoid, (Ref{Int64}, Ref{Float64}, Ptr{Float64}),
            2, model.tmp_coordinates, model.tmp_energy)
    V[1] = austrip(model.tmp_energy[1]*u"eV")
end

function NonadiabaticModels.derivative!(model::H2AgModel, D::AbstractMatrix, R::AbstractMatrix)
    set_coordinates!(model, R)
    ccall((:dpeshon_, h2ag111_pes), Cvoid, (Ref{Int64}, Ref{Float64}, Ptr{Float64}),
            2, model.tmp_coordinates, D)
    @. D = austrip(D*u"eV/Å")
end

function NonadiabaticModels.friction!(model::H2AgModel, F::AbstractMatrix, R::AbstractMatrix)
    set_coordinates!(model, R)
    cd(splitdir(H2AgModel_jll.h2ag111_friction_path)[1]) do
        ccall((:tensor_, h2ag111_friction), Cvoid, (Ref{Float64}, Ptr{Float64}),
              model.tmp_coordinates, model.tmp_friction)
    end
    @. F[1:6,1:6] = austrip(model.tmp_friction*u"ps^-1")
    @. F *= austrip(elements[:H].atomic_mass)
end

function set_coordinates!(model::H2AgModel, R)
    @views model.tmp_coordinates .= ustrip.(auconvert.(u"Å", R[:,model.h2indices]))
    return nothing
end
