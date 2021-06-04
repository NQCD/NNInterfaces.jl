
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

        cd("$(H2AgModel_jll.artifact_dir)/lib") do
            ccall((:pes_init_, h2ag111_pes), Cvoid, ())
        end

        @assert h2indices == [1, 2]
        new{Float64}(h2indices, zeros(3, 2), zeros(3, 2), zeros(6, 6), zeros(1, 1))
    end
end

# function potential!(model::EANN_H₂Ag, V::AbstractVector, R::AbstractMatrix)
#     @views model.tmp_coordinates .= au_to_ang.(R[:,model.h2indices])
#     cd(splitdir(model.potential_path)[1]) do
#         ccall(model.potential_function, Cvoid, (Ref{Int64}, Ref{Float64}, Ptr{Float64}),
#               2, model.tmp_coordinates, model.tmp_energy)
#     end
#     V[1] = eV_to_au(model.tmp_energy[1])
# end

# function derivative!(model::EANN_H₂Ag, D::AbstractMatrix, R::AbstractMatrix)
#     @views model.tmp_coordinates .= au_to_ang.(R[:,model.h2indices])
#     cd(splitdir(model.potential_path)[1]) do
#         ccall(model.force_function, Cvoid, (Ref{Int64}, Ref{Float64}, Ptr{Float64}),
#               2, model.tmp_coordinates, D)
#     end
#     D .= eV_per_ang_to_au.(D)
# end

# function friction!(model::EANN_H₂Ag, F::AbstractMatrix, R::AbstractMatrix)
#     @views model.tmp_friction_coordinates .= au_to_ang.(R[:,model.h2indices])
#     cd(splitdir(model.friction_path)[1]) do
#         ccall(model.friction_function, Cvoid, (Ref{Float64}, Ptr{Float64}),
#               model.tmp_coordinates, model.tmp_friction)
#     end
#     F[1:6,1:6] .= ps_inv_to_au.(model.tmp_friction)
#     F .*= austrip(elements[:H].atomic_mass)
# end
