
using Unitful: @u_str, ustrip
using UnitfulAtomic: austrip, auconvert
using PeriodicTable: PeriodicTable
using NonadiabaticModels: FrictionModels
using NOAuModel_jll

export NOAuModel



# Constants extracted from source code. Comments refer to file locations of constants
const nkpoint = 3 # EFT_input
const outputneuron = 3 # EFT_input
const atomdim = 3 # mod.f90, 3 dofs for each atom
const numatom = 38 # EFT_input, total number of atoms
const neff = 2 # EFT_input, think this is effective number of atoms, so just N and O
const numforce = neff * atomdim # 6, 3 force components for each atom

const symbols = vcat([:O, :N], fill(:Au, 36)) # Correct format for the system, atom

"""
J. Phys. Chem. Lett. 2019, 10, 4962−4967
J. Phys. Chem. C 2020, 124, 186−195
"""
struct NOAuModel{T} <: FrictionModels.AdiabaticFrictionModel
    # h2indices::Vector{Int}
    tmp_coordinates::Matrix{T}
    # tmp_friction_coordinates::Matrix{T}
    # tmp_friction::Matrix{T}
    tmp_energy::Matrix{T}
    tmp_force::Array{T,3}
    function NOAuModel()

        cd(splitdir(NOAuModel_jll.NOAu111_pes_path)[1]) do
            ccall((:init_pes_eft_, NOAu111_pes), Cvoid, ())
        end

        new{Float64}(zeros(atomdim, numatom), zeros(outputneuron, nkpoint),
                     zeros(numforce, outputneuron, nkpoint))
    end
end

ndofs(::NOAuModel) = atomdim

function NonadiabaticModels.potential(model::NOAuModel, R::AbstractMatrix)
    set_coordinates!(model, R)
    ccall(((:get_energy_eft_, NOAu111_pes)), Cvoid, (Ptr{Float64}, Ptr{Float64}),
            model.tmp_coordinates, model.tmp_energy)
    return austrip.(model.tmp_energy.*u"eV")
end

function NonadiabaticModels.derivative!(model::NOAuModel, D::AbstractMatrix, R::AbstractMatrix)
    set_coordinates!(model, R)
    ccall((:get_force_eft_, NOAu111_pes), Cvoid, (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
            model.tmp_coordinates, model.tmp_energy, model.tmp_force)
    # @. D = austrip(D*u"eV/Å")
    return model.tmp_force
end

# function FrictionModels.friction!(model::H2AgModel, F::AbstractMatrix, R::AbstractMatrix)
#     set_coordinates!(model, R)
#     cd(splitdir(H2AgModel_jll.h2ag111_friction_path)[1]) do
#         ccall((:tensor_, h2ag111_friction), Cvoid, (Ref{Float64}, Ptr{Float64}),
#               model.tmp_coordinates, model.tmp_friction)
#     end
#     @. F[1:6,1:6] = austrip(model.tmp_friction*u"ps^-1")
#     @. F *= austrip(PeriodicTable.elements[:H].atomic_mass)
#     return F
# end

function set_coordinates!(model::NOAuModel, R)
    @views model.tmp_coordinates .= ustrip.(auconvert.(u"Å", R))
    return nothing
end
