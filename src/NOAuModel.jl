
using Unitful: @u_str, ustrip
using UnitfulAtomic: austrip, auconvert
using PeriodicTable: PeriodicTable
using NonadiabaticModels: FrictionModels
using NOAuModel_jll

export NOAuModel

struct NOAuModel{T} <: FrictionModels.AdiabaticFrictionModel
    function NOAuModel()

        # cd(splitdir(NOAgModel_jll.NOAu111_pes_path)[1]) do
        #     ccall((:pes_init_, NOAu111_pes), Cvoid, ())
        # end
        @show NOAuModel_jll
        cd(splitdir(NOAuModel_jll.NOAu111_pes_path)[1]) do
            # EANN PES: call init_pes function
            ccall((:init_pes_, NOAu111_pes), Cvoid, ())
            # EANN PES: call init_pes_EFT function
            # ccall((:init_pes_eft_, "NOAu111_pes_EFT.dylib"), Cvoid, ())
        end

        # @assert h2indices == [1, 2]
        new{Float64}()
    end
end

# function NonadiabaticModels.potential(model::H2AgModel, R::AbstractMatrix)
#     set_coordinates!(model, R)
#     ccall(((:pot0_, h2ag111_pes)), Cvoid, (Ref{Int64}, Ref{Float64}, Ptr{Float64}),
#             2, model.tmp_coordinates, model.tmp_energy)
#     return austrip(model.tmp_energy[1]*u"eV")
# end

# function NonadiabaticModels.derivative!(model::H2AgModel, D::AbstractMatrix, R::AbstractMatrix)
#     set_coordinates!(model, R)
#     ccall((:dpeshon_, H2AgModel_jll.h2ag111_pes), Cvoid, (Ref{Int64}, Ref{Float64}, Ptr{Float64}),
#             2, model.tmp_coordinates, D)
#     @. D = austrip(D*u"eV/Å")
#     return D
# end

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

# function set_coordinates!(model::H2AgModel, R)
#     @views model.tmp_coordinates .= ustrip.(auconvert.(u"Å", R[:,model.h2indices]))
#     return nothing
# end
