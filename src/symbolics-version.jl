using LinearAlgebra
using Symbolics

R = 8.3145

# Define molar Gibbs energy for (Al,Ti)2(O)3
@variables YM2O30AL YM2O30TI YM2O31O T

G_M_AL2O3 = -1772163.19+1053.4548*T-156.058*T*log(T)+.00709105*T^2-6.29402E-07*T^3 +12366650*T^(-1)
G_M_TI2O3 = -1581243.06+2395423.68*T^(-1)+940.164783*T-147.673862*T*log(T)-.00173711312*T^2-1.53383348E-10*T^3+19508-2.4214*T

G_M = G_M_AL2O3*YM2O30AL + G_M_TI2O3*YM2O30TI + R*T*(2*YM2O30AL*log(YM2O30AL) + 2*YM2O30TI*log(YM2O30TI))
G_m = G_M/(2*YM2O30AL + 2*YM2O30TI + 3*YM2O31O)

# Formula
mass = [2*YM2O30AL, 3*YM2O31O, 2*YM2O30TI] # moles(Al), moles(O), moles(Ti)

substitute(G_M, Dict(YM2O30AL => 0.5, YM2O30TI => 0.5, YM2O31O => 1.0, T => 2000.0))
substitute(G_M, Dict(YM2O30AL => 1.0, YM2O30TI => 1e-15, YM2O31O => 1.0, T => 2000.0))
substitute(G_m, Dict(YM2O30AL => 0.5, YM2O30TI => 0.5, YM2O31O => 1.0, T => 2000.0))

state_variables = [T]
site_fractions = [YM2O30AL, YM2O30TI, YM2O31O]
vars = vcat(state_variables, site_fractions)
mass_jac = Symbolics.jacobian(mass, site_fractions)

G_M_grad = Symbolics.gradient(G_M, vars)
G_M_Hess = Symbolics.hessian(G_M, vars)



function build_composition_set(G_M, state_variables, site_fractions)

end


function get_phase_matrix(G_M, site_fractions)
    num_internal_constraints = 1  # TODO: hardcoded for now
    num_internal_dof = length(site_fractions)
    N = num_internal_dof + num_internal_constraints
    phase_matrix = Matrix{Num}(undef, N, N)
    G_M_hess = Symbolics.hessian(G_M, site_fractions)
    for i in 1:num_internal_dof
        for j in 1:num_internal_dof
            phase_matrix[i,j] = G_M_hess[i,j]
        end
    end
    phase_matrix[N,:] .= 1.0
    phase_matrix[:,N] .= 1.0
    phase_matrix[N,N] = 0.0
    return phase_matrix
end


# # invert the phase matrix
# pm = get_phase_matrix(G_M, site_fractions)
# pm_inv = LinearAlgebra.inv(pm)

# e_11 = substitute(pm_inv[1,1], Dict(YM2O30AL => 0.5, YM2O30TI => 0.5, YM2O31O => 1.0, T => 2000.0))
# c_11 = pm_inv[1,1]*mass_jac[1,1]

# function c_i_el(i, inverted_phase_matrix, el_index, num_internal_dof)
#     return sum(pm_inv[i,j]*mass_jac[el_index, j] for j in 1:num_internal_dof)
# end


c_yal_al = c_i_el(1, pm_inv, 1, 3)

# struct CompositionSet
#     G_M
#     mass
#     mass_jac
#     ℵ
#     state_variables
#     site_fractions
# end
# cs = CompositionSet(G_M, mass, Symbolics.jacobian(mass, site_fractions), 0.2, state_variables, site_fractions)

# function get_equilibrium_matrix(inverted_phase_matrix, compsets)
#     N = length(compsets[1].mass)  # assume number of elements the same everywhere
#     P = length(compsets)
#     eq_mat = Matrix{Num}(undef, P+N, N+P)
#     for i in 1:P
#         for j in 1:N
#             eq_mat[i,j] = compsets[i].mass[j]
#         end
#         for j in N+1:N+P
#             eq_mat[i,j] = 0.0
#         end
#     end
#     for i in 1:N
#         for j in 1:N
#             term = 0.0
#             for phase_idx in 1:P
#                 compset = compsets[phase_idx]
#                 num_idof = length(compset.site_fractions)
#                 term += compset.ℵ*sum(compset.mass_jac[i,idof]*c_i_el(idof, inverted_phase_matrix, j, num_idof) for idof in 1:num_idof)
#             end
#             eq_mat[P+i,j] = term
#         end
#         for j in 1:P
#             eq_mat[P+i,N+j] = compsets[j].mass[i]
#         end
#     end

#     return eq_mat
# end

get_equilibrium_matrix(pm_inv, [cs])

struct PhaseRecord
    obj
    grad
    hess
    mass
    mass_jac
    state_variables
    site_fractions
    phase_matrix
    inv_phase_matrix
end

struct CompSet
    phase_rec::PhaseRecord
    Y
    ℵ
end

function phase_matrix(hess, num_state_variables, num_site_fractions)
    num_internal_constraints = 1  # TODO: hardcoded for now
    num_internal_dof = num_site_fractions
    phase_mat = Matrix{Num}(undef, num_internal_dof + num_internal_constraints, num_internal_dof + num_internal_constraints)
    for i in 1:num_internal_dof
        for j in 1:num_internal_dof
            phase_mat[i,j] = hess[num_state_variables+i,num_state_variables+j]
        end
    end
    for i in 1:num_internal_constraints
        for j in 1:num_internal_constraints
            phase_mat[num_internal_dof+i,:] .= 1.0
            phase_mat[:,num_internal_dof+j] .= 1.0
            phase_mat[num_internal_dof+i, num_internal_dof+i] = 0.0
        end
    end
    return phase_mat
end

function PhaseRecord(obj, mass, state_variables, site_fractions)
    grad = Symbolics.gradient(obj, [state_variables..., site_fractions...])
    hess = Symbolics.hessian(obj, [state_variables..., site_fractions...])
    mass_jac = Symbolics.jacobian(mass, [state_variables..., site_fractions...])
    pm = phase_matrix(hess, length(state_variables), length(site_fractions))

    return PhaseRecord(obj,
                       grad,
                       hess,
                       mass,
                       mass_jac,
                       state_variables,
                       site_fractions,
                       pm,
                       LinearAlgebra.inv(pm),
                       )
end

prx = PhaseRecord(G_M, mass, state_variables, site_fractions)
compset = CompSet(prx, [0.5, 0.5, 1.0], 0.2)

function c_i_el(i, phase_record, el_index)
    total = 0.0
    for j in 1:length(phase_record.site_fractions)
        total += phase_record.inv_phase_matrix[i,j]*phase_record.mass_jac[el_index, j]
    end
    return total
end

function get_equilibrium_matrix(compsets)
    N = length(compsets[1].phase_rec.mass)  # assume number of elements the same everywhere
    P = length(compsets)
    eq_mat = Matrix{Num}(undef, P+N, N+P)
    for i in 1:P
        for j in 1:N
            eq_mat[i,j] = compsets[i].phase_rec.mass[j]
        end
        for j in N+1:N+P
            eq_mat[i,j] = 0.0
        end
    end
    for i in 1:N
        for j in 1:N
            term = 0.0
            for phase_idx in 1:P
                compset = compsets[phase_idx]
                num_idof = length(compset.phase_rec.site_fractions)
                term += compset.ℵ*sum(compset.phase_rec.mass_jac[i,idof]*c_i_el(idof, compset.phase_rec, j) for idof in 1:num_idof)
            end
            eq_mat[P+i,j] = term
        end
        for j in 1:P
            eq_mat[P+i,N+j] = compsets[j].phase_rec.mass[i]
        end
    end
    return eq_mat
end

get_equilibrium_matrix([compset])
