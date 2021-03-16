using LinearAlgebra
using Symbolics

R = 8.3145

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

function c_i_el(phase_record, i, el_index)
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
                term += compset.ℵ*sum(compset.phase_rec.mass_jac[i,idof]*c_i_el(compset.phase_rec, idof, j) for idof in 1:num_idof)
            end
            eq_mat[P+i,j] = term
        end
        for j in 1:P
            eq_mat[P+i,N+j] = compsets[j].phase_rec.mass[i]
        end
    end
    return eq_mat
end

function get_equilibrium_soln(compsets)
    N = length(compsets[1].phase_rec.mass)  # assume number of elements the same everywhere
    P = length(compsets)
    soln = Array{Num}(undef, P+N)
    for i in 1:P
        soln[i] = compsets[i].phase_rec.obj
    end
    for i in 1:N
        total = 0.0
        for ϕ in 1:P
            cs = compsets[ϕ]
            for j in 1:length(cs.phase_rec.site_fractions)
                total += cs.ℵ * c_i_el(cs.phase_rec, j, i)
            end
        end
        soln[P+i] = total
    end
    return soln
end


######################################################
# Script


# Define molar Gibbs energy for (Al,Ti)2(O)3
@variables YM2O30AL YM2O30TI YM2O31O T

G_M_AL2O3 = -1772163.19+1053.4548*T-156.058*T*log(T)+.00709105*T^2-6.29402E-07*T^3 +12366650*T^(-1)
G_M_TI2O3 = -1581243.06+2395423.68*T^(-1)+940.164783*T-147.673862*T*log(T)-.00173711312*T^2-1.53383348E-10*T^3+19508-2.4214*T

G_M = G_M_AL2O3*YM2O30AL + G_M_TI2O3*YM2O30TI + R*T*(2*YM2O30AL*log(YM2O30AL) + 2*YM2O30TI*log(YM2O30TI) + 3*YM2O31O*log(YM2O31O))
G_m = G_M/(2*YM2O30AL + 2*YM2O30TI + 3*YM2O31O)
# Formula
mass = [2*YM2O30AL, 3*YM2O31O, 2*YM2O30TI] # moles(Al), moles(O), moles(Ti)

# substitute(G_M, Dict(YM2O30AL => 0.5, YM2O30TI => 0.5, YM2O31O => 1.0, T => 2000.0))
# substitute(G_M, Dict(YM2O30AL => 1.0, YM2O30TI => 1e-15, YM2O31O => 1.0, T => 2000.0))
# substitute(G_m, Dict(YM2O30AL => 0.5, YM2O30TI => 0.5, YM2O31O => 1.0, T => 2000.0))

state_variables = [T]
site_fractions = [YM2O30AL, YM2O30TI, YM2O31O]

prx = PhaseRecord(G_M, mass, state_variables, site_fractions)
compset = CompSet(prx, [0.5, 0.5, 1.0], 0.2)

A = get_equilibrium_matrix([compset])
b = get_equilibrium_soln([compset])
x = A \ b
subs_dict = Dict(Dict(zip(compset.phase_rec.site_fractions, compset.Y))..., T => 2000.0)
substitute.(x, (subs_dict,))

# Using lapack least squares, need to convert A and b into floats first

AA = Symbolics.value.(substitute.(A, (subs_dict,)))
bb = Symbolics.value.(substitute.(b, (subs_dict,)))
LinearAlgebra.LAPACK.gelsd!(AA, bb)