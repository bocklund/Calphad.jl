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
    state_vars_offset = length(phase_record.state_variables)
    for j in 1:length(phase_record.site_fractions)
        total += phase_record.inv_phase_matrix[i,j]*phase_record.mass_jac[el_index, state_vars_offset+j]
    end
    return total
end

function c_i_G(phase_record, i)
    total = 0.0
    state_vars_offset = length(phase_record.state_variables)
    for j in 1:length(phase_record.site_fractions)
        total -= phase_record.inv_phase_matrix[i,j]*phase_record.grad[state_vars_offset+j]
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
                state_vars_offset = length(compset.phase_rec.state_variables)
                term += compset.ℵ*sum(compset.phase_rec.mass_jac[i,state_vars_offset+idof]*c_i_el(compset.phase_rec, idof, j) for idof in 1:num_idof)
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
    for A in 1:N
        total = 0.0
        for ϕ in 1:P
            cs = compsets[ϕ]
            state_vars_offset = length(cs.phase_rec.state_variables)
            for _i in 1:length(cs.phase_rec.site_fractions)
                total += cs.ℵ * cs.phase_rec.mass_jac[A,state_vars_offset+_i] * c_i_G(cs.phase_rec, _i)
            end
        end
        soln[P+A] = total
    end
    return soln
end

function get_subs_dict(compsets, statevar_dict)
d = Dict(statevar_dict...)
for compset in compsets
    for (Y_sym, Y_num) in zip(compset.phase_rec.site_fractions, compset.Y)
        d[Y_sym] = Y_num
    end
end
return d
end


function solve(compsets, statevar_dict)
    subs_dict = get_subs_dict(compsets, statevar_dict)
    
    A = get_equilibrium_matrix(compsets);
    b = get_equilibrium_soln(compsets);
    
    AA = Symbolics.value.(substitute.(A, (subs_dict,)));
    bb = Symbolics.value.(substitute.(b, (subs_dict,)));
    
    # Numeric solution
    xx = AA \ bb;
    return xx
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

############################################################

# PHASE LIQUID % 1 1 !
#  CONSTITUENT LIQUID :AL,TI,O : !
#  PARAM G (LIQUID,AL;0) 1 -10000; 10000 N !
#  PARAM G (LIQUID,O;0) 1 -40000; 10000 N !
#  PARAM G (LIQUID,TI;0) 1 -20000; 10000 N !

@variables LIQUID0AL LIQUID0O LIQUID0TI

G_M_LIQUID = -10000*LIQUID0AL + -40000*LIQUID0O + R*T*(LIQUID0AL*log(LIQUID0AL) + LIQUID0O*log(LIQUID0O))
mass = [LIQUID0AL, LIQUID0O]
state_variables = [T]
site_fractions = [LIQUID0AL, LIQUID0O]
prx = PhaseRecord(G_M_LIQUID, mass, state_variables, site_fractions)
compset = CompSet(prx, [0.5, 0.5], 1.0)


############################################################


# Super simple A-B system
@variables Y_BETA_A Y_BETA_B T

G_BETA_A = 8000.0-10.0*T;
G_BETA_B = 12000.0-10.0*T;
G_BETA = Y_BETA_A*G_BETA_A + Y_BETA_B*G_BETA_B + R*T*(Y_BETA_A*log(Y_BETA_A) + Y_BETA_B*log(Y_BETA_B));
mass_BETA = [Y_BETA_A, Y_BETA_B];

state_variables = [T];
site_fractions = [Y_BETA_A, Y_BETA_B];

prx = PhaseRecord(G_BETA, mass_BETA, state_variables, site_fractions);
compset = CompSet(prx, [0.5, 0.5], 1.0);

subs_dict = Dict(Dict(zip(compset.phase_rec.site_fractions, compset.Y))..., T => 300.0)

A = get_equilibrium_matrix([compset]);
b = get_equilibrium_soln([compset]);
# Symbolic solution
x = A \ b;

# Numerical linear least squares soln
AA = Symbolics.value.(substitute.(A, (subs_dict,)));
bb = Symbolics.value.(substitute.(b, (subs_dict,)));

LinearAlgebra.LAPACK.gelsd!(AA, bb)


