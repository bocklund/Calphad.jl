using LinearAlgebra
using Symbolics

R = 8.3145

@variables N P T 
STATE_VARS = (T, P)

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

function c_iA(phase_record, i, A)
    total = 0.0
    state_vars_offset = length(phase_record.state_variables)
    for j in 1:length(phase_record.site_fractions)
        total += phase_record.inv_phase_matrix[i,j]*phase_record.mass_jac[A, state_vars_offset+j]
    end
    return total
end


function c_i_G(phase_record, i)
    total = 0.0
    state_vars_offset = length(phase_record.state_variables)
    for j in 1:length(phase_record.site_fractions)
        total += phase_record.inv_phase_matrix[i,j]*phase_record.grad[state_vars_offset+j]
    end
    return -total
end

function c_iG(phase_record, i)
    total = 0.0
    state_vars_offset = length(phase_record.state_variables)
    for j in 1:length(phase_record.site_fractions)
        total += phase_record.inv_phase_matrix[i,j]*phase_record.grad[state_vars_offset+j]
    end
    return -total
end

# TODO: implement
function c_iPot(phase_record, i, pot_idx)
    total = 0.0
    state_vars_offset = length(phase_record.state_variables)
    for j in 1:length(phase_record.site_fractions)
        total += phase_record.inv_phase_matrix[i,j]*phase_record.grad[state_vars_offset+j]
    end
    return -total
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
                # In Sundman's 2015 paper, this term is positive (Eq. 57), but
                # with this term positive, the chemical potentials are wrong
                # (exactly flipped for a [0.5, 0.5] binary system). I think the
                # error is introduced when he inserts Eq. 43 into Eq. 53 and the
                # signs are mistakenly flipped in Eq. 54. I suspect this is due
                # to redefining the c_iG term as a negative sum in Eq. 44 after
                # Eq. 54 and beyond had been derived. The mistake seems to be
                # propapaged through the rest of the paper after Eq. 54.
                total -= cs.ℵ * cs.phase_rec.mass_jac[A,state_vars_offset+_i] * c_i_G(cs.phase_rec, _i) 
                # TODO: N_A conditions need (N_A - Ñ_A) term added
            end
        end
        soln[P+A] = total
    end
    return soln
end

function get_subs_dict(compsets, conditions_dict)
d = Dict(conditions_dict...)
for compset in compsets
    for (Y_sym, Y_num) in zip(compset.phase_rec.site_fractions, compset.Y)
        d[Y_sym] = Y_num
    end
end
return d
end


"""
get_N_A_row_rhs
"""
# get_N_A_row_rhs
# Equation in LaTeX:
# \sum_\alpha \left[ \aleph^\alpha \sum_i \frac{\partial M_A^\alpha}{\partial y_i^\alpha} \left( \sum_B c_{iB} \mu_B + c_{iT} \Delta T + c_{iP} \Delta P \right) + M_A^\alpha \Delta \aleph^\alpha \right] = \sum_\alpha - c_{iG} + \left( N_A - \tilde{N_A} \right)
# This longer form is more representative of how we do the sum here internally in this function
# In this form, the left-most sum in each term on the left-hand-side corresponds
# to one column of the solution vector i.e.:
# 1. free chemical potentials for pure elements
# 2. free potentials
# 3. free phase amounts
# the right-most term in each sum is the term in the solution vector
# \sum_B \sum_\alpha \sum_i \aleph^\alpha \frac{\partial M_A^\alpha}{\partial y_i^\alpha}  c_{iB} \mu_B + \sum_\mathrm{Pot} \sum_\alpha \sum_i \aleph^\alpha \frac{\partial M_A^\alpha}{\partial y_i^\alpha} c_{i\mathrm{Pot}} \Delta \mathrm{Pot} + \sum_\beta \sum_\alpha M_A^\alpha \Delta \aleph^\beta  = \sum_\alpha - c_{iG} + \left( \sum_\alpha \aleph^\alpha M^\alpha_A - \tilde{N}_A \right)

compsets = [compset]
el_idx = 1
fixed_chempot_idxs = []
free_chempot_idxs = [1, 2]
fixed_pot_idxs = [1, 2]
free_pot_idxs = []
fixed_phase_idxs = []
free_phase_idxs = [1]

soln_size = (length(free_chempot_idxs) + length(free_pot_idxs) + length(free_phase_idxs))
row = Array{Num}(undef, soln_size)
# Chemical potential columns
for B in 1:length(free_chempot_idxs)
    total = 0.0
    for α in 1:length(compsets)
        cs = compsets[α]
        statevar_offset = length(cs.phase_rec.state_variables)
        for i in 1:length(cs.phase_rec.site_fractions)
            total += cs.ℵ * cs.phase_rec.mass_jac[el_idx,statevar_offset+i] * c_iA(cs.phase_rec,i,B)
        end
    end
    row[B] = total
end
# ΔPotential columns
col_offset = length(free_chempot_idxs)
for pot in 1:length(free_pot_idxs)
    total = 0.0
    for α in 1:length(compsets)
        cs = compsets[α]
        statevar_offset = length(cs.phase_rec.state_variables)
        for i in 1:length(cs.phase_recstate_variables)
            total += cs.ℵ * cs.phase_rec.mass_jac[el_idx,statevar_offset+i] * c_iPot(cs.phase_rec,i,pot)
        end
    end
    row[col_offset+pot] = total
end
# Δℵ columns
col_offset += length(free_pot_idxs)
for β in 1:length(free_phase_idxs)
    total = 0.0
    for α in 1:length(compsets)
        total += compsets[α].phase_rec.mass[el_idx]
    end
    row[col_offset+β] = total
end


# Condition on total number of moles, N
# Equation in LaTeX:
# \sum_A \sum_\alpha \left[ \aleph^\alpha \sum_i \frac{\partial M_A^\alpha}{\partial y_i^\alpha} \left( \sum_B c_{iB} \mu_B + c_{iT} \Delta T + c_{iP} \Delta P \right) + M_A^\alpha \Delta \aleph^\alpha \right] = \sum_A \sum_\alpha - c_{iG} + \left( N - \tilde{N} \right)




function solve(compsets, conditions)
    subs_dict = get_subs_dict(compsets, conditions)
    
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

# Super simple A-B system
@variables Y_BETA_A Y_BETA_B

@variables X_A X_B N_A N_B
@variables N_AL N_TI N_O X_AL X_TI X_O

G_BETA_A = 8000.0-10.0*T;
G_BETA_B = 12000.0-10.0*T;
G_BETA = Y_BETA_A*G_BETA_A + Y_BETA_B*G_BETA_B + R*T*(Y_BETA_A*log(Y_BETA_A) + Y_BETA_B*log(Y_BETA_B));
mass_BETA = [Y_BETA_A, Y_BETA_B];
state_variables = [T];
site_fractions = [Y_BETA_A, Y_BETA_B];
prx = PhaseRecord(G_BETA, mass_BETA, state_variables, site_fractions);

condition_dict = Dict(
N => 1.0,
P => 101325.0,
T => 300.0,
N_A => 0.5,
N_B => 0.5,
)
compset = CompSet(prx, [0.5, 0.5], 1.0);
soln = solve([compset], condition_dict)
# TODO: chemical potentials flipped, maybe a sign error
println("Actual chemical potentials: ", soln[1:2]);
println("Expected chemical potentials: ", [3271.04833019, 7271.04833015]);

condition_dict = Dict(
N => 1.0,
P => 101325.0,
T => 300.0,
N_A => 0.25,
N_B => 0.75,
)
compset = CompSet(prx, [0.25, 0.75], 1.0);
soln = solve([compset], Dict(T=>300.0))
# TODO: chemical potentials flipped, maybe a sign error
println("Actual chemical potentials: ", soln[1:2]);
println("Expected chemical potentials: ", [1542.0966603, 8282.42022259]);



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

@variables LIQUID0AL LIQUID0O LIQUID0TI T

G_M_LIQUID = -10000*LIQUID0AL + -40000*LIQUID0O + R*T*(LIQUID0AL*log(LIQUID0AL) + LIQUID0O*log(LIQUID0O))
mass = [LIQUID0AL, LIQUID0O]
state_variables = [T]
site_fractions = [LIQUID0AL, LIQUID0O]
prx = PhaseRecord(G_M_LIQUID, mass, state_variables, site_fractions)
compset = CompSet(prx, [0.5, 0.5], 1.0)

# Solve
A = get_equilibrium_matrix([compset]);
b = get_equilibrium_soln([compset]);
# Symbolic solution
x = A \ b;
subs_dict = Dict(Dict(zip(compset.phase_rec.site_fractions, compset.Y))..., T => 300.0)
substitute.(x, (subs_dict,))


############################################################





# Numerical linear least squares soln
AA = Symbolics.value.(substitute.(A, (subs_dict,)));
bb = Symbolics.value.(substitute.(b, (subs_dict,)));

LinearAlgebra.LAPACK.gelsd!(AA, bb)


