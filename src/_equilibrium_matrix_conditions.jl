# Included by symbolic_solver.jl

function c_iG(phase_record, i)
    # c^\alpha_{iG} = -\sum_j e^\alpha_{ij} \frac{\partial G^\alpha_M}{\partial y^\alpha_j}
    total = 0.0
    state_vars_offset = length(phase_record.state_variables)
    for j in 1:length(phase_record.site_fractions)
        total += phase_record.inv_phase_matrix[i,j]*phase_record.grad[state_vars_offset+j]
    end
    return -total
end

function c_iPot(phase_record, i, pot_idx)
    # c^\alpha_{i\mathrm{Pot}} = -\sum_j e^\alpha_{ij} \frac{\partial^2 G^\alpha_M}{\partial \mathrm{Pot} \partial y^\alpha_j}
    total = 0.0
    state_vars_offset = length(phase_record.state_variables)
    for j in 1:length(phase_record.site_fractions)
        total += phase_record.inv_phase_matrix[i,j]*phase_record.hess[pot_idx,state_vars_offset+j]
    end
    return -total
end

function c_iA(phase_record, i, A)
    # c^\alpha_{iA} = \sum_j e^\alpha_{ij} \frac{\partial M^\alpha_A}{\partial y^\alpha_j}
    total = 0.0
    state_vars_offset = length(phase_record.state_variables)
    for j in 1:length(phase_record.site_fractions)
        total += phase_record.inv_phase_matrix[i,j]*phase_record.mass_jac[A, state_vars_offset+j]
    end
    return total
end

function _ℵ(phase_record)
    # TODO: This will break if two phase records have the same name (miscibility gap)
    str = "ℵ_$(phase_record.phase_name)"
    return Num(Variable(Symbol(str)))
end

# LaTeX equation:
# \sum_{B_{\mathrm{free}}} M_B^\alpha \mu_B + \sum_\mathrm{Pot} -\frac{\partial G_M^\alpha}{\partial \mathrm{Pot}} \Delta \mathrm{Pot} + \sum_\beta 0 = G_M^\alpha + \sum_{B_{\mathrm{fixed}}} -M_B^\alpha \mu_B
"""

# Examples
```
srow, srhs = get_stable_phase_row_rhs([prx], 1, [], [1, 2], [P, T], [], [], [1])
```
"""
function get_stable_phase_row_rhs(phase_records, phase_idx,
                                  fixed_chempot_symbols, free_chempot_idxs,
                                  fixed_pot_symbols, free_pot_idxs,
                                  fixed_phase_symbols, free_phase_idxs,
                                  )
    phase_rec = phase_records[phase_idx]
    # Construct the row in the equilibrium matrix
    soln_size = (length(free_chempot_idxs) + length(free_pot_idxs) + length(free_phase_idxs))
    row = Array{Num}(undef, soln_size)
    for B in 1:length(free_chempot_idxs)
        row[B] = phase_rec.mass[B]
    end
    # ΔPotential columns
    col_offset = length(free_chempot_idxs)
    for pot in 1:length(free_pot_idxs)
        row[col_offset+pot] = -phase_rec.grad[pot]
    end
    # Δℵ columns
    col_offset += length(free_pot_idxs)
    for β in 1:length(free_phase_idxs)
        row[col_offset+β] = 0.0
    end

    # Construct right-hand-side (RHS)
    rhs = 0.0
    rhs += phase_rec.obj
    for B in 1:length(fixed_chempot_symbols)
        rhs -= phase_rec.mass[B] * fixed_chempot_symbols[B]
    end

    return (row, rhs)
end


# Equation in LaTeX:
# This longer form is more representative of how we do the sum here internally in this function
# In this form, the left-most sum in each term on the left-hand-side corresponds
# to one column of the solution vector i.e.:
# 1. free chemical potentials for pure elements
# 2. free potentials
# 3. free phase amounts
# the right-most term in each sum is the term in the solution vector
# \sum_{B_{\mathrm{free}}} \sum_\alpha \sum_i \aleph^\alpha \frac{\partial M_A^\alpha}{\partial y_i^\alpha}  c_{iB} \mu_B + \sum_\mathrm{Pot} \sum_\alpha \sum_i \aleph^\alpha \frac{\partial M_A^\alpha}{\partial y_i^\alpha} c_{i\mathrm{Pot}} \Delta \mathrm{Pot} + \sum_\beta \sum_\alpha M_A^\alpha \Delta \aleph^\beta  \\
# = \sum_\alpha \sum_i - \aleph^\alpha \frac{\partial M_A^\alpha}{\partial y_i^\alpha} c_{iG} + \sum_{B_{\mathrm{fixed}}} \sum_\alpha \sum_i - \aleph^\alpha \frac{\partial M_A^\alpha}{\partial y_i^\alpha}  c_{iB} \mu_B + \left( \sum_\alpha \aleph^\alpha M^\alpha_A - \tilde{N}_A \right)
"""
get_N_A_row_rhs

# Examples
```
using Symbolics
@variables T P N_A
@variables Y_BETA_A Y_BETA_B
G_BETA_A = 8000.0-10.0*T;
G_BETA_B = 12000.0-10.0*T;
G_BETA = Y_BETA_A*G_BETA_A + Y_BETA_B*G_BETA_B + R*T*(Y_BETA_A*log(Y_BETA_A) + Y_BETA_B*log(Y_BETA_B));
mass_BETA = [Y_BETA_A, Y_BETA_B];
state_variables = [P, T];
site_fractions = [Y_BETA_A, Y_BETA_B];
prx = PhaseRecord("BETA", G_BETA, mass_BETA, state_variables, site_fractions);

r, rrhs = get_N_A_row_rhs([prx], 1, N_A, [], [1, 2], [P, T], [], [], [1])
```
"""
function get_N_A_row_rhs(phase_records, el_idx, N_el_sym,
                         fixed_chempot_symbols, free_chempot_idxs,
                         fixed_pot_symbols, free_pot_idxs,
                         fixed_phase_symbols, free_phase_idxs,
                         )
    # Construct the row in the equilibrium matrix
    soln_size = (length(free_chempot_idxs) + length(free_pot_idxs) + length(free_phase_idxs))
    row = Array{Num}(undef, soln_size)
    # Chemical potential columns
    for B in 1:length(free_chempot_idxs)
        total = 0.0
        for α in 1:length(phase_records)
            prx = phase_records[α]
            statevar_offset = length(prx.state_variables)
            for i in 1:length(prx.site_fractions)
                total += _ℵ(prx) * prx.mass_jac[el_idx,statevar_offset+i] * c_iA(prx,i,B)
            end
        end
        row[B] = total
    end
    # ΔPotential columns
    col_offset = length(free_chempot_idxs)
    for pot in 1:length(free_pot_idxs)
        total = 0.0
        for α in 1:length(phase_records)
            prx = phase_records[α]
            statevar_offset = length(prx.state_variables)
            for i in 1:length(prx.state_variables)
                total += _ℵ(prx) * prx.mass_jac[el_idx,statevar_offset+i] * c_iPot(prx,i,pot)
            end
        end
        row[col_offset+pot] = total
    end
    # Δℵ columns
    col_offset += length(free_pot_idxs)
    for β in 1:length(free_phase_idxs)
        total = 0.0
        for α in 1:length(phase_records)
            prx = phase_records[α]
            total += prx.mass[el_idx]
        end
        row[col_offset+β] = total
    end

    # construct the right-hand-side term
    rhs = 0.0
    for α in 1:length(phase_records)
        prx = phase_records[α]
        statevar_offset = length(prx.state_variables)
        for i in 1:length(prx.state_variables)
            rhs -= _ℵ(prx) * prx.mass_jac[el_idx,statevar_offset+i] * c_iG(prx, i)
        end
    end
    for B in 1:length(fixed_chempot_symbols)
        for α in 1:length(phase_records)
            prx = phase_records[α]
            statevar_offset = length(prx.state_variables)
            for i in 1:length(prx.state_variables)
                rhs -= _ℵ(prx) * prx.mass_jac[el_idx,statevar_offset+i] * fixed_chempot_symbols[B]
            end
        end
    end
    for α in 1:length(phase_records)
        prx = phase_records[α]
        rhs += prx.mass[el_idx]
    end
    rhs -= N_el_sym

    return (row, rhs)
end

"""
    get_N_row_rhs

See the function `get_N_A_row_rhs`. This is conceptually the same as an N_A
condition with an inner loop over all the elements in each column.
"""
function get_N_row_rhs(phase_records, N_sym,
                       fixed_chempot_symbols, free_chempot_idxs,
                       fixed_pot_symbols, free_pot_idxs,
                       fixed_phase_symbols, free_phase_idxs,
                       )
    # Construct the row in the equilibrium matrix
    soln_size = (length(free_chempot_idxs) + length(free_pot_idxs) + length(free_phase_idxs))
    N_elements = length(phase_records[1].mass)  # assume mass length is the same for every phase record
    row = Array{Num}(undef, soln_size)
    # Chemical potential columns
    for B in 1:length(free_chempot_idxs)
        total = 0.0
        for A in 1:N_elements
            for α in 1:length(phase_records)
                prx = phase_records[α]
                statevar_offset = length(prx.state_variables)
                for i in 1:length(prx.site_fractions)
                    total += _ℵ(prx) * prx.mass_jac[A,statevar_offset+i] * c_iA(prx,i,B)
                end
            end
        end
        row[B] = total
    end
    # ΔPotential columns
    col_offset = length(free_chempot_idxs)
    for pot in 1:length(free_pot_idxs)
        total = 0.0
        for A in 1:N_elements
            for α in 1:length(phase_records)
                prx = phase_records[α]
                statevar_offset = length(prx.state_variables)
                for i in 1:length(prx.state_variables)
                    total += _ℵ(prx) * prx.mass_jac[A,statevar_offset+i] * c_iPot(prx,i,pot)
                end
            end
        end
        row[col_offset+pot] = total
    end
    # Δℵ columns
    col_offset += length(free_pot_idxs)
    for β in 1:length(free_phase_idxs)
        total = 0.0
        for A in 1:N_elements
            for α in 1:length(phase_records)
                prx = phase_records[α]
                total += prx.mass[A]
            end
        end
        row[col_offset+β] = total
    end

    # construct the right-hand-side term
    rhs = 0.0
    for A in 1:N_elements
        for α in 1:length(phase_records)
            prx = phase_records[α]
            statevar_offset = length(prx.state_variables)
            for i in 1:length(prx.state_variables)
                rhs -= _ℵ(prx) * prx.mass_jac[A,statevar_offset+i] * c_iG(prx, i)
            end
        end
        for B in 1:length(fixed_chempot_symbols)
            for α in 1:length(phase_records)
                prx = phase_records[α]
                statevar_offset = length(prx.state_variables)
                for i in 1:length(prx.state_variables)
                    rhs -= _ℵ(prx) * prx.mass_jac[A,statevar_offset+i] * fixed_chempot_symbols[B]
                end
            end
        end
        for α in 1:length(phase_records)
            prx = phase_records[α]
            rhs += prx.mass[A]
        end
    end
    rhs -= N_sym

    return (row, rhs)
end

# This is the only "public" function of this file.
"""
    cond_row_rhs

This function takes a condition symbol and dispatches on the correct method that
condition that returns the elements in the row and the right hand side.

The number of columns are determined by the `fixed_free_terms` - there's one
column for each free chemical potential, free potential, and free phase amount.
"""
function cond_row_rhs(cond, elements, phases, phase_records, fixed_free_terms)
    # assumes elements, phases, phase_records are sorted
    str_cond = string(cond)
    if str_cond == "N"
        row, rhs = get_N_row_rhs(phase_records, cond, fixed_free_terms...)
    elseif startswith(str_cond, "N_")
        el = str_cond[3:end]
        el_idx = findfirst(x -> x == el, elements)
        row, rhs = get_N_A_row_rhs(phase_records, el_idx, cond, fixed_free_terms...)
    else
        # TODO: implement X(A) condition row
        throw("Condition for $str_cond is not yet implemented")
    end

    return row, rhs
end