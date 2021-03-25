# Included by symbolic_solver.jl

# In this module are methods to compute rows of the equilibrium matrix.
# Each row method has it's equation documented in LaTeX. The more expanded
# notation is designed to better represent the way the code is written. In this
# form, the left-most sum index in each term on the left-hand-side corresponds
# to one column of the solution vector, i.e.:
# 1. free chemical potentials for pure elements
# 2. free potentials
# 3. free phase amounts
# Note that the right-most term in each sum is the term in the solution vector,
# i.e. μ_B, ΔPot, Δℵ

# TODO: all the row methods could use some simplification in the arguments to
# not depend on the fixed/free amounts that are not relevant. It would also be
# nice to make sure they work regardless of whether the inputs are symbolic or
# numeric.

# TODO: The residual terms for N, N_A, and X(B) subtract the current value from the
# prescribed value. It's not clear to me mathematically why this is, but it is easy to
# verify for an ideal system that is not mass balanced that this should be the case.

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

@doc raw"""

# Examples
```
srow, srhs = get_stable_phase_row_rhs([prx], 1, [], [1, 2], [P, T], [], [], [1])
```

# Equation

```math
\sum_{B_{\mathrm{free}}} M_B^\alpha \mu_B + \sum_\mathrm{Pot} -\frac{\partial G_M^\alpha}{\partial \mathrm{Pot}} \Delta \mathrm{Pot} + \sum_\beta 0 = G_M^\alpha + \sum_{B_{\mathrm{fixed}}} -M_B^\alpha \mu_B
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
        row[col_offset+pot] = -phase_rec.grad[free_pot_idxs[pot]]
    end
    # Δℵ columns
    col_offset += length(free_pot_idxs)
    for β in 1:length(free_phase_idxs)
        row[col_offset+β] = 0.0
    end

    # Construct right-hand-side (RHS)
    rhs = 0.0
    # Gibbs energy
    rhs += phase_rec.obj
    # Fixed chemical potentials
    for B in 1:length(fixed_chempot_symbols)
        rhs -= phase_rec.mass[B] * fixed_chempot_symbols[B]
    end

    return (row, rhs)
end


@doc raw"""
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

# Equation
```math
\sum_{B_{\mathrm{free}}} \sum_\alpha \sum_i \aleph^\alpha \frac{\partial M_A^\alpha}{\partial y_i^\alpha}  c_{iB} \mu_B + \sum_\mathrm{Pot} \sum_\alpha \sum_i \aleph^\alpha \frac{\partial M_A^\alpha}{\partial y_i^\alpha} c_{i\mathrm{Pot}} \Delta \mathrm{Pot} + \sum_\beta M_A^\beta \Delta \aleph^\beta  \\
= - \sum_\alpha \sum_i \aleph^\alpha \frac{\partial M_A^\alpha}{\partial y_i^\alpha} c_{iG} - \sum_{B_{\mathrm{fixed}}} \sum_\alpha \sum_i \aleph^\alpha \frac{\partial M_A^\alpha}{\partial y_i^\alpha}  c_{iB} \mu_B + \left( \tilde{N}_A - \sum_\alpha \aleph^\alpha M^\alpha_A \right)
```

"""
function get_N_A_row_rhs(phase_records, el_idx, N_A_prescribed,
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
                total += prx.ℵ * prx.mass_jac[el_idx,statevar_offset+i] * c_iA(prx,i,B)
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
                total += prx.ℵ * prx.mass_jac[el_idx,statevar_offset+i] * c_iPot(prx,i,free_pot_idxs[pot])
            end
        end
        row[col_offset+pot] = total
    end
    # Δℵ columns
    col_offset += length(free_pot_idxs)
    for β in 1:length(free_phase_idxs)
        row[col_offset+β] = phase_records[β].mass[el_idx]
    end

    # construct the right-hand-side term
    rhs = 0.0
    # C_iG terms
    for α in 1:length(phase_records)
        prx = phase_records[α]
        statevar_offset = length(prx.state_variables)
        for i in 1:length(prx.state_variables)
            rhs -= prx.ℵ * prx.mass_jac[el_idx,statevar_offset+i] * c_iG(prx, i)
        end
    end
    # Fixed chemical potentials
    for B in 1:length(fixed_chempot_symbols)
        for α in 1:length(phase_records)
            prx = phase_records[α]
            statevar_offset = length(prx.state_variables)
            for i in 1:length(prx.state_variables)
                rhs -= prx.ℵ * prx.mass_jac[el_idx,statevar_offset+i] * fixed_chempot_symbols[B]
            end
        end
    end
    # Prescribed N_A
    rhs += N_A_prescribed
    # System current N_A
    for α in 1:length(phase_records)
        prx = phase_records[α]
        rhs -= prx.mass[el_idx]
    end

    return (row, rhs)
end

@doc raw"""
    get_x_A_row_rhs

See the function `get_N_A_row_rhs`. This is the corresponding mole fraction condition.

# Equation

```math
% Free chemical potential columns (μ_B)
\sum_{B_{\mathrm{free}}} \sum_\alpha \sum_i \frac{\aleph^\alpha c_{iB}}{N} \left( \frac{\partial M_A^\alpha}{\partial y_i^\alpha} - x_A \sum_C \frac{\partial M_C^\alpha}{\partial y_i^\alpha} \right)  \mu_B & \\
% Free potential columns (ΔP)
+ \sum_\mathrm{Pot} \sum_\alpha \sum_i \frac{\aleph^\alpha c_{i\mathrm{Pot}}}{N} \left( \frac{\partial M_A^\alpha}{\partial y_i^\alpha} - x_A \sum_C \frac{\partial M_C^\alpha}{\partial y_i^\alpha} \right)  \Delta \mathrm{Pot} & \\
% Free phase amount columns (Δℵ)
+ \sum_{\beta_\mathrm{free}} \frac{1}{N} \left( M_A^\beta - x_A \sum_C M_C^\beta \right) \Delta \aleph^\beta &
% Right-hand side (constant)
% c_iG
= - \sum_\alpha \sum_i \frac{\aleph^\alpha c_{iG}}{N} \left( \frac{\partial M_A^\alpha}{\partial y_i^\alpha} - x_A \sum_C \frac{\partial M_C^\alpha}{\partial y_i^\alpha} \right) \\
% Fixed chemical potentials
& - \sum_{B_{\mathrm{fixed}}} \sum_\alpha \sum_i \frac{\aleph^\alpha c_{iB}}{N} \left( \frac{\partial M_A^\alpha}{\partial y_i^\alpha} - x_A \sum_C \frac{\partial M_C^\alpha}{\partial y_i^\alpha} \right)  \mu_B \\
% The term for the condition
& + \left( \tilde{x}_A - x_A \right) \\
```
"""
function get_x_A_row_rhs(phase_records, el_idx, x_A_prescribed,
                         fixed_chempot_symbols, free_chempot_idxs,
                         fixed_pot_symbols, free_pot_idxs,
                         fixed_phase_symbols, free_phase_idxs,
                         )
    num_elements = length(phase_records[1].mass)  # assumes number of elements

    # Constuct N and x_A, which will be used in almost all the expresions below
    N_A = 0.0  # Intermediate, just used for computing x_A
    for α in 1:length(phase_records)
        N_A += phase_records[α].ℵ * phase_records[α].mass[el_idx]
    end
    N = 0.0
    for B in 1:num_elements
        for α in 1:length(phase_records)
            N += phase_records[α].ℵ * phase_records[α].mass[B]
        end
    end
    x_A = N_A / N

    # TODO: conditions other than chemical potential
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
                factor = prx.ℵ * c_iA(prx,i,B) / N
                inner_total = prx.mass_jac[el_idx,statevar_offset+i]
                for C in 1:num_elements
                    inner_total -= x_A * prx.mass_jac[C,statevar_offset+i]
                end
                total += factor * inner_total
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
                factor = prx.ℵ * c_iPot(prx,i,free_pot_idxs[pot]) / N
                inner_total = prx.mass_jac[el_idx,statevar_offset+i]
                for C in 1:num_elements
                    inner_total -= x_A * prx.mass_jac[C,statevar_offset+i]
                end
                total += factor * inner_total
            end
        end
        row[col_offset+pot] = total
    end
    # Δℵ columns
    col_offset += length(free_pot_idxs)
    for β in 1:length(free_phase_idxs)
        prx = phase_records[β]
        factor = 1 / N
        inner_total = prx.mass[el_idx]
        for C in 1:num_elements
            inner_total -= x_A * prx.mass[C]
        end
        row[col_offset+β] = factor * inner_total
    end

    # construct the right-hand-side term
    rhs = 0.0
    # c_iG terms
    for α in 1:length(phase_records)
        prx = phase_records[α]
        statevar_offset = length(prx.state_variables)
        for i in 1:length(prx.state_variables)
            factor = prx.ℵ * c_iG(prx, i) / N
            inner_total = prx.mass_jac[el_idx,statevar_offset+i]
            for C in 1:num_elements
                inner_total -= x_A * prx.mass_jac[C,statevar_offset+i]
            end
            rhs -= factor * inner_total
        end
    end
    # Fixed chemical potentials
    for B in 1:length(fixed_chempot_symbols)
        for α in 1:length(phase_records)
            prx = phase_records[α]
            statevar_offset = length(prx.state_variables)
            for i in 1:length(prx.state_variables)
                factor = prx.ℵ * c_iA(prx, i, B) / N
                inner_total = prx.mass_jac[el_idx,statevar_offset+i]
                for C in 1:num_elements
                    inner_total -= x_A * prx.mass_jac[C,statevar_offset+i]
                end
                rhs -= factor * inner_total
            end
        end
    end
    # Prescribed amount
    rhs += x_A_prescribed
    rhs -= x_A

    return (row, rhs)
end

@doc raw"""
    get_N_row_rhs

See the function `get_N_A_row_rhs`. This is conceptually the same as an N_A
condition with an inner loop over all the elements in each column.

# Equation
```math
\sum_{B_{\mathrm{free}}} \sum_A \sum_\alpha \sum_i \aleph^\alpha \frac{\partial M_A^\alpha}{\partial y_i^\alpha}  c_{iB} \mu_B
+ \sum_\mathrm{Pot} \sum_A \sum_\alpha \sum_i \aleph^\alpha \frac{\partial M_A^\alpha}{\partial y_i^\alpha} c_{i\mathrm{Pot}} \Delta \mathrm{Pot}
+ \sum_\beta \sum_A M_A^\beta \Delta \aleph^\beta  \\
= \sum_A \left(- \sum_\alpha \sum_i \aleph^\alpha \frac{\partial M_A^\alpha}{\partial y_i^\alpha} c_{iG} - \sum_{B_{\mathrm{fixed}}} \sum_\alpha \sum_i \aleph^\alpha \frac{\partial M_A^\alpha}{\partial y_i^\alpha}  c_{iB} \mu_B + \right) + (\tilde{N} - N)
```

"""
function get_N_row_rhs(phase_records, N_prescribed,
                       fixed_chempot_symbols, free_chempot_idxs,
                       fixed_pot_symbols, free_pot_idxs,
                       fixed_phase_symbols, free_phase_idxs,
                       )
    # Construct the row in the equilibrium matrix
    soln_size = (length(free_chempot_idxs) + length(free_pot_idxs) + length(free_phase_idxs))
    num_elements = length(phase_records[1].mass)  # assume mass length is the same for every phase record
    row = Array{Num}(undef, soln_size)
    # Chemical potential columns
    for B in 1:length(free_chempot_idxs)
        total = 0.0
        for A in 1:num_elements
            for α in 1:length(phase_records)
                prx = phase_records[α]
                statevar_offset = length(prx.state_variables)
                for i in 1:length(prx.site_fractions)
                    total += prx.ℵ * prx.mass_jac[A,statevar_offset+i] * c_iA(prx,i,B)
                end
            end
        end
        row[B] = total
    end
    # ΔPotential columns
    col_offset = length(free_chempot_idxs)
    for pot in 1:length(free_pot_idxs)
        total = 0.0
        for A in 1:num_elements
            for α in 1:length(phase_records)
                prx = phase_records[α]
                statevar_offset = length(prx.state_variables)
                for i in 1:length(prx.state_variables)
                    total += prx.ℵ * prx.mass_jac[A,statevar_offset+i] * c_iPot(prx,i,free_pot_idxs[pot])
                end
            end
        end
        row[col_offset+pot] = total
    end
    # Δℵ columns
    col_offset += length(free_pot_idxs)
    for β in 1:length(free_phase_idxs)
        total = 0.0
        for A in 1:num_elements
            total += phase_records[β].mass[A]
        end
        row[col_offset+β] = total
    end

    # construct the right-hand-side term
    rhs = 0.0
    for A in 1:num_elements
        # C_iG terms
        for α in 1:length(phase_records)
            prx = phase_records[α]
            statevar_offset = length(prx.state_variables)
            for i in 1:length(prx.state_variables)
                rhs -= prx.ℵ * prx.mass_jac[A,statevar_offset+i] * c_iG(prx, i)
            end
        end
        # Fixed chemical potential terms
        for B in 1:length(fixed_chempot_symbols)
            for α in 1:length(phase_records)
                prx = phase_records[α]
                statevar_offset = length(prx.state_variables)
                for i in 1:length(prx.state_variables)
                    rhs -= prx.ℵ * prx.mass_jac[A,statevar_offset+i] * fixed_chempot_symbols[B]
                end
            end
        end
    end
    # Prescribed N
    rhs += N_prescribed
    # Current system total N
    for A in 1:num_elements
        for α in 1:length(phase_records)
            prx = phase_records[α]
            rhs -= prx.ℵ * prx.mass[A]
        end
    end

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
    elseif startswith(str_cond, "X_")
        el = str_cond[3:end]
        el_idx = findfirst(x -> x == el, elements)
        row, rhs = get_x_A_row_rhs(phase_records, el_idx, cond, fixed_free_terms...)
    else
        throw("Condition for $str_cond is not yet implemented")
    end

    return row, rhs
end