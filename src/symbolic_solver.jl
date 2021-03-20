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
    str = "ℵ_$(phase_record.phase_name)"
    return Num(Variable(Symbol(str)))
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
compset = CompSet(prx, [0.5, 0.5], 1.0);

r, rrhs = get_N_A_row_rhs([compset], 1, N_A, [], [1, 2], [P, T], [], [], [1])
```
"""
function get_N_A_row_rhs(compsets, el_idx, N_el_sym,
                         fixed_chempot_symbols, free_chempot_idxs,
                         fixed_pot_symbols, free_pot_idxs,
                         fixed_phase_symbols, free_phase_idxs,
                         )

    #compsets = [compset]
    #el_idx = 1
    #N_el_sym = N_A  # symbol for the prescribed amount that will be substituted
    #fixed_chempot_symbols = []
    #free_chempot_idxs = [1, 2]
    #fixed_pot_symbols = [1, 2]
    #free_pot_idxs = []
    #fixed_phase_symbols = []
    #free_phase_idxs = [1]

    # Construct the row in the equilibrium matrix
    soln_size = (length(free_chempot_idxs) + length(free_pot_idxs) + length(free_phase_idxs))
    row = Array{Num}(undef, soln_size)
    # Chemical potential columns
    for B in 1:length(free_chempot_idxs)
        total = 0.0
        for α in 1:length(compsets)
            prx = compsets[α].phase_rec
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
        for α in 1:length(compsets)
            prx = compsets[α].phase_rec
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
        for α in 1:length(compsets)
            prx = compsets[α].phase_rec
            total += prx.mass[el_idx]
        end
        row[col_offset+β] = total
    end

    # construct the right-hand-side term
    rhs = 0.0
    for α in 1:length(compsets)
        prx = compsets[α].phase_rec
        statevar_offset = length(prx.state_variables)
        for i in 1:length(prx.state_variables)
            rhs -= _ℵ(prx) * prx.mass_jac[el_idx,statevar_offset+i] * c_iG(prx, i)
        end
    end
    for B in 1:length(fixed_chempot_symbols)
        for α in 1:length(compsets)
            prx = compsets[α].phase_rec
            statevar_offset = length(prx.state_variables)
            for i in 1:length(prx.state_variables)
                rhs -= _ℵ(prx) * prx.mass_jac[el_idx,statevar_offset+i] * fixed_chempot_symbols[B]
            end
        end
    end
    for α in 1:length(compsets)
        prx = compsets[α].phase_rec
        rhs += prx.mass[el_idx]
    end
    rhs -= N_el_sym

    return (row, rhs)
end


# LaTeX equation:
# \sum_{B_{\mathrm{free}}} M_B^\alpha \mu_B + \sum_\mathrm{Pot} -\frac{\partial G_M^\alpha}{\partial \mathrm{Pot}} \Delta \mathrm{Pot} + \sum_\beta 0 = G_M^\alpha + \sum_{B_{\mathrm{fixed}}} -M_B^\alpha \mu_B
"""

# Examples
```
srow, srhs = get_stable_phase_row_rhs([compset], 1, [], [1, 2], [P, T], [], [], [1])
```
"""
function get_stable_phase_row_rhs(compsets, phase_idx,
                                  fixed_chempot_symbols, free_chempot_idxs,
                                  fixed_pot_symbols, free_pot_idxs,
                                  fixed_phase_symbols, free_phase_idxs,
                                  )
    phase_rec = compsets[phase_idx].phase_rec
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


# Condition on total number of moles, N
# Equation in LaTeX:
# \sum_A \sum_\alpha \left[ \aleph^\alpha \sum_i \frac{\partial M_A^\alpha}{\partial y_i^\alpha} \left( \sum_B c_{iB} \mu_B + c_{iT} \Delta T + c_{iP} \Delta P \right) + M_A^\alpha \Delta \aleph^\alpha \right] = \sum_A \sum_\alpha - c_{iG} + \left( N - \tilde{N} \right)


"""

# Examples
```
using Symbolics
@variables N P T N_A MU_B NP_ALPHA N
elements = ["A", "B", "C"]
phases = ["ALPHA", "BETA", "GAMMA"]
condition_dict = Dict(
   N => 1.0,
   P => 101325.0,
   T => 300.0,
   N_A => 0.5,
   MU_B => -10000,
   NP_ALPHA => 0.8,
)
unpack_indices(elements, phases, condition_dict)
```
"""
function unpack_indices(elements, phases, conditions)
    elements = sort(elements)
    phases = sort(phases)
    POTENTIALS = ("P", "T",)

    fixed_chempot_symbols = []
    fixed_chempot_elements = []
    fixed_pot_symbols = []
    fixed_pot_strings = []
    fixed_phase_symbols = []
    fixed_phase_names = []
    condition_row_symbols = []
    # TODO: move validation outside this function
    for cond in sort(collect(keys(conditions)); by=string)
        str_cond = string(cond)
        if startswith(str_cond, "MU_")
            el = str_cond[4:end]
            if el ∉ elements
                throw("Element $el in condition $str_cond is not in the elements $elements")
            end
            push!(fixed_chempot_symbols, cond)
            push!(fixed_chempot_elements, el)
        elseif str_cond in POTENTIALS
            push!(fixed_pot_symbols, cond)
            push!(fixed_pot_strings, str_cond)
        elseif startswith(str_cond, "NP_")
            phase_name = str_cond[4:end]
            if phase_name ∉ phases
                throw("Phase $phase_name in condition $str_cond is not in the phases $phases")
            end
            push!(fixed_phase_symbols, cond)
            push!(fixed_phase_names, phase_name)
        elseif startswith(str_cond, "N_") | startswith(str_cond, "X_") | (str_cond in ("N", ))
            # TODO: validate N_ and X_ conditions that they contain elements
            push!(condition_row_symbols, cond)
        else
            throw("Unknown condition $cond")
        end
    end

    # The fixed conditions are all in sorted order. Now we have to determine the
    # indices of the free conditions
    free_chempot_idxs = [A for A in 1:length(elements) if elements[A] ∉ fixed_chempot_elements]
    free_pot_idxs = [pot for pot in 1:length(POTENTIALS) if POTENTIALS[pot] ∉ fixed_pot_strings]
    free_phase_idxs = [α for α in 1:length(phases) if phases[α] ∉ fixed_phase_names]

    return (
        condition_row_symbols,
        fixed_chempot_symbols, free_chempot_idxs,
        fixed_pot_symbols, free_pot_idxs,
        fixed_phase_symbols, free_phase_idxs,
    )

end

function cond_row_rhs(cond, val, elements, phases, compsets, fixed_free_terms)
    # assumes elements, phases, compsets are sorted
    str_cond = string(cond)
    if str_cond == "N"
        throw("Condition for $str_cond is not yet implemented")
    elseif startswith(str_cond, "N_")
        el = str_cond[3:end]
        el_idx = findfirst(x -> x == el, elements)
        row, rhs = get_N_A_row_rhs(compsets, el_idx, cond, fixed_free_terms...)
    else
        throw("Condition for $str_cond is not yet implemented")
    end

    return row, rhs
end


function get_solution(compsets, elements, conditions)
    # TODO: in principle, this should be able to take phase records, but I need
    # to be able to symbolically get the phase amount (ℵ) for a phase record
    # symbolically.
    elements = sort(elements)
    compsets = sort(compsets; by = x -> x.phase_rec.phase_name)
    phases = [cs.phase_rec.phase_name for cs in compsets]
    idxs = unpack_indices(elements, phases, conditions)
    condition_row_symbols = first(idxs)
    fixed_free_soln_terms = Base.tail(idxs)

    # construct rows and right-hand-side
    free_chempots = fixed_free_soln_terms[2]
    free_pots = fixed_free_soln_terms[4]
    free_phases = fixed_free_soln_terms[6]

    neqns = length(compsets) + length(condition_row_symbols)
    soln_size = length(free_chempots) + length(free_pots) + length(free_phases)
    @assert neqns == soln_size "The number of matrix equations ($neqns) does not match the number of solution terms ($soln_size)"

    A = Matrix{Num}(undef, neqns, soln_size)
    b = Array{Num}(undef, soln_size)

    # Get the equations
    # One row for each stable phase:
    for phase_idx in 1:length(compsets)
        row, rhs = get_stable_phase_row_rhs(compsets, phase_idx, fixed_free_soln_terms...)
        A[phase_idx, :] = row
        b[phase_idx] = rhs
    end
    # One row for each condition equation
    row_offset = length(compsets)
    for cond_idx in 1:length(condition_row_symbols)
        cond = condition_row_symbols[cond_idx]
        row, rhs = cond_row_rhs(cond, conditions[cond], elements, phases, compsets, fixed_free_soln_terms)
        A[row_offset+cond_idx, :] = row
        b[row_offset+cond_idx] = rhs
    end

    soln = A \ b
    return soln
end


