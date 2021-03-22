# This module contains the abstractions for the solver based on the Symbolics package

include("_equilibrium_matrix_conditions.jl")

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
function unpack_indices(elements, phases, conditions_keys)
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
    for cond in sort(collect(conditions_keys); by=string)
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

function get_solution_parts(phase_records, elements, conditions_keys)
    conditions_keys = collect(conditions_keys)
    elements = sort(elements)
    phase_records = sort(phase_records; by = x -> x.phase_name)
    phases = [prx.phase_name for prx in phase_records]
    idxs = unpack_indices(elements, phases, conditions_keys)
    condition_row_symbols = first(idxs)
    fixed_free_soln_terms = Base.tail(idxs)

    # construct rows and right-hand-side
    free_chempots = fixed_free_soln_terms[2]
    free_pots = fixed_free_soln_terms[4]
    free_phases = fixed_free_soln_terms[6]

    neqns = length(phase_records) + length(condition_row_symbols)
    soln_size = length(free_chempots) + length(free_pots) + length(free_phases)
    @assert neqns == soln_size "The number of matrix equations ($neqns) does not match the number of solution terms ($soln_size)"

    A = Matrix{Num}(undef, neqns, soln_size)
    b = Array{Num}(undef, soln_size)

    # Get the equations
    # One row for each stable phase:
    for phase_idx in 1:length(phase_records)
        row, rhs = get_stable_phase_row_rhs(phase_records, phase_idx, fixed_free_soln_terms...)
        A[phase_idx, :] = row
        b[phase_idx] = rhs
    end
    # One row for each condition equation
    row_offset = length(phase_records)
    for cond_idx in 1:length(condition_row_symbols)
        cond = condition_row_symbols[cond_idx]
        row, rhs = cond_row_rhs(cond, elements, phases, phase_records, fixed_free_soln_terms)
        A[row_offset+cond_idx, :] = row
        b[row_offset+cond_idx] = rhs
    end

    return A, b
end

function get_solution(phase_records, elements, conditions_keys)
    A, b = get_solution_parts(phase_records, elements, conditions_keys)
    soln = A \ b
    return soln
end


@doc raw"""
    get_Delta_y_mat

Notationally, Δy is a vector of length(phase_record.site_fractions) that
updates the site fractions. However in this case, we need to plug in the
results from the solution vector that we do not have symbolic variables for. To
resolve that, we'll design the site fractions to be a matrix that can be used
to take the dot product with the solution vector (padded with a prefixed one to
handle the c_iG term). The usage is therefore

# Examples
```
delta_y_M = get_Delta_y_mat(prx, ["A", "B"], [T, P, N_A, N_B])
# ... compute solution
soln = [3271.04833019, 7271.04833015, 1e-16]
delta_y = delta_y_M * vcat(1, soln...)
```

# Equation
```math
\sum_i \Delta y_i^\alpha = \sum_i \left(c_{iG}^\alpha + \sum_A c_{iA}^\alpha \mu_A + \sum_\mathrm{Pot} c_{i\mathrm{Pot}}^\alpha \Delta \mathrm{Pot} \right)
```
"""
function get_Delta_y_mat(phase_record, elements, conditions_keys)
    conditions_keys = collect(conditions_keys)
    elements = sort(elements)
    idxs = unpack_indices(elements, [], conditions_keys)
    fixed_free_soln_terms = Base.tail(idxs)
    num_free_chempots = length(fixed_free_soln_terms[2])
    num_free_pots = length(fixed_free_soln_terms[4])

    Δys = Matrix{Num}(undef, length(phase_record.site_fractions), 1+num_free_chempots+num_free_pots)
    for i in 1:length(phase_record.site_fractions)
        Δys[i, 1] = c_iG(phase_record, i)
        offset = 1
        for A in 1:num_free_chempots
            Δys[i, offset+A] = c_iA(phase_record, i, A)
        end
        offset += num_free_chempots
        for pot in 1:num_free_pots
            Δys[i, offset+pot] = c_iPot(phase_record, i, pot)
        end
    end
    return Δys
end

function get_subs_dict(compsets, conditions_dict)
    d = Dict{Num, Num}()
    for compset in compsets
        for (Y_sym, Y_num) in zip(compset.phase_rec.site_fractions, compset.Y)
            d[Y_sym] = Num(Y_num)
        end
        d[_ℵ(compset.phase_rec)] = compset.ℵ
    end
    for cond in keys(conditions_dict)
        d[cond] = conditions_dict[cond]
    end
    return d
end

# TODO: if possible use the callable functions here, as below:
# subs_dict = Calphad.get_subs_dict([compset], condition_dict)
# condkeys = collect(keys(subs_dict))
# convalues = map(x -> x.val, collect(values(subs_dict)))
# f_soln_sub = eval(build_function(sym_soln, condkeys)[1]);
# On a simple system, using @btime on the `substitute` version took 1.4ms to
# subs the subs_dict into the sym_soln, but calling the f_soln_sub function
# btimes to 373 ns. HUGE for performance.
# It should be possible to sort the conditions by
# `sort(collect(conditions), by=x->x[1])`, which can ensure that the arguments
# to a collable function are in the correct order. Thus, we can make an
# intermediate function that creates callables for the sym_soln and sym_Delta_y_mats
# and then create a new solve_and_update function can skip the substitution step
# and just do args=map(last(sorted_subs_dict)). The sorting isn't too much more
# expensive to do in the loop:
# 3.295 μs: @btime subs_dict = get_subs_dict([compset], condition_dict)
# 8.994 μs: @btime subs_dict = sort(collect(get_subs_dict([compset], condition_dict)), by=x->string(x[1]))
# If I can avoid collecting and sort in place, it would probably reduce the
# number of allocations (3x as much) and close the speed gap. There's also
# probably algorithmic improvements to make so that the sorting happens outside
# the tight solve loop (i.e. no creating subs_dict insice) and we just get the
# non-condition terms at runtime (just site fractions and phase amounts).
# TODO: test that phase amount updating works properly. Most I've seen so far
# keep Δℵ close to zero (single phase, but I think phase amount should still)
# change because the mass condition RHS shouldn't be satisfied, that is:
# (`(N_A - Ñ_A) != 0`)
function solve_and_update(compsets, conditions, sym_soln, sym_Delta_y_mats, num_free_phases)
    # assumes compsets and sym_Delta_y_mats are sorted in order
    subs_dict = get_subs_dict(compsets, conditions)
    soln = substitute.(sym_soln, (subs_dict,))
    println("Chemical potentials: $(Symbolics.value.(soln[1:2]))")

    for α in 1:length(compsets)
        Δy = substitute.(sym_Delta_y_mats[α] * vcat(1, sym_soln[1:end-num_free_phases]...), (subs_dict,))
        Δℵ = soln[end-num_free_phases+α]
        println("$(compsets[α].phase_rec.phase_name): Δy:$(Δy) Δℵ:$(Δℵ)")
        compsets[α].Y -= Δy
        compsets[α].ℵ -= Δℵ
    end
end