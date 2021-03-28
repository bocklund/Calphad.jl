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
function unpack_indices(elements, phases, conditions_keys; validate_phases=false)
    elements = sort(elements)
    POTENTIALS = ("P", "T",)

    fixed_chempot_symbols = []
    fixed_chempot_elements = []
    fixed_pot_symbols = []
    fixed_pot_strings = []
    fixed_phase_symbols = []
    fixed_phase_names = []
    condition_row_symbols = []
    # TODO: move validation outside this function
    for cond in conditions_keys
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
        elseif startswith(str_cond, "ℵ_")
            phase_name = str_cond[5:end]  # TODO: not sure why it needs to be 5
            # handle case if unpack_indices is used for a single phase Δy
            if validate_phases
                if phase_name ∉ phases
                    throw("Phase $phase_name in condition $str_cond is not in the phases $phases")
                end
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
    idxs = unpack_indices(elements, [phase_record.phase_name], conditions_keys; validate_phases=false)
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

function isphasecond(cond::Num)
    str_cond = string(cond)
    return startswith(str_cond, "ℵ_")
end

# TODO: test this with vectorize_values. It's important that the uniqueness of
# inputs and the order is preserved between inputs and values.
function vectorize_inputs(phase_records::Vector{PhaseRecord}, free_potentials::OrderedDict{Num,Float64}, conditions::OrderedDict{Num,Float64})
    phase_rec_conds = [vcat(prx.site_fractions, [prx.ℵ]) for prx in phase_records]
    return vcat(collect(keys(free_potentials)), [ky for ky in keys(conditions) if !isphasecond(ky)], phase_rec_conds...)
end

function vectorize_values(compsets::Vector{CompSet}, free_potentials::OrderedDict{Num,Float64}, conditions::OrderedDict{Num,Float64})
    compset_conds = [vcat(cs.Y, [cs.ℵ]) for cs in compsets]
    return vcat(collect(values(free_potentials)), [val for (ky, val) in conditions if !isphasecond(ky)], compset_conds...)
end

function build_callable(expr, inputs)
    # The function at index 1 is a callable with a return, the second one is
    # for in-place modification.
    # `expression=Val{false}` compiles the function instead of returning it,
    # avoiding "world age" issues from using eval.
    f = build_function(expr, inputs, expression=Val{false})[1]
    return f
end

function build_delta_y_callables(delta_y_matrices::Vector{Matrix{Num}}, soln::Vector{Num}, inp::Vector{Num}, num_free_phases::Int)
    soln_potentials = vcat([1], soln[1:end-num_free_phases])
    funcs = Vector{Function}(undef, length(delta_y_matrices))
    for i in 1:length(delta_y_matrices)
        delta_y_expr = delta_y_matrices[i] * soln_potentials
        # The function at index 1 is a callable with a return, the second one is for in-place modification
        funcs[i] = build_callable(delta_y_expr, inp)
    end
    return funcs
end

function update(compsets::Vector{CompSet}, x::Vector{Float64}, soln::Vector{Float64}, delta_y_funcs::Vector{Function}, free_phase_idxs::Vector{Int}; step_size=1.0, verbose=false)
    num_free_phases = length(free_phase_idxs)
    # Update Δℵ
    for β in 1:num_free_phases
        α = free_phase_idxs[β]
        Δℵ = soln[end-num_free_phases+β]
        if verbose
            println("$(compsets[α].phase_record.phase_name): Δℵ=$(Δℵ)")
        end
        # Prevent the phase amount from going negative
        compsets[α].ℵ = max(compsets[α].ℵ+Δℵ, 0.0)
    end
    # Update Δy
    for α in 1:length(compsets)
        Δy = delta_y_funcs[α](x)
        if verbose
            println("$(compsets[α].phase_record.phase_name): Δy=$(Δy) (step size=$step_size)")
        end
        compsets[α].Y += step_size*Δy
    end
end

# TODO: benchmark. Maybe investigate caching or other speed workarounds, depending on how long it takes.
"""
    `solution_functions`

This function does not depend on any concrete values and should be valid for
_any_ specific problem so long as the elements, phases, `keys(free_potentials)`
and `keys(conditions)` are the same.

"""
function solution_functions(elements, phase_records, free_potentials, conditions, free_phase_idxs)
    # Solve symbolically and compile solution function
    cond_keys = collect(keys(conditions))
    A, b = get_solution_parts(phase_records, elements, cond_keys);
    sym_soln = A \ b
    sym_Delta_y_mats = get_Delta_y_mat.(phase_records, (elements,), (cond_keys,));
    inp = vectorize_inputs(phase_records, free_potentials, conditions)
    soln_func = Calphad.build_callable(sym_soln, inp)
    # Compile Δy functions
    delta_y_funcs = Calphad.build_delta_y_callables(sym_Delta_y_mats, sym_soln, inp, length(free_phase_idxs))

    return soln_func, delta_y_funcs
end

# TODO: free_phase_idxs could probably be removed if composition sets can keep
# information about whether or not they are fixed amount.
"""

    `find_solution`

Find a solution to a point calculation.

Composition sets and free potenetials are updated in-place.
Free chemical potentials are returned.
"""
function find_solution(soln_func::Function, delta_y_funcs::Vector{Function},
                       compsets::Vector{CompSet},
                       free_potentials::OrderedDict{Num,Float64},
                       conditions::OrderedDict{Num,Float64},
                       free_phase_idxs::Vector{Int};
                       max_iters=1000, verbose=false)
    num_free_potentials = length(free_potentials)

    # Do an initial solve to get the dimension of the solution.
    # This call should be fast and we don't care about the data
    x = vectorize_values(compsets, free_potentials, conditions)
    soln = soln_func(x)
    num_free_chempots = length(soln) - num_free_potentials - length(free_phase_idxs)
    chempots = Vector{Float64}(undef, num_free_chempots)

    for iter in 1:max_iters
        # Step 1: solve
        x = vectorize_values(compsets, free_potentials, conditions)
        soln = soln_func(x)
        if verbose
            println("x: ", x)
            println("equilibrium_soln: ", soln)
        end

        # Check convergence of chemical potentials
        converged = true
        if maximum(abs.(chempots - soln[1:num_free_chempots])) > 1e-12
            converged = false
        end

        # Extract chemical potentials
        chempots[:] = soln[1:num_free_chempots]

        # Check convergence of potentials
        if num_free_potentials > 0
            offset = num_free_chempots
            if maximum(abs.(soln[offset+1:offset+num_free_potentials])) > 1e-12
                converged = false
            end
        end

        # Check convergence of phase amounts
        if length(free_phase_idxs) > 0
            offset = num_free_chempots + num_free_potentials
            if maximum(abs.(soln[offset+1:end])) > 1e-12
                converged = false
            end
        end

        # TODO: add convergence check for Δy without updating?
        # Updating means the chemical potentials need to be solved for again.
        # However, the Δy are contributed to by the change in the solution
        # quantities, so converging them may be enough.
        if converged
            return chempots, iter
        else
            # Step 2: update
            # Update free potentials
            # TODO: assumes free potentials are in the same order as the
            # potentials in the solution vector, is there a way to enforce this
            # structurally without incurring overhead time? Maybe nothing
            # needs to be done if unpack_indices isn't sorting them
            for (i, (ky, val)) in enumerate(free_potentials)
                # offset for number of free chemical potentials
                free_potentials[ky] = val + soln[num_free_chempots+i]
            end
            # Update composition sets
            # TODO: step size control
            update(compsets, x, soln, delta_y_funcs, free_phase_idxs; step_size=0.1, verbose=false)

            # TODO: add/remove composition sets
        end
    end
    return chempots, max_iters
end