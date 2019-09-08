#=
This module exists to set up the solver as a JuMP Model.

The key method is `build_model` which constructs the model.

=#

import JuMP
import Ipopt


function phase_fracion_variable(model::JuMP.Model, identifier)
	return JuMP.variable_by_name(model, "NP$identifier")
end

function internal_dof_variable(model::JuMP.Model, prx::AbstractPhaseRecord, identifier)
    phase_record = phaserecord(prx)
	const_array = phase_record.constituent_array
	num_subl = length(const_array)

	internal_dof_vec_size = 0
	for subl_idx in 1:num_subl
		internal_dof_vec_size += length(const_array[subl_idx])
	end
	return [JuMP.variable_by_name(model, "INTERNAL_DOF$(identifier)[$(i)]") for i in 1:internal_dof_vec_size]
end


function add_phase_internal_constraints(model::JuMP.Model, prx::AbstractPhaseRecord, identifier)
	phase_record = phaserecord(prx)
	const_array = phase_record.constituent_array
	num_subl = length(const_array)

	internal_dof_vec_size = 0
	for subl_idx in 1:num_subl
		internal_dof_vec_size += length(const_array[subl_idx])
	end
	internal_dof = JuMP.@variable(model, [1:internal_dof_vec_size], lower_bound=MIN_SITE_FRACTION, upper_bound=1.0, base_name="INTERNAL_DOF$(identifier)")

	# add constraints for each sublattice
	dof_offset = 0
	for subl_idx in 1:num_subl
		subl_size = length(const_array[subl_idx])
		JuMP.@constraint(model, sum(internal_dof[1+dof_offset:subl_size+dof_offset]) == 1)
		dof_offset += subl_size
	end
end

# requires that phase internal constraints are already specified
function moles_from_phase(model::JuMP.Model, species::String, prx::AbstractPhaseRecord, identifier)
	phase_record = phaserecord(prx)
	const_array = phase_record.constituent_array
	phase_name = phase_record.name
	num_subl = length(const_array)

	moles_vars = []
	internal_dof_idx = 0
	for subl_idx in 1:num_subl
		for const_idx in 1:length(const_array[subl_idx])
			internal_dof_idx += 1
			if const_array[subl_idx][const_idx] == species
				vv = JuMP.variable_by_name(model, "INTERNAL_DOF$(identifier)[$(internal_dof_idx)]")
				mm = JuMP.@NLexpression(model, vv*phase_record.subl_site_ratios[subl_idx])
				push!(moles_vars, mm)
			end
		end
	end

	if length(moles_vars) == 1
		return moles_vars[1]
	else
		return MOLES = JuMP.@NLexpression(model, sum(moles_vars[i] for i in 1:length(moles_vars)))
	end
end

function add_phase_fraction_variable(model::JuMP.Model, identifier)::JuMP.VariableRef
	return JuMP.@variable(model, lower_bound=0.0, upper_bound=1.0, base_name="NP$(identifier)")
end

function moles_array(model::JuMP.Model, phase_records::Array{TY}, comps::Array{String,1}) where {TY<:AbstractPhaseRecord}
	multiphase_problem = length(phase_records) > 1
	comps = sort(comps)
	num_comps = length(comps)

	# build arary of mole variables
	moles = Array{Union{JuMP.NonlinearExpression}, 1}(undef, num_comps)
	for comp_idx in 1:num_comps
		prx_moles = []
		for i in 1:length(phase_records)
			phase_rec = phase_records[i]
			# NP*N(comp)
			mm = moles_from_phase(model, comps[comp_idx], phase_rec, i)
			if multiphase_problem
				npvv = phase_fracion_variable(model, i)
				push!(prx_moles, JuMP.@NLexpression(model, mm*npvv))
			else
				push!(prx_moles, JuMP.@NLexpression(model, mm))
			end # if
		end
		moles[comp_idx] = JuMP.@NLexpression(model, sum(prx_moles[i] for i in 1:length(prx_moles)))
	end
	return moles
end

function add_expressions(current::Union{Expr, Nothing}, new::Expr)::Expr
	if current == nothing
		return new
	else
		return Expr(:call, :+, current, new)
	end
end

function build_model(all_phase_records::Array{TY, 1}, comps::Array{String,1}) where {TY<:AbstractPhaseRecord}
	multiphase_problem = length(all_phase_records) > 1
	comps = sort(comps)
	num_comps = length(comps)

	model = JuMP.Model(JuMP.with_optimizer(Ipopt.Optimizer, print_level=3))

	# State variables
	# TODO: fix implicit P, T conditions
	num_statevars = 2
	JuMP.@variable(model, T >= 0)
	JuMP.@variable(model, P >= 0)

	phase_fractions = Array{JuMP.VariableRef, 1}(undef, length(all_phase_records))
	objectives = Array{Any, 1}(undef, length(all_phase_records))
	current_expression = nothing
	for prx_idx in 1:length(all_phase_records)
		phase_rec = all_phase_records[prx_idx]
		add_phase_internal_constraints(model, phase_rec, prx_idx)
		internal_dof = internal_dof_variable(model, phase_rec, prx_idx)
		obj_sym = Symbol("$(prx_idx)_OBJ")
		if multiphase_problem
			NP =  add_phase_fraction_variable(model, prx_idx)
			phase_fractions[prx_idx] = NP
			obj_expr = Expr(:call, :*, NP, Expr(:call, obj_sym, P, T, internal_dof...)) # TODO: implicit P, T  conditions
		else
			obj_expr = Expr(:call, obj_sym, P, T, internal_dof...) # TODO: implicit P, T  conditions
		end # if
		JuMP.register(model, obj_sym, num_statevars+length(internal_dof), phase_rec.obj, autodiff=true)
		current_expression = add_expressions(current_expression, obj_expr)
	end

	if multiphase_problem
		# phase fractions constraint
		JuMP.@constraint(model, sum(phase_fractions) == 1.0)
	end # if
	# add moles constraint, N=1
	moles = moles_array(model, all_phase_records, comps)
	JuMP.@NLconstraint(model, sum(moles[i] for i in 1:length(moles)) == 1.0)

	JuMP.set_NL_objective(model, JuMP.MOI.MIN_SENSE, current_expression)

	return model
end

function add_composition_constraint(model::JuMP.Model, phase_records::Array{TY, 1}, species::String, constraint_value::Float64, comps::Array{String, 1}) where {TY<:AbstractPhaseRecord}
	comps = sort(comps)
	moles = moles_array(model, phase_records, comps)
	species_idx = findfirst(isequal(species), comps)
	mm = [JuMP.@NLexpression(model, m) for m in moles]
	JuMP.@NLconstraint(model, mm[species_idx]/sum(mm[i] for i in 1:length(mm)) == constraint_value)
end
