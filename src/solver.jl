#=
This module exists to set up the solver as a JuMP Model.

The key method is `build_model` which constructs the model.

=#

using JuMP
import Ipopt


function phase_fracion_variable(model::Model, prx::AbstractPhaseRecord)
	return variable_by_name(model, string("NP", phaserecord(prx).name))
end

function internal_dof_variable(model::Model, prx::AbstractPhaseRecord)
    phase_record = phaserecord(prx)
	const_array = phase_record.constituent_array
	num_subl = length(const_array)

	internal_dof_vec_size = 0
	for subl_idx in 1:num_subl
		internal_dof_vec_size += length(const_array[subl_idx])
	end
	return [variable_by_name(model, string(phase_record.name, "_INTERNAL_DOF[", i, "]")) for i in 1:internal_dof_vec_size]
end


function add_phase_internal_constraints(model::Model, prx::AbstractPhaseRecord)
	phase_record = phaserecord(prx)
	const_array = phase_record.constituent_array
	num_subl = length(const_array)

	internal_dof_vec_size = 0
	for subl_idx in 1:num_subl
		internal_dof_vec_size += length(const_array[subl_idx])
	end
	# internal_dof = @variable(model, [1:internal_dof_vec_size], lower_bound=MIN_SITE_FRACTION, upper_bound=1.0, base_name=string(phase_record.name, "_INTERNAL_DOF"))
	internal_dof = @variable(model, [1:internal_dof_vec_size], lower_bound=MIN_SITE_FRACTION, upper_bound=1.0, base_name=string(phase_record.name, "_INTERNAL_DOF"))

	# add constraints for each sublattice
	dof_offset = 0
	for subl_idx in 1:num_subl
		subl_size = length(const_array[subl_idx])
		@constraint(model, sum(internal_dof[1+dof_offset:subl_size+dof_offset]) == 1)
		dof_offset += subl_size
	end
end

# requires that phase internal constraints are already specified
function moles_from_phase(model::Model, species::String, prx::AbstractPhaseRecord)
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
				vv = variable_by_name(model, string(phase_name, "_INTERNAL_DOF[", internal_dof_idx, "]"))
				mm = @NLexpression(model, vv*phase_record.subl_site_ratios[subl_idx])
				push!(moles_vars, mm)
			end
		end
	end

	if length(moles_vars) == 1
		return moles_vars[1]
	else
		return MOLES = @NLexpression(model, sum(moles_vars[i] for i in 1:length(moles_vars)))
	end
end

function add_phase_fraction_variable(model::Model, prx::AbstractPhaseRecord)::VariableRef
	return @variable(model, lower_bound=0.0, upper_bound=1.0, base_name=string("NP", phaserecord(prx).name))
end

function moles_array(model::Model, phase_records::Array{TY}, comps::Array{String,1}) where {TY<:AbstractPhaseRecord}
	comps = sort(comps)
	num_comps = length(comps)

	# build arary of mole variables
	moles = Array{Union{NonlinearExpression}, 1}(undef, num_comps)
	for comp_idx in 1:num_comps
		prx_moles = []
		for phase_rec in phase_records
			# NP*N(comp)
			mm = moles_from_phase(model, comps[comp_idx], phase_rec)
			npvv = phase_fracion_variable(model, phase_rec)
			push!(prx_moles, @NLexpression(model, mm*npvv))
		end
		moles[comp_idx] = @NLexpression(model, sum(prx_moles[i] for i in 1:length(prx_moles)))
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
	comps = sort(comps)
	num_comps = length(comps)

	model = Model(with_optimizer(Ipopt.Optimizer, print_level=0))

	# State variables
	@variable(model, T >= 0)
	@variable(model, P >= 0)

	phase_fractions = Array{VariableRef, 1}(undef, length(all_phase_records))
	objectives = Array{Any, 1}(undef, length(all_phase_records))
	current_expression = nothing
	for prx_idx in 1:length(all_phase_records)
		phase_rec = all_phase_records[prx_idx]
		add_phase_internal_constraints(model, phase_rec)
		internal_dof = internal_dof_variable(model, phase_rec)
		NP =  add_phase_fraction_variable(model, phase_rec)
		phase_fractions[prx_idx] = NP
		obj_sym = Symbol(string(phase_rec.name, "_OBJ"))
		register(model, obj_sym, 2+length(internal_dof), phase_rec.obj, autodiff=true)
		obj_expr = Expr(:call, :*, NP, Expr(:call, obj_sym, P, T, internal_dof...))
		current_expression = add_expressions(current_expression, obj_expr)
	end

	# phase fractions constraint
	@constraint(model, sum(phase_fractions) == 1.0)
	# add moles constraint, N=1
	moles = moles_array(model, all_phase_records, comps)
	@NLconstraint(model, sum(moles[i] for i in 1:length(moles)) == 1.0)

	set_NL_objective(model, MOI.MIN_SENSE, current_expression)

	return model
end

function add_composition_constraint(model::Model, phase_records::Array{TY, 1}, species::String, constraint_value::Float64, comps::Array{String, 1}) where {TY<:AbstractPhaseRecord}
	comps = sort(comps)
	moles = moles_array(model, phase_records, comps)
	species_idx = findfirst(isequal(species), comps)
	mm = [@NLexpression(model, m) for m in moles]
	@NLconstraint(model, mm[species_idx]/sum(mm[i] for i in 1:length(mm)) == constraint_value)
end
