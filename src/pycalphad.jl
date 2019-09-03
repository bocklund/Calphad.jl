# Interface to pycalphad via PyCall

using PyCall, UUIDs
sympy = pyimport("sympy")
pycalphad = pyimport("pycalphad")
pycalphad.Database
pycalphad.Model

# Regex for binary operations that may not be space-separated
# (bug in SymPy Julia code printing)
fixcodegen(code::String) = replace(code, r"(\.[\*\/\^])" => s" \1 ")

function getfunc(sympyexpr, name::String, func_args::Array{String,1})
    code_str = fixcodegen(sympy.julia_code(sympyexpr))
    func_str =  "function " * name * "(" * join(func_args, ", ") * ")\n" * code_str * "\nend #function"
    return eval(Meta.parse(func_str)), func_args
end # function

"""
    getfunc(model, attr; additional_statevars=Array{String}(undef, 0))

Get a callable function for a Model attribute, with any additional state variables added.
The function is appended with a UUID to ensure the function signature is unique.
"""
function getfunc(model, attr; additional_statevars=Array{String}(undef, 0))
    expr = model.__getattribute__(attr)
    statevars = sort!(vcat([sym.name for sym in model.state_variables], additional_statevars))
    sitefracs = sort!([sym.name for sym in model.site_fractions])
    funcargs = vcat(statevars, sitefracs)
    funcname = attr * "_" * model.phase_name * replace(string(uuid1()), "-" => "_")
    f, funcargs = getfunc(expr, funcname, funcargs)
    return f, funcargs
end # function

"""
    collect_constituents(constituents)

Convert a pycalphad PyObject constituent array to Julia
"""
collect_constituents(constituents) = [[sp.name for sp in subl] for subl in constituents]

"""
    phase_constituents(dbf, phase_name)

Return a constituent array for a phase in a Database.
"""
phase_constituents(dbf, phase_name) = collect_constituents(dbf.phases[phase_name].constituents)

"""
    get_site_ratios(dbf, phase_name)

Return an array of site ratios for a phase in a Database
"""
get_site_ratios(dbf, phase_name) = float.(collect(dbf.phases[phase_name].sublattices))

"""
    get_site_ratios(model)

Return an array of site ratios for a Model instance
"""
get_site_ratios(model) = float.(collect(model.site_ratios))



"""
    active_constituents(model)

Return a constituent array for a model.
Will only contain active species by definition.
"""
active_constituents(model) = collect_constituents(model.constituents)
