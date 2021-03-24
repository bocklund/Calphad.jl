using Symbolics, Calphad, Test

let # Single phase test
elements = ["A", "B"];

# Define Gibbs energy and phase record for BETA
@variables Y_BETA_A Y_BETA_B
G_BETA_A = 8000.0-10.0*T;
G_BETA_B = 12000.0-10.0*T;
G_BETA = Y_BETA_A*G_BETA_A + Y_BETA_B*G_BETA_B + R*T*(Y_BETA_A*log(Y_BETA_A) + Y_BETA_B*log(Y_BETA_B));
mass_BETA = [Y_BETA_A, Y_BETA_B];
prx_BETA = PhaseRecord("BETA", G_BETA, mass_BETA, [P, T], [Y_BETA_A, Y_BETA_B], [2]);
phase_records = [prx_BETA];

# Condition variables
@variables X_B  N
condition_dict = Dict(P => 101325.0, T=>300.0, X_B=>0.5, N=>1.0)
cond_keys = [P, T, X_B, N]

# Get the symbolic solution and Δy matrix
sym_soln = Calphad.get_solution(phase_records, elements, cond_keys);
sym_Delta_y_mats = Calphad.get_Delta_y_mat.(phase_records, (elements,), (cond_keys,));

# Construct a composition set close to the solution
compset_BETA = Calphad.CompSet(prx_BETA, [0.51, 0.49], 1.0)
compsets = [compset_BETA];

# Solve with printing
println(condition_dict)
println(compsets)
Calphad.solve_and_update(compsets, condition_dict, sym_soln, sym_Delta_y_mats, length(compsets); step_size=0.01, doprint=true)
println(compsets)

end


let  # Multi-phase test
elements = ["A", "B"];


# Define phase records for ALPHA and BETA
# ALPHA
@variables Y_ALPHA_A Y_ALPHA_B
G_ALPHA = R*T*(Y_ALPHA_A*log(Y_ALPHA_A) + Y_ALPHA_B*log(Y_ALPHA_B));
mass_ALPHA = [Y_ALPHA_A, Y_ALPHA_B];
prx_ALPHA = PhaseRecord("ALPHA", G_ALPHA, mass_ALPHA, [P, T], [Y_ALPHA_A, Y_ALPHA_B], [2]);

# BETA
@variables Y_BETA_A Y_BETA_B
G_BETA_A = 8000.0-10.0*T;
G_BETA_B = 12000.0-10.0*T;
G_BETA = Y_BETA_A*G_BETA_A + Y_BETA_B*G_BETA_B + R*T*(Y_BETA_A*log(Y_BETA_A) + Y_BETA_B*log(Y_BETA_B));
mass_BETA = [Y_BETA_A, Y_BETA_B];
prx_BETA = PhaseRecord("BETA", G_BETA, mass_BETA, [P, T], [Y_BETA_A, Y_BETA_B], [2]);

phase_records = [prx_BETA, prx_ALPHA];

# Conditions
@variables X_B  N
condition_dict = Dict(P => 101325.0, T=>1000.0, X_B=>0.5, N=>1.0)
cond_keys = [P, T, X_B, N]

# Get the symbolic solution and Δy matrix
sym_soln = Calphad.get_solution(phase_records, elements, cond_keys);
sym_Delta_y_mats = Calphad.get_Delta_y_mat.(phase_records, (elements,), (cond_keys,));

# Composition sets
compset_BETA = Calphad.CompSet(prx_BETA, [0.55911824, 0.44088176], 0.502617267178723)
compset_ALPHA = Calphad.CompSet(prx_ALPHA, [0.44025959, 0.55974041], 0.49738273282127704)
compsets = [compset_BETA, compset_ALPHA]

# Solve with printing
println(condition_dict)
println(compsets)
Calphad.solve_and_update(compsets, condition_dict, sym_soln, sym_Delta_y_mats, length(compsets); step_size=0.01, doprint=true)
println(compsets)


end