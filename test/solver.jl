using Symbolics, Calphad, Test

let # Single phase test
# This case demonstrates the need to have the mass residual be negative.
elements = ["A", "B"];

# Define Gibbs energy and phase record for ALPHA
# ALPHA is an ideal solution
@variables Y_ALPHA_A Y_ALPHA_B
G_ALPHA = R*T*(Y_ALPHA_A*log(Y_ALPHA_A) + Y_ALPHA_B*log(Y_ALPHA_B));
mass_ALPHA = [Y_ALPHA_A, Y_ALPHA_B];
prx_ALPHA = PhaseRecord("ALPHA", G_ALPHA, mass_ALPHA, [P, T], [Y_ALPHA_A, Y_ALPHA_B], [2]);
phase_records = [prx_ALPHA];

# Condition variables
@variables X_B  N
condition_dict = Dict(P => 101325.0, T=>300.0, X_B=>0.5, N=>1.0)
cond_keys = [P, T, X_B, N]

# Get the symbolic solution and Δy matrix
sym_soln = Calphad.get_solution(phase_records, elements, cond_keys);
sym_Delta_y_mats = Calphad.get_Delta_y_mat.(phase_records, (elements,), (cond_keys,));

# Construct a composition set close to the solution
compset_ALPHA = Calphad.CompSet(prx_ALPHA, [0.51, 0.49], 1.0)
compsets = [compset_ALPHA];

# Solve with printing
println(condition_dict)
println(compsets)
# For an ideal solution, Newton's method step should be exact, so we take a full size step.
Calphad.solve_and_update(compsets, condition_dict, sym_soln, sym_Delta_y_mats, [1], Dict(); step_size=1.0, doprint=true)
println(compsets)
@test all(compset_ALPHA.Y .≈ [0.5, 0.5])
@test compset_ALPHA.ℵ ≈ 1.0

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
Calphad.solve_and_update(compsets, condition_dict, sym_soln, sym_Delta_y_mats, [1, 2], Dict(); step_size=0.01, doprint=true)
println(compsets)
# pycalphad solution for step_size = 0.01
@test all(compset_BETA.Y .≈ [0.5591255569183224, 0.4408744430816775])
@test compset_BETA.ℵ ≈ 0.4999626227785548
@test all(compset_ALPHA.Y .≈ [0.4402585398297285, 0.5597414601702715])
@test compset_ALPHA.ℵ ≈ 0.5000373772214453

end

let  # Variable T, fixed phase amount test
# this test "works", but ΔT is really poorly scaled and even the sign doesn't quite converge
# to the correct value (but it's close). Maybe a bug in my equations. We'll keep the test
# around just to be sure that the case runs, but it isn't correct yet.
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
@variables X_B  N ℵ_BETA
condition_dict = Dict(P => 101325.0, ℵ_BETA=>0.0, X_B=>0.5, N=>1.0)
cond_keys = [P, ℵ_BETA, X_B, N]

# Get the symbolic solution and Δy matrix
sym_soln = Calphad.get_solution(phase_records, elements, cond_keys);
sym_Delta_y_mats = Calphad.get_Delta_y_mat.(phase_records, (elements,), (cond_keys,));

# Composition sets
compset_BETA = Calphad.CompSet(prx_BETA, [0.55911824, 0.44088176], 0.0)
compset_ALPHA = Calphad.CompSet(prx_ALPHA, [0.5, 0.5], 1.0)
compsets = [compset_BETA, compset_ALPHA]

# Guess at the free potential value
free_pot_guess = Dict(T => (2, 1000.0))

# Solve with printing
println(condition_dict)
println(compsets)
Calphad.solve_and_update(compsets, condition_dict, sym_soln, sym_Delta_y_mats, [2], free_pot_guess; step_size=0.01, doprint=true)
println(compsets)
println(free_pot_guess)

end


let # pseudo-binary chemical potential conditions
# The phase amount and Oxygen site fraction should remain the same.
elements = ["A", "B", "O"];

@variables Y_M2O3_A Y_M2O3_B Y_M2O3_O
G_M2O3 = R*T*(2*Y_M2O3_A*log(Y_M2O3_A) + 2*Y_M2O3_B*log(Y_M2O3_B) + 3*Y_M2O3_O*log(Y_M2O3_O));
mass_M2O3 = [2*Y_M2O3_A, 2*Y_M2O3_B, 3*Y_M2O3_O];
prx_M2O3 = PhaseRecord("M2O3", G_M2O3, mass_M2O3, [P, T], [Y_M2O3_A, Y_M2O3_B, Y_M2O3_O], [2, 1]);

phase_records = [prx_M2O3];

# Conditions
@variables X_B  N MU_O X_O
condition_dict = Dict(P => 101325.0, T=>1000.0, MU_O=>0.0, X_B=>0.2, N=>1.0)
cond_keys = [P, T, MU_O, X_B, N]

# Get the symbolic solution and Δy matrix
sym_soln = Calphad.get_solution(phase_records, elements, cond_keys);
sym_Delta_y_mats = Calphad.get_Delta_y_mat.(phase_records, (elements,), (cond_keys,));

# Composition sets
compset_M2O3 = Calphad.CompSet(prx_M2O3, [0.51, 0.49, 1.0], 0.2)
compsets = [compset_M2O3]

# Solve with printing
println(condition_dict)
println(compsets)
Calphad.solve_and_update(compsets, condition_dict, sym_soln, sym_Delta_y_mats, [1], Dict(); step_size=0.01, doprint=true)
println(compsets)
Calphad.solve_and_update(compsets, condition_dict, sym_soln, sym_Delta_y_mats, [1], Dict(); step_size=0.01, doprint=true)
println(compsets)
@test !any(isnan.(compsets[1].Y))
@test compsets[1].ℵ ≈ 0.2
@test compsets[1].Y[3] ≈ 1.0
end

