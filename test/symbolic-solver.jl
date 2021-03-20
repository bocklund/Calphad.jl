using Symbolics

function get_subs_dict(compsets, conditions_dict)
    d = Dict(conditions_dict...)
    for compset in compsets
        for (Y_sym, Y_num) in zip(compset.phase_rec.site_fractions, compset.Y)
            d[Y_sym] = Y_num
        end
    end
    return d
end

# Super simple A-B system
@variables Y_BETA_A Y_BETA_B

@variables X_A X_B N_A N_B
@variables N_AL N_TI N_O X_AL X_TI X_O

G_BETA_A = 8000.0-10.0*T;
G_BETA_B = 12000.0-10.0*T;
G_BETA = Y_BETA_A*G_BETA_A + Y_BETA_B*G_BETA_B + R*T*(Y_BETA_A*log(Y_BETA_A) + Y_BETA_B*log(Y_BETA_B));
mass_BETA = [Y_BETA_A, Y_BETA_B];
state_variables = [T];
site_fractions = [Y_BETA_A, Y_BETA_B];
prx = PhaseRecord("BETA", G_BETA, mass_BETA, state_variables, site_fractions);

condition_dict = Dict(
N => 1.0,
P => 101325.0,
T => 300.0,
N_A => 0.5,
N_B => 0.5,
)
compset = CompSet(prx, [0.5, 0.5], 1.0);
soln = solve([compset], condition_dict)
# TODO: chemical potentials flipped, maybe a sign error
println("Actual chemical potentials: ", soln[1:2]);
println("Expected chemical potentials: ", [3271.04833019, 7271.04833015]);

condition_dict = Dict(
N => 1.0,
P => 101325.0,
T => 300.0,
N_A => 0.25,
N_B => 0.75,
)
compset = CompSet(prx, [0.25, 0.75], 1.0);
soln = solve([compset], Dict(T=>300.0))
# TODO: chemical potentials flipped, maybe a sign error
println("Actual chemical potentials: ", soln[1:2]);
println("Expected chemical potentials: ", [1542.0966603, 8282.42022259]);



# # Define molar Gibbs energy for (Al,Ti)2(O)3
# @variables YM2O30AL YM2O30TI YM2O31O T

# G_M_AL2O3 = -1772163.19+1053.4548*T-156.058*T*log(T)+.00709105*T^2-6.29402E-07*T^3 +12366650*T^(-1)
# G_M_TI2O3 = -1581243.06+2395423.68*T^(-1)+940.164783*T-147.673862*T*log(T)-.00173711312*T^2-1.53383348E-10*T^3+19508-2.4214*T

# G_M = G_M_AL2O3*YM2O30AL + G_M_TI2O3*YM2O30TI + R*T*(2*YM2O30AL*log(YM2O30AL) + 2*YM2O30TI*log(YM2O30TI) + 3*YM2O31O*log(YM2O31O))
# G_m = G_M/(2*YM2O30AL + 2*YM2O30TI + 3*YM2O31O)
# # Formula
# mass = [2*YM2O30AL, 3*YM2O31O, 2*YM2O30TI] # moles(Al), moles(O), moles(Ti)

# # substitute(G_M, Dict(YM2O30AL => 0.5, YM2O30TI => 0.5, YM2O31O => 1.0, T => 2000.0))
# # substitute(G_M, Dict(YM2O30AL => 1.0, YM2O30TI => 1e-15, YM2O31O => 1.0, T => 2000.0))
# # substitute(G_m, Dict(YM2O30AL => 0.5, YM2O30TI => 0.5, YM2O31O => 1.0, T => 2000.0))

# state_variables = [T]
# site_fractions = [YM2O30AL, YM2O30TI, YM2O31O]

# prx = PhaseRecord("CORUNDUM", G_M, mass, state_variables, site_fractions)
# compset = CompSet(prx, [0.5, 0.5, 1.0], 0.2)

# A = get_equilibrium_matrix([compset])
# b = get_equilibrium_soln([compset])
# x = A \ b
# subs_dict = Dict(Dict(zip(compset.phase_rec.site_fractions, compset.Y))..., T => 2000.0)
# substitute.(x, (subs_dict,))

# # Using lapack least squares, need to convert A and b into floats first

# AA = Symbolics.value.(substitute.(A, (subs_dict,)))
# bb = Symbolics.value.(substitute.(b, (subs_dict,)))
# LinearAlgebra.LAPACK.gelsd!(AA, bb)

# ############################################################

# # PHASE LIQUID % 1 1 !
# #  CONSTITUENT LIQUID :AL,TI,O : !
# #  PARAM G (LIQUID,AL;0) 1 -10000; 10000 N !
# #  PARAM G (LIQUID,O;0) 1 -40000; 10000 N !
# #  PARAM G (LIQUID,TI;0) 1 -20000; 10000 N !

# @variables LIQUID0AL LIQUID0O LIQUID0TI T

# G_M_LIQUID = -10000*LIQUID0AL + -40000*LIQUID0O + R*T*(LIQUID0AL*log(LIQUID0AL) + LIQUID0O*log(LIQUID0O))
# mass = [LIQUID0AL, LIQUID0O]
# state_variables = [T]
# site_fractions = [LIQUID0AL, LIQUID0O]
# prx = PhaseRecord("LIQUID", G_M_LIQUID, mass, state_variables, site_fractions)
# compset = CompSet(prx, [0.5, 0.5], 1.0)

# # Solve
# A = get_equilibrium_matrix([compset]);
# b = get_equilibrium_soln([compset]);
# # Symbolic solution
# x = A \ b;
# subs_dict = Dict(Dict(zip(compset.phase_rec.site_fractions, compset.Y))..., T => 300.0)
# substitute.(x, (subs_dict,))


# ############################################################





# # Numerical linear least squares soln
# AA = Symbolics.value.(substitute.(A, (subs_dict,)));
# bb = Symbolics.value.(substitute.(b, (subs_dict,)));

# LinearAlgebra.LAPACK.gelsd!(AA, bb)
