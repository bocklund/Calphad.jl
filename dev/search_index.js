var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = Calphad","category":"page"},{"location":"#Calphad","page":"Home","title":"Calphad","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for Calphad.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Calphad]","category":"page"},{"location":"#Calphad.c_iA-Tuple{Any, Any, Any}","page":"Home","title":"Calphad.c_iA","text":"c_iA\n\nEquation\n\nc^alpha_iA = sum_j e^alpha_ij fracpartial M^alpha_Apartial y^alpha_j\n\n\n\n\n\n","category":"method"},{"location":"#Calphad.c_iG-Tuple{Any, Any}","page":"Home","title":"Calphad.c_iG","text":"c_iG\n\nEquation\n\nc^alpha_iG = -sum_j e^alpha_ij fracpartial G^alpha_Mpartial y^alpha_j\n\n\n\n\n\n","category":"method"},{"location":"#Calphad.c_iPot-Tuple{Any, Any, Any}","page":"Home","title":"Calphad.c_iPot","text":"c_iPot\n\nEquation\n\nc^alpha_imathrmPot = -sum_j e^alpha_ij fracpartial^2 G^alpha_Mpartial mathrmPot partial y^alpha_j\n\n\n\n\n\n","category":"method"},{"location":"#Calphad.cond_row_rhs-NTuple{5, Any}","page":"Home","title":"Calphad.cond_row_rhs","text":"cond_row_rhs\n\nThis function takes a condition symbol and dispatches on the correct method that condition that returns the elements in the row and the right hand side.\n\nThe number of columns are determined by the fixed_free_terms - there's one column for each free chemical potential, free potential, and free phase amount.\n\n\n\n\n\n","category":"method"},{"location":"#Calphad.get_Delta_y_mat-Tuple{Any, Any, Any}","page":"Home","title":"Calphad.get_Delta_y_mat","text":"get_Delta_y_mat\n\nNotationally, Δy is a vector of length(phaserecord.sitefractions) that updates the site fractions. However in this case, we need to plug in the results from the solution vector that we do not have symbolic variables for. To resolve that, we'll design the site fractions to be a matrix that can be used to take the dot product with the solution vector (padded with a prefixed one to handle the c_iG term). The usage is therefore\n\nExamples\n\ndelta_y_M = get_Delta_y_mat(prx, [\"A\", \"B\"], [T, P, N_A, N_B])\n# ... compute solution\nsoln = [3271.04833019, 7271.04833015, 1e-16]\ndelta_y = delta_y_M * vcat(1, soln...)\n\nEquation\n\nsum_i Delta y_i^alpha = sum_i left(c_iG^alpha + sum_A c_iA^alpha mu_A + sum_mathrmPot c_imathrmPot^alpha Delta mathrmPot right)\n\n\n\n\n\n","category":"method"},{"location":"#Calphad.get_N_A_row_rhs-NTuple{9, Any}","page":"Home","title":"Calphad.get_N_A_row_rhs","text":"get_N_A_row_rhs\n\nExamples\n\nusing Symbolics\n@variables T P N_A\n@variables Y_BETA_A Y_BETA_B\nG_BETA_A = 8000.0-10.0*T;\nG_BETA_B = 12000.0-10.0*T;\nG_BETA = Y_BETA_A*G_BETA_A + Y_BETA_B*G_BETA_B + R*T*(Y_BETA_A*log(Y_BETA_A) + Y_BETA_B*log(Y_BETA_B));\nmass_BETA = [Y_BETA_A, Y_BETA_B];\nstate_variables = [P, T];\nsite_fractions = [Y_BETA_A, Y_BETA_B];\nprx = PhaseRecord(\"BETA\", G_BETA, mass_BETA, state_variables, site_fractions);\n\nr, rrhs = get_N_A_row_rhs([prx], 1, N_A, [], [1, 2], [P, T], [], [], [1])\n\nEquation\n\nsum_B_mathrmfree sum_alpha sum_i aleph^alpha fracpartial M_A^alphapartial y_i^alpha  c_iB mu_B + sum_mathrmPot sum_alpha sum_i aleph^alpha fracpartial M_A^alphapartial y_i^alpha c_imathrmPot Delta mathrmPot + sum_beta M_A^beta Delta aleph^beta  \n= - sum_alpha sum_i aleph^alpha fracpartial M_A^alphapartial y_i^alpha c_iG - sum_B_mathrmfixed sum_alpha sum_i aleph^alpha fracpartial M_A^alphapartial y_i^alpha  c_iB mu_B + left( tildeN_A - sum_alpha aleph^alpha M^alpha_A right)\n\n\n\n\n\n","category":"method"},{"location":"#Calphad.get_N_row_rhs-NTuple{8, Any}","page":"Home","title":"Calphad.get_N_row_rhs","text":"get_N_row_rhs\n\nSee the function get_N_A_row_rhs. This is conceptually the same as an N_A condition with an inner loop over all the elements in each column.\n\nEquation\n\nsum_B_mathrmfree sum_A sum_alpha sum_i aleph^alpha fracpartial M_A^alphapartial y_i^alpha  c_iB mu_B\n+ sum_mathrmPot sum_A sum_alpha sum_i aleph^alpha fracpartial M_A^alphapartial y_i^alpha c_imathrmPot Delta mathrmPot\n+ sum_beta sum_A M_A^beta Delta aleph^beta  \n= sum_A left(- sum_alpha sum_i aleph^alpha fracpartial M_A^alphapartial y_i^alpha c_iG - sum_B_mathrmfixed sum_alpha sum_i aleph^alpha fracpartial M_A^alphapartial y_i^alpha  c_iB mu_B + right) + (tildeN - N)\n\n\n\n\n\n","category":"method"},{"location":"#Calphad.get_stable_phase_row_rhs-NTuple{8, Any}","page":"Home","title":"Calphad.get_stable_phase_row_rhs","text":"get_stable_phase_row_rhs\n\nExamples\n\nsrow, srhs = get_stable_phase_row_rhs([prx], 1, [], [1, 2], [P, T], [], [], [1])\n\nEquation\n\nsum_B_mathrmfree M_B^alpha mu_B + sum_mathrmPot -fracpartial G_M^alphapartial mathrmPot Delta mathrmPot + sum_beta 0 = G_M^alpha + sum_B_mathrmfixed -M_B^alpha mu_B\n\n\n\n\n\n","category":"method"},{"location":"#Calphad.get_x_A_row_rhs-NTuple{9, Any}","page":"Home","title":"Calphad.get_x_A_row_rhs","text":"get_x_A_row_rhs\n\nSee the function get_N_A_row_rhs. This is the corresponding mole fraction condition.\n\nEquation\n\nsum_B_mathrmfree sum_alpha sum_i fracaleph^alpha c_iBN left( fracpartial M_A^alphapartial y_i^alpha - x_A sum_C fracpartial M_C^alphapartial y_i^alpha right)  mu_B  \n+ sum_mathrmPot sum_alpha sum_i fracaleph^alpha c_imathrmPotN left( fracpartial M_A^alphapartial y_i^alpha - x_A sum_C fracpartial M_C^alphapartial y_i^alpha right)  Delta mathrmPot  \n+ sum_beta_mathrmfree frac1N left( M_A^beta - x_A sum_C M_C^beta right) Delta aleph^beta \n= - sum_alpha sum_i fracaleph^alpha c_iGN left( fracpartial M_A^alphapartial y_i^alpha - x_A sum_C fracpartial M_C^alphapartial y_i^alpha right) \n - sum_B_mathrmfixed sum_alpha sum_i fracaleph^alpha c_iBN left( fracpartial M_A^alphapartial y_i^alpha - x_A sum_C fracpartial M_C^alphapartial y_i^alpha right)  mu_B \n + left( tildex_A - x_A right) \n\n\n\n\n\n","category":"method"},{"location":"#Calphad.unpack_indices-Tuple{Any, Any, Any}","page":"Home","title":"Calphad.unpack_indices","text":"Examples\n\nusing Symbolics\n@variables N P T N_A MU_B NP_ALPHA N\nelements = [\"A\", \"B\", \"C\"]\nphases = [\"ALPHA\", \"BETA\", \"GAMMA\"]\ncondition_dict = Dict(\n   N => 1.0,\n   P => 101325.0,\n   T => 300.0,\n   N_A => 0.5,\n   MU_B => -10000,\n   NP_ALPHA => 0.8,\n)\nunpack_indices(elements, phases, condition_dict)\n\n\n\n\n\n","category":"method"}]
}
