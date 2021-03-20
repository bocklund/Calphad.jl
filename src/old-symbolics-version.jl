function c_i_el(phase_record, i, el_index)
    total = 0.0
    state_vars_offset = length(phase_record.state_variables)
    for j in 1:length(phase_record.site_fractions)
        total += phase_record.inv_phase_matrix[i,j]*phase_record.mass_jac[el_index, state_vars_offset+j]
    end
    return total
end




function c_i_G(phase_record, i)
    total = 0.0
    state_vars_offset = length(phase_record.state_variables)
    for j in 1:length(phase_record.site_fractions)
        total += phase_record.inv_phase_matrix[i,j]*phase_record.grad[state_vars_offset+j]
    end
    return -total
end

function get_equilibrium_matrix(compsets)
    N = length(compsets[1].phase_rec.mass)  # assume number of elements the same everywhere
    P = length(compsets)
    eq_mat = Matrix{Num}(undef, P+N, N+P)
    for i in 1:P
        for j in 1:N
            eq_mat[i,j] = compsets[i].phase_rec.mass[j]
        end
        for j in N+1:N+P
            eq_mat[i,j] = 0.0
        end
    end
    for i in 1:N
        for j in 1:N
            term = 0.0
            for phase_idx in 1:P
                compset = compsets[phase_idx]
                num_idof = length(compset.phase_rec.site_fractions)
                state_vars_offset = length(compset.phase_rec.state_variables)
                term += compset.ℵ*sum(compset.phase_rec.mass_jac[i,state_vars_offset+idof]*c_i_el(compset.phase_rec, idof, j) for idof in 1:num_idof)
            end
            eq_mat[P+i,j] = term
        end
        for j in 1:P
            eq_mat[P+i,N+j] = compsets[j].phase_rec.mass[i]
        end
    end
    return eq_mat
end

function get_equilibrium_soln(compsets)
    N = length(compsets[1].phase_rec.mass)  # assume number of elements the same everywhere
    P = length(compsets)
    soln = Array{Num}(undef, P+N)
    for i in 1:P
        soln[i] = compsets[i].phase_rec.obj
    end
    for A in 1:N
        total = 0.0
        for ϕ in 1:P
            cs = compsets[ϕ]
            state_vars_offset = length(cs.phase_rec.state_variables)
            for _i in 1:length(cs.phase_rec.site_fractions)
                # In Sundman's 2015 paper, this term is positive (Eq. 57), but
                # with this term positive, the chemical potentials are wrong
                # (exactly flipped for a [0.5, 0.5] binary system). I think the
                # error is introduced when he inserts Eq. 43 into Eq. 53 and the
                # signs are mistakenly flipped in Eq. 54. I suspect this is due
                # to redefining the c_iG term as a negative sum in Eq. 44 after
                # Eq. 54 and beyond had been derived. The mistake seems to be
                # propapaged through the rest of the paper after Eq. 54.
                total -= cs.ℵ * cs.phase_rec.mass_jac[A,state_vars_offset+_i] * c_i_G(cs.phase_rec, _i)
                # TODO: N_A conditions need (N_A - Ñ_A) term added
            end
        end
        soln[P+A] = total
    end
    return soln
end


function solve(compsets, conditions)
    subs_dict = get_subs_dict(compsets, conditions)

    A = get_equilibrium_matrix(compsets);
    b = get_equilibrium_soln(compsets);

    AA = Symbolics.value.(substitute.(A, (subs_dict,)));
    bb = Symbolics.value.(substitute.(b, (subs_dict,)));

    # Numeric solution
    xx = AA \ bb;
    return xx
end

