%%
load('E2_1a.mat')

msa_aa = msa;
msa_aa_ex = msa_bin;
ind_non_conserve = setdiff(1:363,conserved);
energy_vs_fitness_E2(1, msa_aa_ex, msa_aa, num_mutants_combine_array, amino_single_combine_array, ...
    J_MPF_BML, ind_non_conserve)



%%

load('data_fitnessCosts_E2_99900.mat')

energy_vs_fitness_E2(1, msa_aa_ex, msa_aa, num_mutants_combine_array, amino_single_combine_array, ...
    J_MPF_BML, ind_non_conserve)


