%% Joint_model

load('E1E2_1a.mat')

energy_vs_fitness_E1E2_new(1, msa_aa_ex, msa_aa, num_mutants_combine_array, amino_single_combine_array, J_MPF_BML, ind_non_conserve,0)


%% indept_model

load('E1_1a.mat', 'J_MPF_BML')
load('E1_1a.mat', 'num_mutants_combine_array')
load('E1_1a.mat', 'amino_single_combine_array')
load('E1_1a.mat', 'conserved')
J_MPF_BML_E1 = J_MPF_BML;
num_mutants_combine_array_E1 = num_mutants_combine_array;
amino_single_combine_array_E1 = amino_single_combine_array;
ind_non_conserve_E1 = setdiff(1:192,conserved);


% load('data_fitnessCosts_E2_99900.mat', 'J_MPF_BML')
% load('data_fitnessCosts_E2_99900.mat', 'num_mutants_combine_array')
% load('data_fitnessCosts_E2_99900.mat', 'amino_single_combine_array')
% load('data_fitnessCosts_E2_99900.mat', 'ind_non_conserve')

load('E2_1a.mat', 'J_MPF_BML')
load('E2_1a.mat', 'num_mutants_combine_array')
load('E2_1a.mat', 'amino_single_combine_array')
load('E2_1a.mat', 'conserved')
ind_non_conserve = setdiff(1:363,conserved);

J_MPF_BML_E2 = J_MPF_BML;
num_mutants_combine_array_E2 = num_mutants_combine_array;
amino_single_combine_array_E2 = amino_single_combine_array;
ind_non_conserve_E2 = ind_non_conserve+192;

load('E1E2_1a.mat', 'msa_aa')
load('E1E2_1a.mat', 'msa_aa_ex')





J_MPF_BML = blkdiag(J_MPF_BML_E1,J_MPF_BML_E2);
num_mutants_combine_array = [num_mutants_combine_array_E1 num_mutants_combine_array_E2];
amino_single_combine_array = [amino_single_combine_array_E1;amino_single_combine_array_E2];
ind_non_conserve = [ind_non_conserve_E1 ind_non_conserve_E2];

energy_vs_fitness_E1E2_new(1, msa_aa_ex, msa_aa, num_mutants_combine_array, amino_single_combine_array, J_MPF_BML, ind_non_conserve,1)
