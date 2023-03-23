% form new binary MSA based on unique sequences
load('data_fitnessCosts_E2_99900.mat')
load('E1E2_1a.mat', 'sequences')
load('E1E2_1a.mat', 'weight_seq')
msa_aa =cell2mat(sequences');


msa_aa = msa_aa(:,193:end);

msa_bin = Binary_Seq(msa_aa,amino_single_combine_array,ind_conserve);
[msa_bin_unique ind1 ind2]= unique(msa_bin,'rows');

% find new patient weighting WEIGHT_SEQ_UNIQUE based on unique sequences
for indi_bin = 1:length(ind1)
    num_term = ind2(ind1(indi_bin));
    ind_values = find(ind2==num_term);
    weight_seq_unique(indi_bin) = sum(weight_seq(ind_values));
end





[r1,r2,r3,p1,p2,p3]=verify_param_im(J_MPF_BML,msa_bin_unique,weight_seq_unique,num_mutants_combine_array);