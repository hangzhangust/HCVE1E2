%%
% *HCV E1 influences the fitness landscape of E2 and 
%  may enhance  escape from E2-specific antibodies *

% Code for developing the model and re-generating the figures in the paper

%% Setting up paths (of functions and data files required) and necessary parameters

clear;
close all;
clc;

addpath data
addpath functions
addpath 3rd_party_code
addpath figures_scripts

%% Data preparation and model training

% HCV E1E2 1a sequences were downloaded from GLUE data base
% (http://hcv-glue.cvr.gla.ac.uk; accessed March 20, 2020)
% load sequences, patient id and headers
load('cleaned_E1E2_1a.mat')




% re-weight sequences according to number of sequences from each patient
weight_seq = get_weight_seq(patient{1});

% feed the sequences and weight_seq to the MPF-BML-GUI to train the model
% Code for running MPF-BML is freely available at <https://github.com/ahmedaq/MPF-BML-GUI>. 

%% Fig. 1 Statistical validation of the inferred E1E2 fitness 
% landscape (joint model).

clear;
clc;
run figure1.m

%% Fig. 2 Comparison of the fraction of the correlated structure in
% E1E2 protein  captured by the joint model (JM) 
% and the independent model (IM).

clear;
clc;
run figure2.m

%% Fig. 3 Comparison of the E1E2 fitness prediction by JM and IM.

clear;
clc;
run figure3.m

%% Fig. 4 Strong E1E2 inter-protein interactions are largely compensatory.

clear;
clc;
run figure4.m

%% Fig. 5 Role of E1 in facilitating viral escape from E2-specific HmAbs.
% (a) Distribution of escape times of E2 residues using the inferred JM.
% (b) Comparison of escape times of E2 residues inferred from
% the JM and the E2-only model for the known E2 escape mutations (left panel) 
% and the remaining E2 residues (right panel)
% (c) Circos plot displaying the interactions between strongly-coupled residues

clear;
clc;
run figure5.m

%% Fig. 6 Evaluation of known HmAbs using the escape
% times inferred from the JM and the E2-only model.

clear;
clc;
run figure6.m

%% Supp. Fig. 1 Statistical validation of the inferred IM for the E1E2 protein.

clear;
clc;
run suppfigure1.m
%% Supp. Fig. 2 Robustness of the fraction of compensatory 
% E1E2 inter-protein interactions (Fig. 4) in the
% number of top inter-protein couplings selected

clear;
clc;
run suppfigure2.m

%% Supp. Fig. 3 Correlation between infectivity measurements and 
% predictions obtained from a site-independent model.

clear;
clc;
run suppfigure3.m

%% Supp. Fig. 4 E2-only model inferred in our previous study reproduces 
% statistics of the MSA based on the latest E2 sequence data.

clear;
clc;
run suppfigure4.m

%% Supp. Fig. 5 Comparison of the fitness prediction of the E2-only model 
% inferred in this work (left panel) and in our previous study (right panel).

clear;
clc;
run suppfigure5.m

%% Supp. Fig. 6 Estimate of the true entropy S_true


clear;
clc;
run suppfigure6.m
%% Supp. Fig. 7 Entropy calculated from different models for the E1E2 protein.

clear;
clc;
run suppfigure7.m




%% Supp. Fig. 8 Binary classifier designed to determine the optimal cut-off for 
% escape time based on experimentally or clinically identified escape mutations.


clear;
clc;
run suppfigure8.m