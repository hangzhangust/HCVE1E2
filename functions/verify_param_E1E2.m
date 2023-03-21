 function z = verify_param_E1E2(Jstore_mat,msa_bin_unique,weight_seq_unique,num_mutants_combine_array)
% verify_param(Jstore,msa_bin_unique,weight_seq_unique,num_mutants_combine_array)
% 
% Verify the couplings

% Inputs:
%       Jstore - fields/couplings
%       msa_bin_unique - extended unique binary matrix
%       weight_seq_unique -  sequence weighting of sequences in msa_bin_unique
%       num_mutants_combine_array - number of mutants per residue
%       options_Verify -  options
%
% Outputs:
%       out - 1
%        
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run startup.m
Jstore = ((Jstore_mat + Jstore_mat')/2);
Jstore = Jstore(:);

%% Process Options
if nargin < 5
    options_Verify = [];
end


[no_seeds,parOpt,thin,burnin,nosim] = ...
    myProcessOptions(options_Verify,'no_seeds',1,'parOpt',0,'thin',3e3,'burnin',1e4,'nosim',1e7);


% out =          helper_verify_MCMCv3(J_MPF,msa_bin_unique,num_mutants_combine_array,thin,burnin,nosim,weight_seq_unique,num_patients,no_seeds);

num_patients = sum(weight_seq_unique);
mut_mat_MSA = full(((msa_bin_unique')*diag(weight_seq_unique)*msa_bin_unique))/num_patients; % mutant probability matrix

[delta_cij delta_cij_bound] = helper_covar(mut_mat_MSA,num_patients);


% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
msa_bin_unique = full(msa_bin_unique);
num_residues_binary = size(msa_bin_unique,2);
cumul_num_mutants_combine_array  = cumsum(num_mutants_combine_array);
protein_length_aa = length(num_mutants_combine_array);

J_MINFLOW_mat_array = full(Jstore);
J_MINFLOW = reshape(Jstore,num_residues_binary,num_residues_binary);

%%%%%%%%


number_samples = ceil((nosim-burnin)/thin);

totalnosample_par = zeros(1,no_seeds);
double_mutantsum_par = zeros(no_seeds,num_residues_binary*num_residues_binary);
num_mutant_array_par = zeros(no_seeds,number_samples);

randvalue = randi([1 size(msa_bin_unique,1)],1, no_seeds);


t_samp = tic();
% noiterations
for ite_seeds=1:no_seeds
    
    curr_vector = msa_bin_unique(randvalue(ite_seeds),:);
    
    [doublemutant nosample number_mutants]= gibbs_potts_mex(curr_vector,J_MINFLOW_mat_array,num_residues_binary,nosim,cumul_num_mutants_combine_array,num_mutants_combine_array,burnin,thin,number_samples,protein_length_aa);
    

 
    double_mutantsum_par(ite_seeds,:) =doublemutant;
    totalnosample_par(ite_seeds)= nosample;
    num_mutant_array_par(ite_seeds,:) =number_mutants;
    
    
end

if (no_seeds>1)
    double_mutantsum = sum(double_mutantsum_par);
    totalnosample = sum(totalnosample_par);
    num_mutant_array=num_mutant_array_par(:);
    
    mut_mat_MCMC_array=double_mutantsum/totalnosample;
    mut_mat_MCMC = reshape(mut_mat_MCMC_array,num_residues_binary,num_residues_binary);
    mut_mat_MCMC = (mut_mat_MCMC+mut_mat_MCMC') - diag(diag(mut_mat_MCMC));
    
else
    mut_mat_MCMC_array=double_mutantsum_par/nosample;
    mut_mat_MCMC = reshape(mut_mat_MCMC_array,num_residues_binary,num_residues_binary);
    mut_mat_MCMC = (mut_mat_MCMC+mut_mat_MCMC') - diag(diag(mut_mat_MCMC));
    num_mutant_array=num_mutant_array_par;
end


t_samp = toc(t_samp);
time_MCMC = t_samp;
fprintf( 'MCMC in %f seconds \n', t_samp );

cross_thres=1/num_patients;
[single_error double_error epsilon_max covar_error] = helper_eps(mut_mat_MSA,mut_mat_MCMC,delta_cij,delta_cij_bound,num_patients,cumul_num_mutants_combine_array,cross_thres);
single_error
double_error
set(0,'DefaultTextFontSize',8)
set(0,'DefaultAxesFontSize',8)

color = [0.6000    0.6000    0.6000];%gray
markersize = 3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot 1: Single mutant probabilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)
single_MCMC = diag(mut_mat_MCMC);
single_MSA = diag(mut_mat_MSA);




arrayline_min = min([single_MCMC-single_MSA]);
arrayline_max = max([single_MCMC-single_MSA]);
arrayline = arrayline_min:0.01:arrayline_max;
FIG = figure;
stem(single_MSA-single_MCMC,'o','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',markersize,'LineStyle','--','Color',color,'Linewidth',0.05);hold;grid off;
plot([0 length(single_MSA)],[0 0],'k')
xlabel({'Single mutant probability'})
ylabel({'Difference between', 'the SIM and the MSA'})
xticks([])
% title(['Average epsilon=' num2str(single_error)])
FIG.Name = 'single_3a';
box off
FIG.Units = 'centimeters';
set(gcf,'Position',[10 10 5 3.75]);
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(gca,'TickDir','out')
set(gca,'TickLength',[0.02, 0.03])
xlim([0 length(single_MSA)])
ylim([-0.1 0.1])
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.38, 0.5, 0]);
set(get(gca,'XLabel'), 'Units', 'Normalized', 'Position', [0.5, -0.2, 0]);
set(gca,'Position',[.3 .24 .62 .72]);  %调整 XLABLE和YLABLE不会被切掉
set(findobj('FontSize',8),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);

P = immse(single_MSA,single_MCMC)^0.5

text(50,-0.08,sprintf('RMSE = %.4f',P),'FontSize',8)





try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600');
catch 
    print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');
   
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot 2: Double mutant probabilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find location of two points in flattened array
zeros_remove = (ones(num_residues_binary,num_residues_binary));
phi_extend = [1 cumul_num_mutants_combine_array+1];
for aba=1:length(phi_extend)-1
    array_sites = phi_extend(aba):phi_extend(aba+1)-1;
    for indi=1:length(array_sites)
        for indj=1:length(array_sites)
            if (indi~=indj)
                zeros_remove(array_sites(indi),array_sites(indj))=0;
            end
        end
    end
end




z = zeros_remove;


zeros_remove = sparse(zeros_remove);

ones_double = sparse(triu(ones(num_residues_binary,num_residues_binary),1));
double_pick_final = zeros_remove.*ones_double;
double_pick_final_flat = sparse(double_pick_final(:));
ind_double_pick_final = find(double_pick_final_flat==1); % location of two points

double_MCMC_flat = mut_mat_MCMC(:);
double_MSA_flat = mut_mat_MSA(:);

double_MCMC_flat = double_MCMC_flat(ind_double_pick_final);
double_MSA_flat = double_MSA_flat(ind_double_pick_final);


arrayline_min = min([double_MCMC_flat ;double_MSA_flat]);
arrayline_max = max([double_MCMC_flat; double_MSA_flat]);
arrayline = arrayline_min:0.01:arrayline_max;
FIG = figure;
stem(double_MCMC_flat-double_MSA_flat,'o','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',markersize,'LineStyle','--','Color',color,'Linewidth',0.05);hold;grid off;
plot([0 length(double_MCMC_flat)],[0 0],'k')
xlabel({'Double mutant probability'})
ylabel({'Difference between', 'the SIM and the MSA'})
xticks([])
% title(['Average epsilon=' num2str(single_error)])
FIG.Name = 'double_3a';
box off
FIG.Units = 'centimeters';
set(gcf,'Position',[10 10 5 3.75]);
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(gca,'TickDir','out')
set(gca,'TickLength',[0.02, 0.03])
xlim([0 length(single_MSA)])
ylim([-0.1 0.1])
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.38, 0.5, 0]);
set(get(gca,'XLabel'), 'Units', 'Normalized', 'Position', [0.5, -0.2, 0]);
set(gca,'Position',[.3 .24 .62 .72]);  %调整 XLABLE和YLABLE不会被切掉
set(findobj('FontSize',8),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);

P = immse(double_MCMC_flat,double_MSA_flat)^0.5

text(50,-0.08,sprintf('RMSE = %.4f',P),'FontSize',8)
try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600');
catch 
   print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot 3: Connected correlation (covariance)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
connected_MCMC = mut_mat_MCMC - single_MCMC*single_MCMC';
connected_MSA = mut_mat_MSA - single_MSA*single_MSA';

connected_MCMC_flat = connected_MCMC(:);
connected_MSA_flat = connected_MSA(:);

connected_MCMC_flat = connected_MCMC_flat(ind_double_pick_final);
connected_MSA_flat = connected_MSA_flat(ind_double_pick_final);


arrayline_min = min(connected_MCMC_flat);
arrayline_max = max(connected_MCMC_flat);
arrayline = arrayline_min:0.01:arrayline_max;
FIG = figure;
stem(connected_MCMC_flat-connected_MSA_flat,'o','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',markersize,'LineStyle','--','Color',color,'Linewidth',0.05);hold;grid off;
plot([0 length(connected_MCMC_flat)],[0 0],'k')
xlabel({'Connected correlation'})
ylabel({'Difference between', 'the SIM and the MSA'})
xticks([])
% title(['Average epsilon=' num2str(single_error)])
FIG.Name = 'connect_3a';
box off
FIG.Units = 'centimeters';
set(gcf,'Position',[10 10 5 3.75]);
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(gca,'TickDir','out')
set(gca,'TickLength',[0.02, 0.03])
xlim([0 length(single_MSA)])
ylim([-0.1 0.1])
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.38, 0.5, 0]);
set(get(gca,'XLabel'), 'Units', 'Normalized', 'Position', [0.5, -0.2, 0]);
set(gca,'Position',[.3 .24 .62 .72]);  %调整 XLABLE和YLABLE不会被切掉
set(findobj('FontSize',8),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);

P = immse(connected_MCMC_flat,connected_MSA_flat)^0.5

text(50,-0.08,sprintf('RMSE = %.4f',P),'FontSize',8)
try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600');
catch 
   print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot 4: Probability of number of mutations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps_MSA = helper_PN(weight_seq_unique,msa_bin_unique);
max_mutant = length(eps_MSA);

for indi=1:max_mutant
    eps_MCMC(indi) = length(find(num_mutant_array==indi))/length(num_mutant_array);
end

FIG=figure
p1=semilogy(1:max_mutant,eps_MSA,'Color','k','LineWidth',1);hold;grid;
p=semilogy(1:max_mutant,eps_MCMC,'Color',color,'LineWidth',1);
% h = legend('MSA','MCMC','Location ','best');



% h.Position = [0.80 0.8324 0.1080 0.0488];
legend boxoff
% legendflex(p,{'MSA','MCMC'},'xscale',0.5 ,'anchor', {'ne','ne'}, ...
% 'buffer', [5 -5], ...
% 'ncol', 3, ...
% 'fontsize', 8, ...
% 'xscale', 0.8, ...
% 'box', 'off');
xlabel('Number of mutants, x')
ylabel('Frequency of x')
xlim([10 70])
% xlim([-0.1 no_of_nonzero_mutations(end)])
ylim([1e-4 1])
set(gca,'YTick', [1e-4 1e-2 1e-1 1])
rho = corr(eps_MSA',eps_MCMC')
epsilon_max
set(gca,'TickLength',[0.02, 0.03])
FIG.Name = 'num_3a';
grid off
box off
FIG.Units = 'centimeters';
set(gcf,'Position',[10 10 5 3.75]);
set(gca,'Position',[.25 .25 .61 .61]);  %调整 XLABLE和YLABLE不会被切掉

figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(gca,'TickDir','out')
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.31, 0.5, 0]);
set(gca,'Position',[.26 .24 .65 .72]);  %调整 XLABLE和YLABLE不会被切掉
set(findobj('FontSize',8),'FontSize',figure_FontSize);
legendflex([p1,p], {'MSA','SIM'}, 'ref', gcf, ...
                       'anchor', {'ne','ne'}, ...
                       'buffer',[-10 0], ...
                       'nrow',2, ...
                       'fontsize',8,'box','off','xscale',0.5);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(get(gca,'XLabel'), 'Units', 'Normalized', 'Position', [0.5, -0.2, 0]);
try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600');
catch 
   print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');
end
