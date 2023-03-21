
load('E1E2_1a.mat', 'weight_seq')
load('E1E2_1a.mat', 'msa_aa')

N = size(msa_aa,1);

all_S = [];
num_iter = 100;

for i =[500:500:N N]
    S=[];
    for j =1:num_iter
        rng(j+i);
        N = size(msa_aa,1);
        I = randperm(N);
        
        
        weight_seq = weight_seq(I);
        msa_aa_tmp = msa_aa(I,:);
        [msa_unique ind1 ind2]= unique(msa_aa_tmp(1:i,:),'rows');
        weight_seq_unique=[];
        for indi_bin = 1:length(ind1)
            num_term = ind2(ind1(indi_bin));
            ind_values = find(ind2==num_term);
            weight_seq_unique(indi_bin) = sum(weight_seq(ind_values));
        end
        weight_seq_unique = weight_seq_unique/sum(weight_seq_unique);
        
        Stmp = sum(-weight_seq_unique.*log(weight_seq_unique));
        S = [S;Stmp];

    end
    all_S =[all_S S];
end
%%
run startup.m
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)

p = polyfit(1./[500:500:N N],mean(all_S,1),2);
FIG = figure;
scatter(1./[500:500:N N]',mean(all_S,1), 20,blue,'filled','MarkerFaceAlpha',0.6); 
xlabel('1/M')
ylabel('\langle S_{naive} \rangle')

x = 0:max(1./[500:500:N N])/100:max(1./[500:500:N N]);
y = x.^2*p(1)+x*p(2)+p(3);

hold on ;

plot(x,y,'--','Color',blue)


FIG.Units = 'centimeters';
set(gca,'Position',[.15 .2 .8 .77]);  %调整 XLABLE和YLABLE不会被切掉
% set(gca,'Position',[.1 .13 .2 .83]);  %调整 XLABLE和YLABLE不会被切掉

set(gcf,'Position',[10 10 8 6]);
% set(gca,'Position',[.1 .13 .67 .83]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(gca,'TickDir','out')
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
% set(gca,'Position',[.1 .12 .74 .87]);  %调整 XLABLE和YLABLE不会被切掉
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.145, 0.5, 0]);
set(get(gca,'XLabel'), 'Units', 'Normalized', 'Position', [0.5, -0.16, 0]);
set(get(gca,'title'), 'Units', 'Normalized', 'Position', [0.42, 1.1, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca,'TickLength',[0.02, 0.03])
