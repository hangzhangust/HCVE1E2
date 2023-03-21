run startup.m
load('E1E2_1a.mat', 'J_MPF_BML')
load('E1E2_1a.mat', 'num_mutants_combine_array')
load('E1E2_1a.mat', 'amino_single_combine_array')
load('E1E2_1a.mat', 'ind_non_conserve')
num_mutants_combine_array_acc = cumsum(num_mutants_combine_array);


num_mutants_combine_array_acc_all = [0 num_mutants_combine_array_acc];


J = triu(J_MPF_BML,1);
[ ~, Ind ] = sort(abs(J(:)),1,'descend');

[ ind_row, ind_col ] = ind2sub(size(J),Ind); % fetch indices
ia = ind_col>ind_row;
ind_col =ind_col(ia);
ind_row =ind_row(ia);
all_pair = [];
for i =1:length(ind_row)
    all_pair = [all_pair;J(ind_row(i),ind_col(i))];
end
all_res=[];
ind_col_residue=zeros(length(ind_row),1);
ind_row_residue=zeros(length(ind_row),1);
for i =1:length(ind_row)
    pos = find(num_mutants_combine_array_acc>=ind_row(i),1);
    Start_site = num_mutants_combine_array_acc_all(pos)+1;
    End_site = num_mutants_combine_array_acc_all(pos+1);
    p = ind_row(i)-Start_site+1;
    res = amino_single_combine_array{pos,1};
    res = flip(res(2:end));
    r = res(p);
    ind_row_coverted{i} = [num2str(ind_non_conserve(pos)) r];
    ind_row_residue(i) = ind_non_conserve(pos);
end
for i =1:length(ind_col)
    pos = find(num_mutants_combine_array_acc>=ind_col(i),1);
    Start_site = num_mutants_combine_array_acc_all(pos)+1;
    End_site = num_mutants_combine_array_acc_all(pos+1);
    p = ind_col(i)-Start_site+1;
    res = amino_single_combine_array{pos,1};
    res = flip(res(2:end));
    r = res(p);
    ind_col_coverted{i} = [num2str(ind_non_conserve(pos)) r];
    ind_col_residue(i) = ind_non_conserve(pos);
end

ia = ind_col_residue>192 & ind_row_residue<=192;
pos_neg = all_pair>0;
ind_col_residue_inter = ind_col_residue(ia)+191;
ind_row_residue_inter = ind_row_residue(ia)+191;
pos_neg_inter = pos_neg(ia);

all_pair = all_pair(ia);

all_pair = all_pair(1:300);
pos_neg_inter = pos_neg_inter(1:300);
ind_col_residue_inter = ind_col_residue_inter(1:300);
ind_row_residue_inter = ind_row_residue_inter(1:300);
%%
run startup.m
ratio = [];

for nums=1:300
markersize = 6;
line_width = 1;

names={};
for i =192:746
    names = [names;num2str(i)];
end
angles =[];
for i =1:length(names)
    angles = [angles;(i-1)/length(names)*360];
end
pos_ind=[];
neg_ind =[];
name_ind=[];




pair_i =  ind_row_residue_inter;
pair_j = ind_col_residue_inter; 
pair_sign = pos_neg_inter;


name_ind = [pair_i(1:nums);pair_j(1:nums)];

for i =1:nums
   resi = pair_i(i);
   resj = pair_j(i);

   if  pair_sign(i)>0

       pos_ind = [pos_ind ;[find(strcmp(num2str(resi),names))    find(strcmp(num2str(resj),names))]];
   else   
       neg_ind = [neg_ind ;[find(strcmp(num2str(resi),names))    find(strcmp(num2str(resj),names))]];
   end


end
ratio = [ratio;size(neg_ind,1)/(size(neg_ind,1)+size(pos_ind,1))];

end



% l = categorical({'1a','1b'});
% l = reordercats(l,{'1a','1b'});
% b = bar(l,[length(cum_fre_1a) length(cum_fre_1b)],0.2,'LineWidth',0.1,'BarWidth',1);
% b.FaceColor = 'flat';
% b.CData(1,:)=purple;
% b.CData(2,:)=orange;


FIG = figure;

plot(ratio,'Color',[0.2,0.2,0.2],'LineWidth',1)
box off
set(gca,'TickDir','out')
FIG.Name = 'top_x';
set(gca,'TickLength',[0.02, 0.01])
ylabel({'Fraction of compensatory interactions'})
xlabel({'Top x pairs'})
% legend('Subtype 1a', 'Subtype 1b','Location','best')
FIG.Units = 'centimeters';


set(gcf,'Position',[10 10 8 6]);
set(gca,'Position',[.18 .17 .77 .8]);  %调整 XLABLE和YLABLE不会被切掉
% set(gca,'Position',[.22 .25 .88 .71]);  %调整 XLABLE和YLABLE不会被切掉
% set(gcf,'Position',[10 10 7.84 6]);
% set(gca,'Position',[.18 .17 .76 .74]);  %调整 XLABLE和YLABLE不会被切掉

figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(gca,'TickDir','out')

set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.18, 0.5, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
ylim([0 1])
% set(gca,'YTick', [0:3:9])
set(gca,'TickLength',[0.035, 0.03])
% print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600')
