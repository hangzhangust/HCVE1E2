run startup.m
load('escape_time_1a_500.mat', 'mean_escape_time')

load('mean_escape_time_E1E2.mat', 'mean_escape_time_E1E2')
escape_mutations = [384;386;388;390;391;393;394;...
    395;396;397;398;399;400;401;402;403;404;405;407;408;410;415;416;417;422;424;431;433;...
    434;435;438;442;444;446;453;456;461;466;475;482;501;524;528;531;533;538;557;558;560;580;608;610;636;713;520]';

remaining_sites = setdiff(1:363,escape_mutations-383);

% 

data = [mean_escape_time_E1E2(escape_mutations-191) mean_escape_time_E1E2(remaining_sites+192)];


G = [zeros(size(escape_mutations)) 2*ones(size(remaining_sites)) ];


set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)
box_lineWidth = 0.75;
box_widths_value = 0.5;
black = [0 0 0];
box_color = [black;black;];
box_color_transparency = 0; %faceAlpha
median_lineWidth = 0.75;
median_color = 'k';
whisker_value = 1.5;
outlier_marker = '';
outlier_markerSize = 3.5;
outlier_marker_edgeWidth = 0.001;
outlier_marker_edgeColor = 'w';
outlier_jitter_value = 0;
label_xaxis_data = {'Escape',sprintf('Remaining')};
text_ylabel = 'Escape time';
text_xlabel = '';
text_title = '';%'E2-escape mutations [Keck2009],[Morin2012],[Bailey2015]';
label_orientation_choice = 'horizontal'; %'horizontal'
ylim_min = 0;
ylim_max = 0.2;
savefig = 0;
savefig_name = 'escape_mutations';
fig_width_cm = 4;
fig_height_cm = 5;
FIG=figure;
set(gcf,'renderer','Painters')


hold on
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)
size_marker=10;
dots = mean_escape_time_E1E2(escape_mutations-191) ;
nbins =55;
max_range = 0.65;
center = 1;
interval = 15;
scale =1;
[x_data, y_data] = dot_boxplot(dots,nbins,center,max_range,scale,interval);
% f1=scatter(x1,all_E(polymorphisms_associated_with_neutralization_resistance-383) ,'o','MarkerEdgeColor','w','MarkerFaceColor',blue,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on 
f1=scatter(x_data,y_data ,'o','MarkerEdgeColor','w','MarkerFaceColor',blue,'SizeData',size_marker);f1.MarkerFaceAlpha = 0.5;hold on 


hold on;
dots = mean_escape_time_E1E2(remaining_sites+192) ;
nbins =100;
max_range = 0.45;
center = 2;
interval = 20;
scale =1;
[x_data, y_data] = dot_boxplot(dots,nbins,center,max_range,scale,interval);
% f1=scatter(x1,all_E(polymorphisms_associated_with_neutralization_resistance-383) ,'o','MarkerEdgeColor','w','MarkerFaceColor',blue,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on 

ind = (x_data<0.86+2 | x_data>1.14+2) & y_data ==max(y_data);
x_data = x_data(~ind);
y_data = y_data(~ind);
f3 = scatter(x_data,y_data ,'o','MarkerEdgeColor','w','MarkerFaceColor',red,'SizeData',size_marker);f3.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on 

figure_boxplot(data,G,...
    box_lineWidth,box_widths_value,box_color,box_color_transparency,...
    median_lineWidth,median_color,...
    whisker_value,...
    outlier_marker,outlier_markerSize,outlier_marker_edgeWidth,outlier_marker_edgeColor,outlier_jitter_value,...
    label_xaxis_data,text_ylabel,text_xlabel,text_title,label_orientation_choice,...
    ylim_min,ylim_max,...
    savefig,savefig_name,fig_width_cm,fig_height_cm);

xtick =set(gca,'XTick',1:2,'XTickLabel',...
    {'Escape','Remaining'},'FontName','Arial');

set(gca,'YTick',50:150:600)
yt = get(gca, 'YTick');
xlim([0.5 2.5])
axis([xlim    50  600])
xt = get(gca, 'XTick');
hold on
% plot(xt([3 4]), [1 1]*550, '-k','LineWidth',0.5)
% plot(xt([1 2]), [1 1]*600, '-k','LineWidth',0.5)
% plot(xt([1 3]), [1 1]*600, '-k','LineWidth',0.5)
% plot(xt([2 4]), [1 1]*650, '-k','LineWidth',0.5)
% plot(xt([1 1]), [0.95 1]*max(yt)*1.05, '-k','LineWidth',0.5)
% plot(xt([2 2]), [0.95 1]*max(yt)*1.05, '-k','LineWidth',0.5)
% P = ranksum(mean_escape_time_E1E2(remaining_sites+192),mean_escape_time(remaining_sites),'tail','left')
P = ranksum(mean_escape_time_E1E2(escape_mutations-191),mean_escape_time_E1E2(remaining_sites+192),'tail','left')
% P = ranksum(mean_escape_time(escape_mutations-383),mean_escape_time(remaining_sites),'tail','left')



% text(1.1,max(yt)*1.1,'$$ P = 1.3 \times 10^{-18} $$','interpreter','latex','FontSize',12)
% text(0.8,-100,'residues','FontSize',8)


FIG.Name = 'Escape_mutation';
ylabel({'Escape time'})

FIG.Units = 'centimeters';
% set(gcf,'Position',[6.53 6.53 5 14.35]);
box off
% text(150,60,sprintf('r = %.2f',r),'FontSize',8)
% xlabel({'Escape time of the joint model'})
% ylabel({'Escape time of the','E2-only model'})
% set(gcf,'Position',[6.53 6.53 5.22 3.92]);
% set(gca,'Position',[.23 .08 .65 .88]);  %调整 XLABLE和YLABLE不会被切掉
% set(gca,'Position',[.3 .2 .6 .78]);  %调整 XLABLE和YLABLE不会被切掉

% set(gca,'Position',[.27 .2 .7 .75]);  %调整 XLABLE和YLABLE不会被切掉
set(gcf,'Position',[6.53 6.53 4.5 6]);
box off
% text(150,60,sprintf('r = %.2f',r),'FontSize',8)
% xlabel({'Escape time of the joint model'})
% ylabel({'Escape time of the','E2-only model'})
% set(gcf,'Position',[6.53 6.53 5.22 3.92]);
% set(gca,'Position',[.23 .08 .65 .88]);  %调整 XLABLE和YLABLE不会被切掉
% set(gca,'Position',[.3 .2 .6 .78]);  %调整 XLABLE和YLABLE不会被切掉

set(gca,'Position',[.27 .2 .7 .73]);  %调整 XLABLE和YLABLE不会被切掉
set(gca,'Position',[.27 .2 .7 .73]);  %调整 XLABLE和YLABLE不会被切掉
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.27, 0.5, 0]);
set(get(gca,'XLabel'), 'Units', 'Normalized', 'Position', [0.5, -0.13, 0]);
set(gca,'TickLength',[0.02, 0.01])
FIG.Name = 'Escape_mutation';
set(gca,'TickDir','out')
%%
% 
data = [mean_escape_time_E1E2(escape_mutations-191) mean_escape_time(escape_mutations-383)];


G = [zeros(size(escape_mutations)) ones(size(escape_mutations))];


set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)
box_lineWidth = 0.75;
box_widths_value = 0.5;
black = [0 0 0];
box_color = [black;black;];
box_color_transparency = 0; %faceAlpha
median_lineWidth = 0.75;
median_color = 'k';
whisker_value = 1.5;
outlier_marker = '';
outlier_markerSize = 3.5;
outlier_marker_edgeWidth = 0.001;
outlier_marker_edgeColor = 'w';
outlier_jitter_value = 0;
label_xaxis_data = {'JM',sprintf('E2-only')};
text_ylabel = 'Escape time';
text_xlabel = '';
text_title = '';%'E2-escape mutations [Keck2009],[Morin2012],[Bailey2015]';
label_orientation_choice = 'horizontal'; %'horizontal'
ylim_min = 0;
ylim_max = 0.2;
savefig = 0;
savefig_name = 'escape_mutations';
fig_width_cm = 4;
fig_height_cm = 5;
FIG=figure;
set(gcf,'renderer','Painters')


hold on
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)
size_marker=10;
dots = mean_escape_time_E1E2(escape_mutations-191) ;
nbins =55;
max_range = 0.65;
center = 1;
interval = 15;
scale =1;
[x_data, y_data] = dot_boxplot(dots,nbins,center,max_range,scale,interval);
% f1=scatter(x1,all_E(polymorphisms_associated_with_neutralization_resistance-383) ,'o','MarkerEdgeColor','w','MarkerFaceColor',blue,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on 
f1=scatter(x_data,y_data ,'o','MarkerEdgeColor','w','MarkerFaceColor',blue,'SizeData',size_marker);f1.MarkerFaceAlpha = 0.5;hold on 


hold on;

dots =  mean_escape_time(escape_mutations-383);
center = 2;
[x_data, y_data] = dot_boxplot(dots,nbins,center,max_range,scale,interval);
% f2=scatter(x2,all_E(setdiff(1:L,polymorphisms_associated_with_neutralization_resistance-383)),'o','MarkerEdgeColor','w','MarkerFaceColor',red,'SizeData',size_marker,'LineWidth',0.01);f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on
f2=scatter(x_data,y_data,'o','MarkerEdgeColor','w','MarkerFaceColor',blue,'SizeData',size_marker);f2.MarkerFaceAlpha = 1;hold on
hold on


figure_boxplot(data,G,...
    box_lineWidth,box_widths_value,box_color,box_color_transparency,...
    median_lineWidth,median_color,...
    whisker_value,...
    outlier_marker,outlier_markerSize,outlier_marker_edgeWidth,outlier_marker_edgeColor,outlier_jitter_value,...
    label_xaxis_data,text_ylabel,text_xlabel,text_title,label_orientation_choice,...
    ylim_min,ylim_max,...
    savefig,savefig_name,fig_width_cm,fig_height_cm);

xtick =set(gca,'XTick',1:2,'XTickLabel',...
    {'JM','E2-only'},'FontName','Arial');

set(gca,'YTick',50:50:200)
yt = get(gca, 'YTick');
xlim([0.5 2.5])
axis([xlim    50  250])
xt = get(gca, 'XTick');
hold on
% plot(xt([3 4]), [1 1]*550, '-k','LineWidth',0.5)
% plot(xt([1 2]), [1 1]*220, '-k','LineWidth',0.5)
% plot(xt([1 3]), [1 1]*600, '-k','LineWidth',0.5)
% plot(xt([2 4]), [1 1]*650, '-k','LineWidth',0.5)
% plot(xt([1 1]), [0.95 1]*max(yt)*1.05, '-k','LineWidth',0.5)
% plot(xt([2 2]), [0.95 1]*max(yt)*1.05, '-k','LineWidth',0.5)
% P = ranksum(mean_escape_time_E1E2(remaining_sites+192),mean_escape_time(remaining_sites),'tail','left')
% P = ranksum(mean_escape_time_E1E2(escape_mutations-191),mean_escape_time_E1E2(remaining_sites+192),'tail','left')
% P = ranksum(mean_escape_time(escape_mutations-383),mean_escape_time(remaining_sites),'tail','left')
% 
P = ranksum(mean_escape_time_E1E2(escape_mutations-191),mean_escape_time(escape_mutations-383),'tail','left')


% text(1.1,max(yt)*1.1,'$$ P = 1.3 \times 10^{-18} $$','interpreter','latex','FontSize',12)
% text(0.8,-100,'residues','FontSize',8)
hold on
text(1.7,-25,'model','FontName','Arial','FontSize',8)

FIG.Name = 'Escape_compare';
ylabel({'Escape time'})

FIG.Units = 'centimeters';
set(gcf,'Position',[6.53 6.53 4.5 6]);
box off
% text(150,60,sprintf('r = %.2f',r),'FontSize',8)
% xlabel({'Escape time of the joint model'})
% ylabel({'Escape time of the','E2-only model'})
% set(gcf,'Position',[6.53 6.53 5.22 3.92]);
% set(gca,'Position',[.23 .08 .65 .88]);  %调整 XLABLE和YLABLE不会被切掉
% set(gca,'Position',[.3 .2 .6 .78]);  %调整 XLABLE和YLABLE不会被切掉

set(gca,'Position',[.27 .2 .7 .73]);  %调整 XLABLE和YLABLE不会被切掉
set(gca,'Position',[.27 .2 .7 .73]);  %调整 XLABLE和YLABLE不会被切掉
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.27, 0.5, 0]);
set(get(gca,'XLabel'), 'Units', 'Normalized', 'Position', [0.5, -0.13, 0]);
set(gca,'TickLength',[0.02, 0.01])
set(gca,'TickLength',[0.02, 0.01])
set(gca,'TickDir','out')
%%
data = [mean_escape_time_E1E2(remaining_sites+192) mean_escape_time(remaining_sites)];


G = [2*ones(size(remaining_sites)) 3*ones(size(remaining_sites))];


set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)
box_lineWidth = 0.75;
box_widths_value = 0.5;
black = [0 0 0];
box_color = [black;black;];
box_color_transparency = 0; %faceAlpha
median_lineWidth = 0.75;
median_color = 'k';
whisker_value = 1.5;
outlier_marker = '';
outlier_markerSize = 3.5;
outlier_marker_edgeWidth = 0.001;
outlier_marker_edgeColor = 'w';
outlier_jitter_value = 0;
label_xaxis_data = {'JM',sprintf('E2-only')};
text_ylabel = 'Escape time';
text_xlabel = '';
text_title = '';%'E2-escape mutations [Keck2009],[Morin2012],[Bailey2015]';
label_orientation_choice = 'horizontal'; %'horizontal'
ylim_min = 0;
ylim_max = 0.2;
savefig = 0;
savefig_name = 'escape_mutations';
fig_width_cm = 4;
fig_height_cm = 5;
FIG=figure;
set(gcf,'renderer','Painters')


hold on
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)
size_marker=10;
hold on;
dots = mean_escape_time_E1E2(remaining_sites+192) ;
nbins =100;
max_range = 0.45;
center = 1;
interval = 20;
scale =1;
[x_data, y_data] = dot_boxplot(dots,nbins,center,max_range,scale,interval);
% f1=scatter(x1,all_E(polymorphisms_associated_with_neutralization_resistance-383) ,'o','MarkerEdgeColor','w','MarkerFaceColor',blue,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on 

ind = (x_data<0.86 | x_data>1.14) & y_data ==max(y_data);
x_data = x_data(~ind);
y_data = y_data(~ind);
f3 = scatter(x_data,y_data ,'o','MarkerEdgeColor','w','MarkerFaceColor',red,'SizeData',size_marker);f3.MarkerFaceAlpha = 0.5;hold on 

dots =  mean_escape_time(remaining_sites);
center = 2;
[x_data, y_data] = dot_boxplot(dots,nbins,center,max_range,scale,interval);
ind = (x_data<0.86+1 | x_data>1.14+1) & y_data ==max(y_data);
x_data = x_data(~ind);
y_data = y_data(~ind);
% f2=scatter(x2,all_E(setdiff(1:L,polymorphisms_associated_with_neutralization_resistance-383)),'o','MarkerEdgeColor','w','MarkerFaceColor',red,'SizeData',size_marker,'LineWidth',0.01);f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on
f4 =scatter(x_data,y_data,'o','MarkerEdgeColor','w','MarkerFaceColor',red,'SizeData',size_marker);f4.MarkerFaceAlpha = 1;hold on
hold on

figure_boxplot(data,G,...
    box_lineWidth,box_widths_value,box_color,box_color_transparency,...
    median_lineWidth,median_color,...
    whisker_value,...
    outlier_marker,outlier_markerSize,outlier_marker_edgeWidth,outlier_marker_edgeColor,outlier_jitter_value,...
    label_xaxis_data,text_ylabel,text_xlabel,text_title,label_orientation_choice,...
    ylim_min,ylim_max,...
    savefig,savefig_name,fig_width_cm,fig_height_cm);

xtick =set(gca,'XTick',1:2,'XTickLabel',...
    {'JM','E2-only'},'FontName','Arial','FontSize',8);

set(gca,'YTick',50:150:600)
yt = get(gca, 'YTick');
xlim([0.5 2.5])
axis([xlim    50  600])
xt = get(gca, 'XTick');
hold on
% plot(xt([3 4]), [1 1]*550, '-k','LineWidth',0.5)
% plot(xt([1 2]), [1 1]*550, '-k','LineWidth',0.5)
% plot(xt([1 3]), [1 1]*600, '-k','LineWidth',0.5)
% plot(xt([2 4]), [1 1]*650, '-k','LineWidth',0.5)
% plot(xt([1 1]), [0.95 1]*max(yt)*1.05, '-k','LineWidth',0.5)
% plot(xt([2 2]), [0.95 1]*max(yt)*1.05, '-k','LineWidth',0.5)
P = ranksum(mean_escape_time_E1E2(remaining_sites+192),mean_escape_time(remaining_sites),'tail','left')
% P = ranksum(mean_escape_time_E1E2(escape_mutations-191),mean_escape_time_E1E2(remaining_sites+192),'tail','left')
% P = ranksum(mean_escape_time(escape_mutations-383),mean_escape_time(remaining_sites),'tail','left')



% text(1.1,max(yt)*1.1,'$$ P = 1.3 \times 10^{-18} $$','interpreter','latex','FontSize',12)
text(1.7,-50,'model','FontName','Arial','FontSize',8)


FIG.Name = 'remain_compare';
ylabel({'Escape time'})

FIG.Units = 'centimeters';
set(gcf,'Position',[6.53 6.53 4.5 6]);
box off
% text(150,60,sprintf('r = %.2f',r),'FontSize',8)
% xlabel({'Escape time of the joint model'})
% ylabel({'Escape time of the','E2-only model'})
% set(gcf,'Position',[6.53 6.53 5.22 3.92]);
% set(gca,'Position',[.23 .08 .65 .88]);  %调整 XLABLE和YLABLE不会被切掉
% set(gca,'Position',[.3 .2 .6 .78]);  %调整 XLABLE和YLABLE不会被切掉

set(gca,'Position',[.27 .2 .7 .73]);  %调整 XLABLE和YLABLE不会被切掉
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.27, 0.5, 0]);
set(get(gca,'XLabel'), 'Units', 'Normalized', 'Position', [0.5, -0.13, 0]);
set(gca,'TickLength',[0.02, 0.01])
set(gca,'TickDir','out')


%%
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


nums=[300];
markersize = 6;
line_width = 2;

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

   if  isempty(intersect(resj,escape_mutations)) && isempty(intersect(resi,escape_mutations)) 

       pos_ind = [pos_ind ;[find(strcmp(num2str(resi),names))    find(strcmp(num2str(resj),names))]];
   else   
       neg_ind = [neg_ind ;[find(strcmp(num2str(resi),names))    find(strcmp(num2str(resj),names))]];
   end


end
r=1.1;
ang1 = [0:1:124];
ang2 = [125:1:359];
% x= r*sin([0:1:360]/360*2*pi);
% y = r*cos([0:1:360]/360*2*pi);
x1= r*sin(ang1/180*pi);
y1 = r*cos(ang1/180*pi);
x2= r*sin(ang2/180*pi);
y2 = r*cos(ang2/180*pi);


set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)



    
FIG=figure;
plot(y1,x1,'Color',[0.2 0.2 0.2],'LineWidth',5);
hold on
plot(y2,x2,'Color',[0.8 0.8 0.8],'LineWidth',5);
hold on;

line_width = 1;

for i = 1:size(pos_ind,1)
    alpha = angles(pos_ind(i,1));
    beta = angles(pos_ind(i,2));
    color = [[220 213 185]/256 0.5];
    if abs(alpha-beta)~=180
    plot_arc(alpha,beta,color,line_width)
    else
        u  = [cos(alpha/180*pi);sin(alpha/180*pi)];
        v  = [cos(beta/180*pi);sin(beta/180*pi)];
        plot(u,v,'Color',color,'LineWidth', line_width)
    end    
    
end
for i = 1:size(neg_ind,1)
    alpha = angles(neg_ind(i,1));
    beta = angles(neg_ind(i,2));
    color = [blue 0.5];
    if abs(alpha-beta)~=180
        plot_arc(alpha,beta,color,line_width)
    else
        u  = [cos(alpha/180*pi),cos(beta/180*pi)];
        v  = [sin(alpha/180*pi);sin(beta/180*pi)];
        plot(u,v,'Color',color,'LineWidth', line_width)
    end
    
end
name_ind =[192;383;384;746];
% name_ind =intersect(name_ind,region_E2{7,1});
name_ind = unique(name_ind)-191;
for i=1:length(name_ind)
    ang = angles(name_ind(i));
    v  = 1.05*r*[cos(ang/180*pi);sin(ang/180*pi)];
    align ='left';
    if name_ind(i) == 384-191

        ang = ang+3;
        ang_E1 = ang;
        v  = 1.05*r*[cos(ang/180*pi);sin(ang/180*pi)];
        txt=text(v(1),v(2),names{name_ind(i)},'HorizontalAlignment',align,'VerticalAlignment', 'middle','rotation',ang);

        continue
    end
    if name_ind(i) == 383-191

        ang = ang-3;
        ang_E1 = ang;
        v  = 1.05*r*[cos(ang/180*pi);sin(ang/180*pi)];
        txt=text(v(1),v(2),names{name_ind(i)},'HorizontalAlignment',align,'VerticalAlignment', 'middle','rotation',ang);

        continue
    end
    if name_ind(i) == 192-191

        ang = ang+3;
        ang_E1 = ang;
        v  = 1.05*r*[cos(ang/180*pi);sin(ang/180*pi)];
        txt=text(v(1),v(2),names{name_ind(i)},'HorizontalAlignment',align,'VerticalAlignment', 'middle','rotation',ang);

        continue
    end

    if name_ind(i) == 746-191

        ang = ang-3;
        ang_E1 = ang;
        v  = 1.05*r*[cos(ang/180*pi);sin(ang/180*pi)];
        txt=text(v(1),v(2),names{name_ind(i)},'HorizontalAlignment',align,'VerticalAlignment', 'middle','rotation',ang);

        continue
    end
    txt=text(v(1),v(2),names{name_ind(i)},'HorizontalAlignment',align,'VerticalAlignment', 'middle','rotation',ang);
    
end

ang = 45;
v  = 1.05*r*[cos(ang/180*pi);sin(ang/180*pi)];
txt=text(v(1),v(2),'E1','HorizontalAlignment',align,'VerticalAlignment', 'middle','rotation',ang);

ang = 225;
v  = 1.05*r*[cos(ang/180*pi);sin(ang/180*pi)];
txt=text(v(1),v(2),'E2','HorizontalAlignment',align,'VerticalAlignment', 'middle','rotation',ang);


ang1 = [10:1:40];
ang2 = [135:220 ];
% x= r*sin([0:1:360]/360*2*pi);
% y = r*cos([0:1:360]/360*2*pi);
x1= 1.1*r*sin(ang1/180*pi);
y1 = 1.1*r*cos(ang1/180*pi);
x2= 1.1*r*sin(ang2/180*pi);
y2 = 1.1*r*cos(ang2/180*pi);
plot(y1,x1,'Color',[0 0 0],'LineWidth',1);
hold on
plot(y2,x2,'Color',[0 0 0],'LineWidth',1);

ang1 = [50:114];
ang2 = [230:350 ];
% x= r*sin([0:1:360]/360*2*pi);
% y = r*cos([0:1:360]/360*2*pi);
x1= 1.1*r*sin(ang1/180*pi);
y1 = 1.1*r*cos(ang1/180*pi);
x2= 1.1*r*sin(ang2/180*pi);
y2 = 1.1*r*cos(ang2/180*pi);
plot(y1,x1,'Color',[0 0 0],'LineWidth',1);
hold on
plot(y2,x2,'Color',[0 0 0],'LineWidth',1);





pos_res =[];
for i =1:length(name_ind)
    pos_res = [pos_res;str2num(names{name_ind(i)})];
end
pos_res = unique(pos_res);


set(gca,'XColor', 'none','YColor','none')
xlim([-1.4,1.4])
ylim([-1.4,1.4])
FIG.Name = join(["escape_mutation_",num2str(nums)],'');
FIG.Units = 'centimeters';
set(gcf,'Position',[5 5 6 6]);
set(gca,'Position',[.01 .01 .98 .98]);  %è°æ´ XLABLEåYLABLEä¸?ä¼è¢«åæ?
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.16, 0.5, 0]);
set(findobj('FontSize',8),'FontSize',figure_FontSize);
set(gca,'TickLength',[0.02, 0.03])
hold on
% text(0.5,-1.3,sprintf('p = 7.7e-20'),'FontSize',8)
text(0.5,-1.3,sprintf('p = 5.4e-19'),'FontSize',8)