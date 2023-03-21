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

nums=300;
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
for i = 1:size(pos_ind,1)
    alpha = angles(pos_ind(i,1));
    beta = angles(pos_ind(i,2));
    color = [purple 0.3];
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
    color = [orange 0.3];
    if abs(alpha-beta)~=180
        plot_arc(alpha,beta,color,line_width)
    else
        u  = [cos(alpha/180*pi),cos(beta/180*pi)];
        v  = [sin(alpha/180*pi);sin(beta/180*pi)];
        plot(u,v,'Color',color,'LineWidth', line_width)
    end
    
end


name_ind =[192;383;384;746];

name_ind = unique(name_ind)-191;
for i=1:length(name_ind)
    ang = angles(name_ind(i));
    v  = 1.05*r*[cos(ang/180*pi);sin(ang/180*pi)];
    align ='left';
    if name_ind(i) == 384-191

        ang = ang+2;
        ang_E1 = ang;
        v  = 1.05*r*[cos(ang/180*pi);sin(ang/180*pi)];
        txt=text(v(1),v(2),names{name_ind(i)},'HorizontalAlignment',align,'VerticalAlignment', 'middle','rotation',ang);

        continue
    end
    if name_ind(i) == 383-191

        ang = ang-2;
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




size(neg_ind,1)/(size(neg_ind,1)+size(pos_ind,1))

pos_res =[];
for i =1:length(name_ind)
    pos_res = [pos_res;str2num(names{name_ind(i)})];
end
pos_res = unique(pos_res);


set(gca,'XColor', 'none','YColor','none')
xlim([-1.4,1.4])
ylim([-1.4,1.4])
FIG.Name = join(["Top_",num2str(nums)],'');
FIG.Units = 'centimeters';
set(gcf,'Position',[5 5 10 10]);
set(gca,'Position',[.01 .01 .98 .98]);  %è°æ´ XLABLEåYLABLEä¸?ä¼è¢«åæ?
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.16, 0.5, 0]);
set(findobj('FontSize',8),'FontSize',figure_FontSize);
set(gca,'TickLength',[0.02, 0.03])

