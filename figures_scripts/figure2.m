rng(0);
run startup.m
 load('entropy.mat')

m_IM =mean(S_IM);
m_JM = mean(S_JM);
frac_IM = (S_ind-m_IM)/(S_ind-S_true);
frac_JM = (S_ind-m_JM)/(S_ind-S_true);
Frac = [frac_JM frac_IM ];


% old
% Frac = [ 0.65407 0.37611];


FIG=figure;


hold on;

% l = categorical({'1a','1b'});
% l = reordercats(l,{'1a','1b'});
% b = bar(l,[length(cum_fre_1a) length(cum_fre_1b)],0.2,'LineWidth',0.1,'BarWidth',1);
% b.FaceColor = 'flat';
% b.CData(1,:)=purple;
% b.CData(2,:)=orange;

bar(1,Frac(1), 0.5, 'FaceColor',color_scheme(1,:),'LineWidth',0.2,'FaceAlpha',0.6,'EdgeColor','none');

bar(2,Frac(2), 0.5, 'FaceColor',color_scheme(2,:),'LineWidth',0.2,'FaceAlpha',0.6,'EdgeColor','none');

hold on;

ER1 = errorbar(2,frac_IM,(m_IM-59.9471)/(S_ind-S_true),(65.6489-m_IM)/(S_ind-S_true));
hold on;
ER1.Color = [0 0 0]; 
% ER1.LineWidth = 0.25; 
ER1.Color = color_scheme(2,:);      
ER1.LineWidth = 0.5; 
ER1.LineStyle = 'none'; 





ER2 = errorbar(1,frac_JM,(m_JM-39.2233)/(S_ind-S_true),(45.3939-m_JM)/(S_ind-S_true));
ER2.Color = color_scheme(1,:);   
ER2.LineWidth = 0.5; 
ER2.LineStyle = 'none'; 



f1 = (S_ind - S_JM)/(S_ind-S_true);
f2 = (S_ind - S_IM)/(S_ind-S_true);
ranksum(f2,f1,'tail','left')

plot([1 2], [1 1]*0.8, '-k','LineWidth',0.5)

xlim([0.5 2.5])
set(gca,'XTick',1:2,'XTickLabel',...
    {'JM','IM'});
set(gca,'TickDir','out')
FIG.Name = 'peak_tract';
set(gca,'TickLength',[0.02, 0.01])
% ylabel({'I_{model}/I'})
ylabel({'Fraction of the correlated structure (FCS)', 'captured by the model'})

set(gca,'YTick',0:0.2:1)

% legend('Subtype 1a', 'Subtype 1b','Location','best')
FIG.Units = 'centimeters';


set(gcf,'Position',[10 10 6 8]);
set(gca,'Position',[.24 .17 .72 .65]);  %调整 XLABLE和YLABLE不会被切掉
% set(gca,'Position',[.22 .25 .88 .71]);  %调整 XLABLE和YLABLE不会被切掉
% set(gcf,'Position',[10 10 7.84 6]);
% set(gca,'Position',[.18 .17 .76 .74]);  %调整 XLABLE和YLABLE不会被切掉

figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.24, 0.5, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
ylim([0 1])
set(gca,'YTick', [0 0.5 1])
set(gca,'TickLength',[0.035, 0.03])
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');     end