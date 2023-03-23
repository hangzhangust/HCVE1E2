%%
run startup.m
load('entropy.mat')
% old 67.4
m_IM = mean(S_IM);
m_JM = mean(S_JM);

Entropy = [S_ind m_JM  m_IM  S_true];
% Entropy = [103.2 40.9 67.4 7.9];
FIG=figure;


hold on;


bar(1,Entropy(1), 0.8, 'FaceColor',color_scheme(9,:),'LineWidth',0.2,'FaceAlpha',0.6,'EdgeColor','none');

bar(3,Entropy(2), 0.8, 'FaceColor',color_scheme(1,:),'LineWidth',0.2,'FaceAlpha',0.6,'EdgeColor','none');

bar(2,Entropy(3), 0.8, 'FaceColor',color_scheme(2,:),'LineWidth',0.2,'FaceAlpha',0.6,'EdgeColor','none');

bar(4,Entropy(4), 0.8, 'FaceColor',color_scheme(4,:),'LineWidth',0.2,'FaceAlpha',0.6,'EdgeColor','none');


ER1 = errorbar(2,m_IM,m_IM-59.9471,65.6489-m_IM);
hold on;
ER1.Color = [0 0 0]; 
% ER1.LineWidth = 0.25; 
ER1.Color = color_scheme(2,:);      
ER1.LineWidth = 0.5; 
ER1.LineStyle = 'none'; 





ER2 = errorbar(3,m_JM,m_JM-39.2233,45.3939-m_JM);
ER2.Color = color_scheme(1,:);   
ER2.LineWidth = 0.5; 
ER2.LineStyle = 'none'; 
% b.CData(2,:)=orange;
% for i =1:4
%     bar(i,Entropy(i), 0.8, 'FaceColor',color_scheme(i,:),'LineWidth',0.2);
% end
xlim([0.5 4.5])
set(gca,'XTick',1:4,'XTickLabel',...
    {'S_{ind}','S_{IM}','S_{JM}','S_{true}'});
set(gca,'TickDir','out')
FIG.Name = 'peak_tract';
set(gca,'TickLength',[0.02, 0.01])
ylabel({'Entropy (nats)'})
% legend('Subtype 1a', 'Subtype 1b','Location','best')
FIG.Units = 'centimeters';


set(gcf,'Position',[10 10 8 10]);
set(gca,'Position',[.18 .17 .8 .74]);  %调整 XLABLE和YLABLE不会被切掉
% set(gca,'Position',[.22 .25 .88 .71]);  %调整 XLABLE和YLABLE不会被切掉
% set(gcf,'Position',[10 10 7.84 6]);
% set(gca,'Position',[.18 .17 .76 .74]);  %调整 XLABLE和YLABLE不会被切掉

figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
ylim([0 110])
set(gca,'YTick', [0 50 100])
set(gca,'TickLength',[0.035, 0.03])