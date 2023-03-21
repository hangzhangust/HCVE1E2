%%
HmAbs_all=[];
load('mean_escape_time_E1E2.mat')
HmAb_Pierce2016 = [];
run startup.m
%CBH-4D
HmAb_Pierce2016{1} = [494 497 502 504 506:509 511 537 539 540 542:545 547 549:552 554 556 559 561 562 564 565 584 585 592 594 598 600 602 603 607:611 614 617:619 621 623 624 626 627 629:633 638 640 642:644];
%CBH-4G
HmAb_Pierce2016{2} = [494 497 503 505:509 511 517 537 539 540 542:545 547 549:552 554 556 559 561 564 565 584 592:594 597 598 600 602 603 607:611 614 617:619 621 623 624 627 629 631:633 638 640 642:644];
%CBH-4B
HmAb_Pierce2016{3} = [494 497 503 505:509 537 539 540 550:552 554 559 561 564 565 600 602 603 607:611 614 617:619 621 624 627 629 631:633 638 640 642:644];
%CBH-20
% HmAb_Pierce2016{4} = [494 497 503 505:509 537 539 540 542:545 547 549:552 554 559 561 564 565 569 581 583:585 592 594 597 598 600 602 603 607 608 611 614 617:619 621 623 624 627 631 632 638 640 642:644];
%included the residues with close to threshold (RB = 21,22) as well
HmAb_Pierce2016{4} = [494 497 503 505:509 537 539 540 542:545 547 549:552 554 559 561 564 565 569 581 583:585 592 594 597 598 600 602 603 607 608 610 611 614 617:619 621 623 624 627 631 632 638 640 642:644]; 
%CBH-21
HmAb_Pierce2016{5} = [452 494 497 503 505:509 537 539 540 542:545 547 549:552 554 559 561 564 565 569 581 583:585 592:594 597 598 600 602 603 607:611 614 617:619 621 623 624 627 631:633 638 640 642:644];
%CBH-22
HmAb_Pierce2016{6} = [452 494 497 503 505:509 537 539 540 542:545 547 549:552 554 559 561 564 565 569 581 583:585 592:594 597 598 600 602 603 607:611 614 617:619 621 623 624 627 631:633 638 640 642:644];

%HC-1
HmAb_Pierce2016{7} = [429 494 503:506 508 509 529 530 535 537 539 552 554 559 564 607 611 614 617 644];
%HC-11
HmAb_Pierce2016{8} = [425 428 429 436:438 442 443 494 497 502:504 506:509 511 520 530 535 537 539 550:552 554 556 559 564 565 602 603 607 608 611 614 617:621 624 640 643 644];
%A27
HmAb_Pierce2016{9} = [424:429 437 438 494 497 499 502:504 506:509 511 520 529 530 535 537 539 540 550:552 554 556 559 564 565 602 603 607:611 614 616:619 621 623 624 638 640 642:644];

%CBH-23
HmAb_Pierce2016{10} = [494 508 509 537 539 549 552 554 564 611 614 644]; %same as CBH7
%CBH-7
HmAb_Pierce2016{11} = [494 506 508 509 537 539 549 552 554 564 611 614 617 621 644];

%HC84-20
HmAb_Pierce2016{12} = [429 441 494 497 502:509 511 537 539 552 554 559 564 607 608 611 613 614 616:619 621 640 643 644];
%HC84-24
HmAb_Pierce2016{13} = [429 442 443 494 497 502:509 511 537 539 552 554 559 564 607 608 611 614 617:619 621 643 644];
%HC84-26
HmAb_Pierce2016{14} = [441 442 494 497 502:509 511 537 539 552 554 559 564 603 607 608 611 614 617:619 621 640 643 644];

%HC33-1
HmAb_Pierce2016{15} = [413 418 420];
%HC33-4
HmAb_Pierce2016{16} = [408 413 420];
% %CD81_bs
% HmAb_Pierce2016{17} = [420 421 424 427 430 436:438 440:443 523 527 529 530 535 540 613 614 616:618];



%HCV1
HmAbs_Gopal2017{1} = [413 415 417 418 420 422 624 639];
%AR1A
HmAbs_Gopal2017{2} = [417 485 494 497 502 504 506:511 514 516 519 537:540 542 544 545 547 549:552 554 559 564 565 602 603 607:611 614 617:619 621 623 624 632 638 640 642 643 644]; %559 is RB=21%
%AR1B
HmAbs_Gopal2017{3} = [494 497 504 506:509 511 537 539 544 545 547:552 554 564 607:608 610 611 614 617:619 621 640 644]; %610 is RB=21%
%AR2A
HmAbs_Gopal2017{4} = [552 607:611 614 617:619 621 623:625 628 638 640 643 644];
%AR3A
HmAbs_Gopal2017{5} = [424 425 428 429 436:438 440:442 485 494 496 497 502:509 511 516 523 525 529 530 535 537 539 552 554 559 564 602 603 607:611 614 616:619 621 623 624 625 640 642:644];
%AR3B
HmAbs_Gopal2017{6} = [424 425 427:429 432 436 437 440:442 494 496 497 502:509 511 516 518 520 523 529 530 535 537 539 540 552 554 555 556 558 559 564 602 603 607:611 614 616:619 621 623 624 625 630 638 640:644];
%AR3C
HmAbs_Gopal2017{7} = [424 425 428:429 437 438 441:443 494 496 497 502:509 511 516 525 529 530 535 537 539 540 552 554 559 564 602 603 607:611 614 616:619 621 623 624 625 630 638 640 641 643 644];
%AR3D
HmAbs_Gopal2017{8} = [424 425 427:430 436:438 440:442 459 494 496 497 499 502:509 511 516 518 520 523 529 530 535 537 539 540 550:552 554 555 556 558 559 564 565 600 602 603 607:612 614 616:619 621 623 624 625 638 640:644];
 
 

%CBH-5
HmAbs_Keck2019{1} = [424:429 436:438 441:443 535 638];
%212.1.1
HmAbs_Keck2019{2} = [425:429 433 434 529 530 535];
%212.1.10
HmAbs_Keck2019{3} = [426 428 429 442];
%212.15
HmAbs_Keck2019{4} = [544 547 549 637 638 639];
%212.25
HmAbs_Keck2019{5} = [549 636 638 639];

 
% Bailey2017

%HEPC3
HmAbs_Bailey2017{1} = [425 427 428 437 499 520 530 535];
%HEPC43
HmAbs_Bailey2017{2} = [425 427 428 432 436 437 438 442 443 499 517 520 527 529 530 535 616];
%HEPC74
HmAbs_Bailey2017{3} = [425 428 436 437 530 535];
%HEPC46
HmAbs_Bailey2017{4} = [541:546 548 549 594 598 633];
%HEPC50
HmAbs_Bailey2017{5} = [543 544 545 549 594 597 598];
%HEPC98
HmAbs_Bailey2017{6} = [402 405 408];



% % A4
% HmAbs_new{1} = [199 201 202];


% IGH526
HmAbs_new{1} = [316, 320, 321, 323, 324];

%AR4A
HmAbs_new{2} = [201, 205, 207, 212, 222, 226, 228, 229, 238, 239, 304, 306, 459, 486, 487, 494, 497, 504, 506, 507, 508, 509, 511, 537, 539, 543, 545, 550, 551, 552, 554, 559, 564, 565, 569, 585, 594, 597, 600, 602, 607, 608, 611, 614, 617, 618, 619, 621, 635, 640, 643, 644, 652, 657, 677, 679, 698];


%AR5A
HmAbs_new{3} = [201, 205, 207, 212, 220, 222, 226, 228, 229, 238, 304, 306, 459, 486, 487, 494, 497, 504, 506, 507, 508, 509, 511, 513, 537, 539, 543, 544, 550, 551, 552, 554, 559, 564, 565, 569, 573, 585, 594, 597, 600, 602, 603, 607, 608, 610, 611, 614, 617, 618, 619, 621, 623, 635, 638, 639, 640, 643, 644, 652, 657, 658, 665, 677, 679];


% HmAbs_all=[HmAbs_all  HmAbs_Bailey2017 HmAbs_Gopal2017 HmAbs_Keck2019 HmAb_Pierce2016 HmAbs_new];

HmAbs_all=[HmAbs_all  HmAbs_Gopal2017 HmAbs_Keck2019 HmAb_Pierce2016 HmAbs_new];

for kk = 1:length(HmAbs_all)
    data_mean_HmAb_E1E2{kk} = mean_escape_time_E1E2(HmAbs_all{kk}-191); 
    [min_data_mean_HmAb_E1E2(kk),I] =min(data_mean_HmAb_E1E2{kk});
    r = HmAbs_all{kk};
    res_E1E2(kk) = r(I);
end
load('escape_time_1a_500.mat', 'mean_escape_time')
for kk = 1:length(HmAbs_all)
    if all(HmAbs_all{kk}-383>0)
    
        data_mean_HmAb_E2{kk} = mean_escape_time(HmAbs_all{kk}-383); 
        [min_data_mean_HmAb_E2(kk),I] =min(data_mean_HmAb_E2{kk});
         r = HmAbs_all{kk};
         res_E2(kk) = r(I);
    else
        min_data_mean_HmAb_E2(kk) =0;
        res_E2(kk)=0;
    end
end


FIG=figure;

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)

% names = {'HEPC3',...
%     'HEPC43','HEPC74','HEPC46','HEPC50','HEPC98','HCV1','AR1A','AR1B','AR2A','AR3A','AR3B',...
%     'AR3C','AR3D','CBH-5','212.1.1',' 212.10','212.15','212.25','CBH-4D','CBH-4G','CBH-4B','CBH-20','CBH-21','CBH-22',...
%     'HC-1','HC-11',' A27','CBH-23','CBH-7','HC84-20','HC84-24','HC84-26',...
%     'HC33-1','HC33-4','IGH526','AR4A','AR5A'};
names = {'HCV1','AR1A','AR1B','AR2A','AR3A','AR3B',...
    'AR3C','AR3D','CBH-5','212.1.1',' 212.10','212.15','212.25','CBH-4D','CBH-4G','CBH-4B','CBH-20','CBH-21','CBH-22',...
    'HC-1','HC-11',' A27','CBH-23','CBH-7','HC84-20','HC84-24','HC84-26',...
    'HC33-1','HC33-4','IGH526','AR4A','AR5A'};
marker_size=30;
hold on
dx = 1;
dy =1*-0.015;


not_E2=[];
not_E1E2=[];
bad_antibody = [];
for kk = 1:length(HmAbs_all)
 if min_data_mean_HmAb_E2(kk)<100 && min_data_mean_HmAb_E1E2(kk)<100 && min_data_mean_HmAb_E2(kk)>0
        scatter(  min_data_mean_HmAb_E1E2(kk), min_data_mean_HmAb_E2(kk),marker_size,blue,'filled','MarkerFaceAlpha',0.6);
        hold on;
 end
  if min_data_mean_HmAb_E2(kk)>100 && min_data_mean_HmAb_E1E2(kk)<100
        scatter( min_data_mean_HmAb_E1E2(kk), min_data_mean_HmAb_E2(kk), marker_size,orange,'filled','MarkerFaceAlpha',0.6);
       
        hold on;
         text(  min_data_mean_HmAb_E1E2(kk)+dy',min_data_mean_HmAb_E2(kk)+dx', names(kk))
         hold on
            
         not_E2  = [not_E2 ; res_E2(kk)];
         not_E1E2  = [not_E1E2 ; res_E1E2(kk)];
         bad_antibody = [bad_antibody ; kk];
  end
   if min_data_mean_HmAb_E2(kk)>100 &&  min_data_mean_HmAb_E1E2(kk)>100
        scatter( min_data_mean_HmAb_E1E2(kk),min_data_mean_HmAb_E2(kk),  marker_size,red,'filled','MarkerFaceAlpha',0.6);
       
        hold on;
         text( min_data_mean_HmAb_E1E2(kk)+dy', min_data_mean_HmAb_E2(kk)+dx', names(kk))
         hold on
         
   end
      if min_data_mean_HmAb_E2(kk)==0
        scatter(min_data_mean_HmAb_E1E2(kk),  min_data_mean_HmAb_E2(kk), marker_size,[0.6 0.6 0.6 ],'filled','MarkerFaceAlpha',0.6);
       
        hold on;
         text( min_data_mean_HmAb_E1E2(kk)+dy',min_data_mean_HmAb_E2(kk)+dx',  names(kk))
         hold on
      end
end 

for i =1:length(bad_antibody)
    time = mean_escape_time_E1E2(HmAbs_all{bad_antibody(i)}-191);
    r = HmAbs_all{bad_antibody(i)};
    bad_residue{i} = r(time<100);
    
    
end

xlim([0 400]);
ylim([0 400]);

hold on;



 
 plot([0 400],[100 100],'--','Color','k');
 hold on;
 plot([96 96],[0 400],'--','Color','k');
 hold on;
 xlabel('Minimum escape time predicted by the JM')
 ylabel('Minimum escape time predicted by the E2-only model')


set(gca,'TickDir','out'); 
set(get(gca,'XLabel'),'FontSize',8,'Vertical','top');
% set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',8,'Vertical','middle');
grid off
FIG.Name = 'protein_length';
FIG.Units = 'centimeters';
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.08, 0.5, 0]);
set(get(gca,'XLabel'), 'Units', 'Normalized', 'Position', [0.5, -0.08, 0]);
FIG.Units = 'centimeters';
FIG.Name = 'RBall_1a';
set(gcf,'Position',[5 5  20 15]);
% ylabel('Escape time threshold, \tau')
set(gca,'Position',[ .09 .1  .87 .87]);  %调整 XLABLE和YLABLE不会被切掉
% print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');
% print(['C:\Users\27909\Desktop\1' ],'-dpng','-r600');
