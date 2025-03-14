%% behavioral_plot
% 基线和学习条件分别做单因素方差分析的柱状图

%% It's always good to start with a clean sheet
clear,clc,close all

%% Load data to plot
readdir = 'E:\SSVEP\eegcode\plot\dat2plot'; 
filename =  'behavioral_plot.mat';
outputfilename = [readdir filesep filename];
load(outputfilename)
color_list = {[26,26,26]/355,[255,68,0]/255}; % baseline condition, learning condition

%% E1
x_list =[1,2];
y_list = [mean(E1_RT_noDis),mean(E1_RT_Dis)];
figure(); 
for i = 1:length(x_list)
    b1 = bar(x_list(i),y_list(i),0.5,'FaceColor',[1 1 1],'EdgeColor',color_list{1},'LineWidth',2.5);hold on
    if i == 2 
        hatchfill2(b1,'cross','HatchAngle',45,'HatchDensity',20,'HatchColor',color_list{1},'HatchLineWidth',1.5);       
    end
end
axis([0.4 2.6 700 1500]);
set(gca, 'xTick', [1, 2]);
set(gca,'XTickLabel',''); 
set(gca, 'yTick', [700:200:1500]);
set(gca,'YTickLabel',''); 

% plot error bar
nSubjects = 28;
e1 = std(E1_RT_noDis)/sqrt(nSubjects-1);
e2 = std(E1_RT_Dis)/sqrt(nSubjects-1);
E = [e1,e2];
errorbar(1, y_list(1),E(1),'linewidth',3.5,'LineStyle','none','Color',color_list{1});
errorbar(2, y_list(2),E(2),'linewidth',3.5,'LineStyle','none','Color',color_list{1});

% 依次画每个被试的点并相连
for i = 1:nSubjects
    x_list =[1,2];
    per_mean = [E1_RT_noDis(i),E1_RT_Dis(i)]; 
    Meandot_1 = scatter(x_list(1),per_mean(1),20, [0.5,0.5,0.5],'fill');hold on
    Meandot_2 = scatter(x_list(2),per_mean(2),20, [0.5,0.5,0.5],'fill');hold on
    L = line(x_list,per_mean);hold on
    L.LineWidth = 1.1;
    L.Color = [0.5,0.5,0.5];
end

% 调整坐标轴等
set(gca,'Box','off')
set(gca,'LineWidth',1.5,'xcolor','k','ycolor','k')
set(gca,'Layer','bottom')
set(gca,'tickdir','out') 



%% E2
x_list =[1,2,3];
y_list = [mean(E2_RT_noDis),mean(E2_RT_H),mean(E2_RT_L)];

figure(); 
for i = 1:length(x_list)
    b1 = bar(x_list(i),y_list(i),0.5,'FaceColor',[1 1 1],'EdgeColor',color_list{2},'LineWidth',2.5);hold on
     if i == 2 
        hatchfill2(b1,'single','HatchAngle',45,'HatchDensity',35,'HatchColor',color_list{2},'HatchLineWidth',2,'HatchLineStyle',':');  
     end
     if i == 3
        hatchfill2(b1,'single','HatchAngle',-45,'HatchDensity',30,'HatchColor',color_list{2},'HatchLineWidth',2);         
     end
end
axis([0.4 3.6 700 1500]);    %修改坐标轴显示范围，[x-min  x-max  y-min  y-max]
set(gca, 'xTick', [1, 2, 3]);
set(gca,'XTickLabel','');
set(gca, 'yTick', [700:200:1700]);
set(gca,'YTickLabel','');


% plot error bar
nSubjects = 28;
e1 = std(E2_RT_noDis)/sqrt(nSubjects-1);
e2 = std(E2_RT_H)/sqrt(nSubjects-1);
e3 = std(E2_RT_L)/sqrt(nSubjects-1);
E = [e1,e2,e3];
errorbar(1, y_list(1),E(1),'linewidth',3.5,'LineStyle','none','Color',color_list{2});
errorbar(2, y_list(2),E(2),'linewidth',3.5,'LineStyle','none','Color',color_list{2});
errorbar(3, y_list(3),E(3),'linewidth',3.5,'LineStyle','none','Color',color_list{2});

% 依次画每个被试的点并相连
for i = 1:nSubjects
    x_list =[1,2,3];
    per_mean = [E2_RT_noDis(i),E2_RT_H(i),E2_RT_L(i)]; 
    Meandot_1 = scatter(x_list(1),per_mean(1),20, [0.5,0.5,0.5],'fill');hold on
    Meandot_2 = scatter(x_list(2),per_mean(2),20, [0.5,0.5,0.5],'fill');hold on
    Meandot_3 = scatter(x_list(3),per_mean(3),20, [0.5,0.5,0.5],'fill');hold on
    L = line(x_list,per_mean);hold on
    L.LineWidth = 1.1;
    L.Color = [0.5,0.5,0.5];
end

% 调整坐标轴等
set(gca,'Box','off')
set(gca,'LineWidth',1.5,'xcolor','k','ycolor','k')
set(gca,'Layer','bottom')
set(gca,'tickdir','out') 

%% Save the plot
% check output dir
savedir = 'E:\01\figure\summarize_2\fig';
if ~exist(savedir,'dir')
    mkdir(savedir)
end
% 设置保存参数
resolution = 600;  % 分辨率（dpi）
figure(1);  print( [savedir filesep 'E1_RTs.tif'], '-dtiff', ['-r' num2str(resolution)]); 
figure(2); print([savedir filesep 'E2_RTs.tif'], '-dtiff', ['-r' num2str(resolution)]); 