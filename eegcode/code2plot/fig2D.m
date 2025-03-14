%% Time_domain_plot

%% It's always good to start with a clean sheet
clear,clc,close all

%% Load data to plot 
readdir = 'E:\01\eegcode\summarize_plot\plot\dat2plot'; 
filename =  'time_domain_output.mat';
outputfilename = [readdir filesep filename];
load(outputfilename)

%% low-high 差异线
% 准备数据
X = dif_L-dif_H;
tftime = -500:2:3250;
figure
set(legend,'AutoUpdate','off')
plot(tftime,squeeze(mean(X)),'Color','#8E91F5','LineWidth',1.2),hold on 
nSubjects = size(X,1);
L1= shadedErrorBar(tftime,squeeze(mean(X)),std(X)/sqrt(nSubjects-1),'k','1');hold on
L1.patch.FaceColor = '#8E91F5';
L1.patch.FaceAlpha = 0.15;  
L1.mainLine.Visible = 'off';
L1.edge(1).Visible = 'off';
L1.edge(2).Visible = 'off';
% 作图细节
set(gca,'Box','off') % 去掉边框线
set(gca,'XAxisLocation','origin');    %将x轴的位置设置在y=0处
set(gca,'YAxisLocation','origin');     %将y轴的位置设置在x=0处
set(gca,'xlim',[-500 3300]);  % x轴的范围
set(gca,'xTick',[-400,-200,0,200,400,600,800,1000,1200,1250,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200]);  
set(gca,'xticklabel',[])     %去掉x轴上的数字
set(gca,'TickDir','out') % 刻度朝外
set(gca,'LineWidth',1.5,'xcolor','k','ycolor','k')
set(gca,'ylim',[-1.1 0.55]);  % y轴的范围，考虑到向上留出图注的位置，向下留出差异线的位置
set(gca,'yTick',[-2 -1 0 1]); % y轴的刻度，取整数不要留小数点
set(gca,'yticklabel',[])   %set(gca,'ydir','reverse');

%% 再画显著差异线
% 准备数据  E:\01\figure\summarize_2\eegcode\plot\dat2plot\sig_time_output.mat
for pi = 4
    tftime = 0:2:3250;
    conn = 8;
    plotclust = bwconncomp(sig_time(pi,:),conn);
    for blob=1:plotclust.NumObjects
        L3 = plot([tftime(plotclust.PixelIdxList{blob}(1)) tftime(plotclust.PixelIdxList{blob}(end))],[-0.8 -0.8],'color','#8E91F5','linewidth',1.5);
        a = tftime(plotclust.PixelIdxList{blob}(1));
        b = tftime(plotclust.PixelIdxList{blob}(end));
        disp(a)
        disp(b)
    end 
end

%% Save the plot
savedir = 'E:\01\figure\summarize_2\fig';
if ~exist(savedir,'dir')
    mkdir(savedir)
end
% 设置保存参数
savedir = 'E:\01\figure\summarize_2\fig';
resolution = 600;  % 分辨率（dpi）
figure(1);  print( [savedir filesep 'time_domain_L-H.tif'], '-dtiff', ['-r' num2str(resolution)]); 


