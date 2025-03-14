%% time_domain_plot2

% 准备数据，superior组 和 inferior组 low-high difference
datset = {dif_L_s-dif_H_s,dif_L_i-dif_H_i}; 
color_list = {'#E68C7C','#35BFE6'};
tftime = -500:2:3250;

%% 先画两组 low-high difference
figure
for pi = 1:2    
    X = datset{pi};
    set(legend,'AutoUpdate','off')
    plot(tftime,squeeze(mean(X)),'Color',color_list{pi},'LineWidth',1.2),hold on 
    nSubjects = size(X,1);
    L1= shadedErrorBar(tftime,squeeze(mean(X)),std(X)/sqrt(nSubjects-1),'k','1');hold on
    L1.patch.FaceColor = color_list{pi};
    L1.patch.FaceAlpha = 0.15;  

    L1.mainLine.Visible = 'off';
    L1.edge(1).Visible = 'off';
    L1.edge(2).Visible = 'off';
    
end
% 作图细节 
set(gca,'Box','off')
set(gca,'XAxisLocation','origin');    %将x轴的位置设置在y=0处
set(gca,'YAxisLocation','origin');     %将y轴的位置设置在x=0处
set(gca,'xlim',[-500 3300]);
set(gca,'xTick',[-400,-200,0,200,400,600,800,1000,1200,1250,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200]);  
set(gca,'xticklabel',[])     %去掉x轴上的数字
set(gca,'TickDir','out')
set(gca,'LineWidth',1.2,'xcolor','k','ycolor','k')
set(gca,'ylim',[-1.1 0.55]); 
set(gca,'yTick',[-2 -1 0 1]);
set(gca,'yticklabel',[]);

%% 再画 superior组 low-high difference 的显著差异线
for pi = 3
    tftime=0:2:3250;
    conn = 8;
    plotclust = bwconncomp(sig_time(pi,:),conn);
    for blob=1:plotclust.NumObjects
        L2 = plot([tftime(plotclust.PixelIdxList{blob}(1)) tftime(plotclust.PixelIdxList{blob}(end))],[-0.9 -0.9],'color',color_list{1},'linewidth',1.5);
    end 
end

%% inferior组 low-high difference 之间无显著差异

%% superior组 和 inferior组 之间的比较, 无差异

%% Save the plot
savedir = 'E:\01\figure\summarize_2\fig';
if ~exist(savedir,'dir')
    mkdir(savedir)
end
% 设置保存参数
savedir = 'E:\01\figure\summarize_2\fig';
resolution = 600;  % 分辨率（dpi）
figure(1);  print( [savedir filesep 'group_L-H.tif'], '-dtiff', ['-r' num2str(resolution)]); 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 准备数据，superior组 和 inferior组 alpha

datset = {dif_A_s,dif_A_i};
color_list = {'#E68C7C','#35BFE6'};
tftime = -500:2:3250;

figure
for pi = 1:2    
    X = datset{pi};
    set(legend,'AutoUpdate','off')
    plot(tftime,squeeze(mean(X)),'Color',color_list{pi},'LineWidth',1.2),hold on 
    nSubjects = size(X,1);
    L1= shadedErrorBar(tftime,squeeze(mean(X)),std(X)/sqrt(nSubjects-1),'k','1');hold on
    L1.patch.FaceColor = color_list{pi};
    L1.patch.FaceAlpha = 0.15;  

    L1.mainLine.Visible = 'off';
    L1.edge(1).Visible = 'off';
    L1.edge(2).Visible = 'off';
    
end
% 作图细节  
set(gca,'Box','off')
set(gca,'XAxisLocation','origin');    %将x轴的位置设置在y=0处
set(gca,'YAxisLocation','origin');     %将y轴的位置设置在x=0处
set(gca,'xlim',[-500 3300]);
set(gca,'xTick',[-400,-200,0,200,400,600,800,1000,1200,1250,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200]);  
set(gca,'xticklabel',[])     %去掉x轴上的数字
set(gca,'TickDir','out')
set(gca,'LineWidth',1.2,'xcolor','k','ycolor','k')
set(gca,'ylim',[-2.1 1.1]);
set(gca,'yTick',[-2 -1 0 1]);
set(gca,'yticklabel',[]);

%% 再画两组 alpha 各自的差异线

for pi = 1:2
    tftime = 0:2:3250;
    conn = 8;
    plotclust = bwconncomp(sig_time(pi+6,:),conn);
    for blob=1:plotclust.NumObjects
       L2 = plot([tftime(plotclust.PixelIdxList{blob}(1)) tftime(plotclust.PixelIdxList{blob}(end))],[-1.8 -1.8],'color',color_list{pi},'linewidth',1.5);
       a = tftime(plotclust.PixelIdxList{blob}(1));
       b = tftime(plotclust.PixelIdxList{blob}(end));
       disp(a)
       disp(b)       
    end 
    disp -------------------------------------------------
end

%% Save the plot
savedir = 'E:\01\figure\summarize_2\fig';
if ~exist(savedir,'dir')
    mkdir(savedir)
end
% 设置保存参数
savedir = 'E:\01\figure\summarize_2\fig';
resolution = 600;  % 分辨率（dpi）
figure(2);  print( [savedir filesep 'group_alpha.tif'], '-dtiff', ['-r' num2str(resolution)]); 