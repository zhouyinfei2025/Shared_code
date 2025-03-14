%% It's always good to start with a clean sheet
clear,clc,close all

%% Get all the data file names
E1 = 'E:\01\eegdata\proc_data\TF\E1'; % size(tf_pow): 1 59 13 1876
E2 = 'E:\01\eegdata\proc_data\TF\E2';

readdir = E2; %%%%%%%%%%%% E1:baseline condition; E2: learning condition %%%%%%%%%%%%%%%%
sublist=dir(readdir);
sublist={sublist.name};
sublist=sublist(3:30);

%% Inatialize output matrix
dat2plot_E1 = zeros(28,12);%基线条件下，28名被试12个频率对应的Power值（平均电极和时间[0 3250]）
dat2plot_E2 = zeros(28,12);%学习条件下，28名被试12个频率对应的Power值（平均电极和时间[0 3250]）

%% Analysis dat2plot
for subno = 1:length(sublist)
    %% Load data
    dname = sublist{subno};
    fprintf('Loading subject %s to plot ...\n',dname);
    load([readdir filesep dname],'tf_pow','dim')
    tfdata = tf_pow; 
    
    %% Set parameters   
    time = [0 3250]; % whole_epoch[0 3250], prestimulus[0 1250],poststimulus[1250 3250]
    t2plot=dsearchn(dim.times',time')';
    
    chan = {'O1','PO7','PO3','O2','PO8','PO4'};
    chan2plot=[];
    for ch=1:length(chan)
        chan2plot(ch) = find(strcmpi({dim.chans.labels},chan(ch)));
    end
     
    dat2plot = squeeze(mean(mean(tfdata(:,chan2plot,:,t2plot(1):t2plot(2)),2),4)); % average chans & time
    
    if strcmp(readdir,E1)
        dat2plot_E1(subno,:) = dat2plot;
    else
        dat2plot_E2(subno,:) = dat2plot;
    end  
end

%% Save some output matrix
E1_average_tf = dat2plot_E1;
E2_average_tf = dat2plot_E2;

%% Save
% check output dir
savedir = 'G:\01\eegcode\plot\dat2plot'; 
if ~exist(savedir,'dir')
    mkdir(savedir)
end

filename =  'freq_matrix_output.mat';
outputfilename = [savedir filesep filename];
save(outputfilename,'E1_average_tf','E2_average_tf','dim')

%% frequency_domain_plot
% 用分组条形图表示两因素重复测量方差分析的结果
% 折现图表示频谱图

%% It's always good to start with a clean sheet
clear,clc,close all

%% Load data to plot 
readdir = 'G:\01\eegcode\plot\dat2plot'; 
filename =  'freq_matrix_output.mat';
outputfilename = [readdir filesep filename];
load(outputfilename,'E1_average_tf','E2_average_tf','dim')

%% 2 (learning, baseline) *3 (control, flashed, alpha) 的重复测量方差分析
E1_selected = mean(E1_average_tf(:,1:4),2);
E2_selected = mean(E2_average_tf(:,1:4),2);
E1_alpha = mean(E1_average_tf(:,5:9),2);
E2_alpha = mean(E2_average_tf(:,5:9),2);
E1_control = mean(E1_average_tf(:,10:12),2);
E2_control = mean(E2_average_tf(:,10:12),2);

%% 水平多的变量 (frequency type) 放在横轴，且实现颜色按横轴水平分类，一组两个柱状图之间留一点间隔
% 准备数据
x_list = [1,2.6,5,6.6,9,10.6];
y_list = [mean(E1_selected),mean(E2_selected),mean(E1_alpha),mean(E2_alpha),mean(E1_control),mean(E2_control)];
color_list = {[244,37,66]/255,[244,37,66]/255,[68,136,245]/255,[68,136,245]/255,[46,230,104]/255,[46,230,104]/255}; % flashed,alpha,control

figure()
for i = 1:length(x_list)
    if i == 1 || i == 3 || i ==5 % 基线条件 'FaceColor',[1 1 1]
        b1 = bar(x_list(i),y_list(i),1,'FaceColor',[1 1 1],'EdgeColor',color_list{i},'LineWidth',2.5);hold on
    elseif i == 2
        b1 = bar(x_list(i),y_list(i),1,'FaceColor',color_list{i},'FaceAlpha', 0.8,'EdgeColor',color_list{i},'LineWidth',2.5);hold on
    else  % 学习条件填充
        b1 = bar(x_list(i),y_list(i),1,'FaceColor',color_list{i},'FaceAlpha', 0.9,'EdgeColor',color_list{i},'LineWidth',2.5);hold on
    end
end

% 画误差条
nSubjects = 28;
e1 = errorbar(x_list(1),mean(E1_selected),std(E1_selected)/sqrt(nSubjects-1),'linewidth',3.5,'LineStyle','none','Color',color_list{1}); hold on
e2 = errorbar(x_list(2),mean(E2_selected),std(E2_selected)/sqrt(nSubjects-1),'linewidth',3.5,'LineStyle','none','Color',color_list{2}); hold on
e3 = errorbar(x_list(3),mean(E1_alpha),std(E1_alpha)/sqrt(nSubjects-1),'linewidth',3.5,'LineStyle','none','Color',color_list{3}); hold on
e4 = errorbar(x_list(4),mean(E2_alpha),std(E2_alpha)/sqrt(nSubjects-1),'linewidth',3.5,'LineStyle','none','Color',color_list{4}); hold on
e5 = errorbar(x_list(5),mean(E1_control),std(E1_control)/sqrt(nSubjects-1),'linewidth',3.5,'LineStyle','none','Color',color_list{5}); hold on
e6 = errorbar(x_list(6),mean(E2_control),std(E2_control)/sqrt(nSubjects-1),'linewidth',3.5,'LineStyle','none','Color',color_list{6}); hold on

% 依次画每个被试的点，且相同条件的被试相连
% selected
for i = 1:nSubjects
    x_list =[1,2.6];
    per_mean = [E1_selected(i),E2_selected(i)]; 
    Meandot_1 = scatter(x_list(1),per_mean(1),15, [0.5,0.5,0.5],'fill');hold on
    Meandot_2 = scatter(x_list(2),per_mean(2),15, [0.5,0.5,0.5],'fill');hold on
    L = line(x_list,per_mean);hold on
    L.LineWidth = 1;
    L.Color = [0.5,0.5,0.5];
end
% alpha
for i = 1:nSubjects
    x_list =[5,6.6];
    per_mean = [E1_alpha(i),E2_alpha(i)]; 
    Meandot_1 = scatter(x_list(1),per_mean(1),15, [0.5,0.5,0.5],'fill');hold on
    Meandot_2 = scatter(x_list(2),per_mean(2),15, [0.5,0.5,0.5],'fill');hold on
    L = line(x_list,per_mean);hold on
    L.LineWidth = 1;
    L.Color = [0.5,0.5,0.5];
end
% control
for i = 1:nSubjects
    x_list =[9,10.6];
    per_mean = [E1_control(i),E2_control(i)]; 
    Meandot_1 = scatter(x_list(1),per_mean(1),15, [0.5,0.5,0.5],'fill');hold on
    Meandot_2 = scatter(x_list(2),per_mean(2),15, [0.5,0.5,0.5],'fill');hold on
    L = line(x_list,per_mean);hold on
    L.LineWidth = 1;
    L.Color = [0.5,0.5,0.5];
end

% 修改横坐标名称、字体
set(gca,'xticklabel',[])     %去掉x轴上的标签
set(gca,'xtick',[])   % 去掉x轴上的刻度
% y 轴刻度  
ylim([-9,5]);
set(gca,'ytick',[-8,-4,0,4])  
% 作图细节
set(gca,'Box','off')
set(gca,'LineWidth',1.5,'xcolor','k','ycolor','k')
set(gca,'XAxisLocation','origin');    %将x轴的位置设置在y=0处
set(gca,'tickdir','out') %把刻度线朝外
set(gca,'yticklabel',[])     %去掉y轴上的数字
set(gca,'Layer','bottom')

%% Save the plot
% check output dir
savedir = 'E:\01\figure\summarize_2\fig';
if ~exist(savedir,'dir')
    mkdir(savedir)
end
% 设置保存参数
resolution = 600;  % 分辨率（dpi）
figure(1);  print( [savedir filesep 'two-way AVOVA.tif'], '-dtiff', ['-r' num2str(resolution)]); 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the tf_pow
figure()
set(gcf,'Position',[100 100 400 280]);  % [left bottom width height]
% tight_subplot(Nh, Nw, gap, marg_h, marg_w)
% Nh, Nw用法同subplot(row, col)表示行数和列数
% gap（如[0.01, 0.1]）表示子图之间垂直方向和水平方向的间隔
% marg_h表示的是全部子图到figure上下边界的距离，marg_w则表示的是全部子图到figure左右边界
ha = tight_subplot(2, 1, [0.04 0.04], [0.1 0.1], [0.1 0.1]);
set(ha(1),'position',[0.15 0.43 0.75 0.5])
set(ha(2),'position',[0.15 0.13 0.75 0.25])
color_list = {[26,26,26]/355,[255,68,0]/255};

axes(ha(1))
set(legend,'AutoUpdate','off')

E1_Line = plot(dim.freqs,mean(E1_average_tf),'k - o','LineWidth',1.1);hold on
E1_Line.Color = color_list{1};
E1_Line.MarkerSize = 3;
E1_Line.MarkerEdgeColor = color_list{1};
E1_Line.MarkerFaceColor = color_list{1};

E2_Line = plot(dim.freqs,mean(E2_average_tf),'r - d','LineWidth',1.1);hold on
E2_Line.Color = color_list{2};
E2_Line.MarkerSize = 3;
E2_Line.MarkerEdgeColor = color_list{2};
E2_Line.MarkerFaceColor = color_list{2};

%shadedErrorBar(x,y,errBar,lineProps,transparent)
nSubjects = 28;
L1=shadedErrorBar(dim.freqs,squeeze(mean(E1_average_tf)),std(E1_average_tf)/sqrt(nSubjects-1),'k',1.5);hold on
L2=shadedErrorBar(dim.freqs,squeeze(mean(E2_average_tf)),std(E2_average_tf)/sqrt(nSubjects-1),'r',1.5);hold on
L1.patch.FaceColor = color_list{1};
L2.patch.FaceColor = color_list{2};
L1.patch.FaceAlpha = 0.1;    % 透明度
L2.patch.FaceAlpha = 0.1;
L1.edge(1).Visible = 'off';
L1.edge(2).Visible = 'off';
L2.edge(1).Visible = 'off';
L2.edge(2).Visible = 'off'; 
set(ha(1),'XTickLabel',''); 
set(gca,'ylim',[-7.5 3.5])
set(gca,'yticklabel',[])  
set(gca,'Box','off')
set(gca,'XAxisLocation','origin');    %将x轴的位置设置在y=0处
set(gca,'YAxisLocation','origin');     %将y轴的位置设置在x=0处
set(gca,'TickDir','out')
set(gca,'LineWidth',1.5,'xcolor','k','ycolor','k')

%% plot the difference
color_dif =[0,176,179]/255;
axes(ha(2))
set(legend,'AutoUpdate','off')
dif_Line = plot(dim.freqs,mean(E2_average_tf-E1_average_tf),'b - s','LineWidth',1.5);hold on
dif_Line.Color = color_dif;
dif_Line.MarkerSize = 2.5;
dif_Line.MarkerEdgeColor = color_dif;
dif_Line.MarkerFaceColor = color_dif;

nSubjects = 28;
L3=shadedErrorBar(dim.freqs,squeeze(mean(E2_average_tf - E1_average_tf)),std(E2_average_tf - E1_average_tf)/sqrt(nSubjects-1),'b',1);hold on
L3.patch.FaceColor = color_dif;
L3.patch.FaceAlpha = 0.15;    % 透明度
L3.mainLine.Visible = 'off';
L3.edge(1).Visible = 'off';
L3.edge(2).Visible = 'off';
set(gca,'Box','off')
set(gca,'XAxisLocation','origin');    %将x轴的位置设置在y=0处
set(gca,'YAxisLocation','origin');     %将y轴的位置设置在x=0处
set(gca,'TickDir','out')
set(gca,'ylim',[-1.1 0.1])
set(gca,'yTick',[-1 0])
set(gca,'yticklabel',[]) 
set(gca,'xticklabel',[])  
set(gca,'LineWidth',1.5,'xcolor','k','ycolor','k')

%% Save the plot
savedir = 'E:\01\figure\summarize_2\fig';
if ~exist(savedir,'dir')
    mkdir(savedir)
end
% 设置保存参数
resolution = 600;  % 分辨率（dpi）
figure(2);  print( [savedir filesep 'frequency_domain.tif'], '-dtiff', ['-r' num2str(resolution)]); 
