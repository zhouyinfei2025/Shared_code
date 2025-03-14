%% topoplot


%% It's always good to start with a clean sheet
clear,clc,close all

%% Get all the data file names
E1 = 'E:\01\eegdata\proc_data\TF\E1'; % size(tf_pow): 1 59 13 1876
E2 = 'E:\01\eegdata\proc_data\TF\E2';

readdir = E1; %%%%%%%%%%%% E1:baseline condition; E2: learning condition %%%%%%%%%%%%%%%%
sublist=dir(readdir);
sublist={sublist.name};
sublist=sublist(3:30);
%% 
E1_H = deal(zeros(28,56)); %28个被试前56个电极的数据，平均[0 1250]时间段
E2_H = deal(zeros(28,56));
E1_L = deal(zeros(28,56)); %28个被试前56个电极的数据，平均[0 1250]时间段
E2_L = deal(zeros(28,56));
E1_A = deal(zeros(28,56)); %28个被试前56个电极的数据，平均[0 1250]时间段
E2_A = deal(zeros(28,56));

%% Analyse
%high-prob dist-loc
 for subno = 1:length(sublist)
    %% Load data
    dname = sublist{subno};
    fprintf('Loading subject %s to plot ...\n',dname);
    load([readdir filesep dname],'tf_pow','dim')
    tfdata = tf_pow;   % size(tf_pow): 1 59 13 1876
   
    %% Parameters for time_domain analysis
  
    % flashed freqs to the high-prob dist_loc of subject 1-28 (corresponding to subno 1-28)
    freq = [5.45,7.5,2.4,4.29,5.45,7.5,2.4,4.29,5.45,7.5,2.4,4.29,5.45,7.5,2.4,4.29,5.45,7.5,2.4,4.29,5.45,7.5,2.4,4.29,5.45,7.5,2.4,4.29];
    f2plot=dsearchn(dim.freqs',freq(subno)')';
    
    time = [0 1250]; % whole_epoch[0 3250], prestimulus[0 1250],poststimulus[1250 3250]
    t2plot=dsearchn(dim.times',time')';
    
    dat2plot = squeeze(mean(tfdata(:,1:56,f2plot,t2plot(1):t2plot(2)),4)); % time average
      
    %% Save the data
    if strcmp(readdir,E1) 
        E1_H(subno,:) = dat2plot; 
    else
        E2_H(subno,:) = dat2plot;
    end
 end   

 %%  low-prob dist-loc
for subno = 1:length(sublist)
    %% Load data
    dname = sublist{subno};
    fprintf('Loading subject %s to plot ...\n',dname);
    load([readdir filesep dname],'tf_pow','dim')
    tfdata = tf_pow;  % size(tf_pow): 1 59 13 1876
   
    %% Parameters for time_domain analysis
    
    % flashed freqs to the low-prob dist_loc of subject 1-28 (corresponding to subno 1-28)
    freq = [2.4,4.29,7.5;2.4,4.29,5.45;4.29,5.45,7.5;2.4,5.45,7.5;2.4,4.29,7.5;2.4,4.29,5.45;4.29,5.45,7.5;2.4,5.45,7.5;2.4,4.29,7.5;2.4,4.29,5.45;4.29,5.45,7.5;2.4,5.45,7.5;2.4,4.29,7.5;2.4,4.29,5.45;4.29,5.45,7.5;2.4,5.45,7.5;2.4,4.29,7.5;2.4,4.29,5.45;4.29,5.45,7.5;2.4,5.45,7.5;2.4,4.29,7.5;2.4,4.29,5.45;4.29,5.45,7.5;2.4,5.45,7.5;2.4,4.29,7.5;2.4,4.29,5.45;4.29,5.45,7.5;2.4,5.45,7.5];
    f2plot=dsearchn(dim.freqs',freq(subno,:)')';
    
    time = [0 1250]; % whole_epoch[0 3250], prestimulus[0 1250],poststimulus[1250 3250]
    t2plot=dsearchn(dim.times',time')';
    
    dat2plot = squeeze(mean(mean(tfdata(:,1:56,f2plot(1):f2plot(2),t2plot(1):t2plot(2)),3),4)); % freq average and then time average
      
    %% Save the data
    if strcmp(readdir,E1) 
        E1_L(subno,:,:) = dat2plot; 
    else
        E2_L(subno,:,:) = dat2plot;
    end
end   
 
%%  Alpha
for subno = 1:length(sublist)
    %% Load data
    dname = sublist{subno};
    fprintf('Loading subject %s to plot ...\n',dname);
    load([readdir filesep dname],'tf_pow','dim')
    tfdata = tf_pow; 
   
    %% Parameters for time_domain analysis
    freq = [8 12]; 
    f2plot=dsearchn(dim.freqs',freq')';
    
    time = [0 1250]; % whole_epoch[0 3250], prestimulus[0 1250],poststimulus[1250 3250]
    t2plot=dsearchn(dim.times',time')';
    
    dat2plot = squeeze(mean(mean(tfdata(:,1:56,f2plot(1):f2plot(2),t2plot(1):t2plot(2)),3),4)); % freq average and then time average
      
    %% Save the data
    if strcmp(readdir,E1) 
        E1_A(subno,:,:) = dat2plot; 
    else
        E2_A(subno,:,:) = dat2plot;
    end
end

%% 计算条件间差异，平均被试
dif_H = E2_H-E1_H;
dif_L = E2_L -E1_L;
dif_A = E2_A-E1_A;
mean_dif_H = mean(dif_H); 
mean_dif_L = mean(dif_L);
mean_dif_A = mean(dif_A);

%% Save the output matrix
savedir = 'E:\01\eegcode\plot\dat2plot'; 
if ~exist(savedir,'dir')
    mkdir(savedir)
end
filename =  'topoplot.mat';
outputfilename = [savedir filesep filename];
save(outputfilename, 'dif_H', 'dif_L', 'dif_A','mean_dif_H','mean_dif_L','mean_dif_A','dim')

%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%% Now ready to plot %%%%%%%%%%%%%%%%%%%%%%%%%


%% Load data
readdir = 'E:\01\eegcode\plot\dat2plot'; 
filename =  'topoplot.mat';
outputfilename = [readdir filesep filename];
load(outputfilename)

%% 横排放三个，alpha, low, high; 且三个图共享一个colorbar
figure()
set(gcf,'Position',[100 100 480 240]);  % figure的位置
t=tiledlayout(1,3,'TileSpacing','Compact');  % TileSpacing可以调节图之间的留白
nexttile
topoplot(mean_dif_A,dim.chans(1:56),'style','map', 'maplimits', [-1 1]) %将颜色缩放到指定范围
nexttile
topoplot(mean_dif_L,dim.chans(1:56),'style','map', 'maplimits', [-1 1]) %将颜色缩放到指定范围
nexttile
topoplot(mean_dif_H,dim.chans(1:56),'style','map', 'maplimits', [-1 1]) %将颜色缩放到指定范围
colormap cool;

% t.TileSpacing = 'tight';%由宽至窄 'loose'|'compact'|'tight'|'none'
% t.Padding = 'tight';%由大面积填充到紧致填充 'loose'|'compact'|'tight'

cb = colorbar;
cb.Layout.Tile = 'east';
set(cb,'tickdir','out')  % 朝外
set(cb,'Box','off') 
set(cb,'LineWidth',1.5,'xcolor','k','ycolor','k')
set(cb,'YTick',[-1 0 1]); %色标值范围及显示间隔
set(cb,'TickLabels',[]); 
ax = gca;
axpos = ax.Position;
cb.Position(3) = 0.5*cb.Position(3);
ax.Position = axpos;

%% Save the plot
savedir = 'E:\01\figure\topoplot';
if ~exist(savedir,'dir')
    mkdir(savedir)
end
% 设置保存参数
resolution = 600;  % 分辨率（dpi）
figure(1);  print( [savedir filesep 'topoplot.tif'], '-dtiff', ['-r' num2str(resolution)]);
