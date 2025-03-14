%% Power correlation across participants
% In specific, correlation beteen Δalpha [learning-baseline] before the stimulus onset (i.e. 0-1250ms)
% and difference between tagging response to high- and low-prob.location [low minus high]
% 用 cluster-based permutation test 进行校正

%% It's always good to start with a clean sheet
clear,clc,close all

%% Load the data to analyse
readdir = 'E:\01\eegcode\plot\dat2plot'; 
filename =  'time_domain_output.mat';
outputfilename = [readdir filesep filename];
load(outputfilename,'dif_A','dif_H','dif_L')
dif_HL = dif_L - dif_H;  % 考虑到alpha是一个负值，使dif_HL也是负值
nSubj = size(dif_A,1);

%% Sliding window and correlation
win = 25; % sampling points, 滑动窗口包含25个采样点，长50ms（采样率=500Hz）
slide = 10; % sampling points, 移动步长包含10个采样点，长20ms

%% Prepare the data to analyze
% Initialize the output matrix
Alpha_pow = zeros(nSubj,61);  % 28*61matrix
HL_pow= zeros(nSubj,61); 

% Ready to analyze
m = 1;
for i = 1:slide:626-win  % 0:2:1250时间段包含626个采样点
    Alpha_pow(:,m) = mean(dif_A(:,i:i+win),2);
    HL_pow(:,m) = mean(dif_HL(:,i:i+win),2); 
    m = m + 1;
end

% 把61个对应的时间窗算出来
i = 1:slide:626-win;
tftime = 0:2:1250;
t_begin = tftime(i);
t_end = tftime(i+win);

%%  Set some parameters
nperm = 1000;
pval = 0.05;
cluster_pval = 0.05;
conn=8;
tail = 2;  
r_thresh = 0.375; % 双尾检验，p=0.05, n=28, 只有当r值 ≥ 0.375时，相关显著

%% real r-values
corr_matrix = zeros(size(Alpha_pow,2),size(HL_pow,2)); % Initialize the corr_matrix
pval_matrix = zeros(61,61); % Initialize the pval_matrix

for a = 1:size(Alpha_pow,2)
    for b = 1:size(HL_pow,2)
        [R,P] = corr(Alpha_pow(:,a),HL_pow(:,b),'type','Spearman');  % 斯皮尔曼相关系数不需要假设数据呈现线性关系，对于一些非线性关系也能较好地反映出相关性
        corr_matrix(a,b) = R;
        pval_matrix(a,b) = P;
    end    
end

%%
% 对 corr_matrix 矩阵中小于阈值的元素进行逻辑判断，生成一个由逻辑值（True或False）组成的矩阵 thresh_matrix
% 并将满足条件的元素对应位置设为True，不满足条件的元素对应位置设为False。 thresh_matrix 矩阵可以用来筛选出显著性水平内的元素。
thresh_matrix = corr_matrix;
thresh_matrix(abs(corr_matrix) < r_thresh) = 0; % 对 Z_corr_matrix 矩阵中小于阈值的元素进行逻辑判断，生成一个由逻辑值（True或False）组成的矩阵，并将满足条件的元素对应位置设为True，不满足条件的元素对应位置设为False。这个矩阵 thresh_matrix 可以用来筛选出显著性水平内的元素。

%% permutation
% initialize null hypothesis matrices
max_clust_info   = zeros(nperm,1);  % nperm = 1000    

fprintf('Performing %i permutations:\n',nperm);

for permi=1:nperm
    if mod(permi,100)==0, fprintf('..%i\n',permi); end
    %% fake r-values
    fake_corr_matrix = zeros(size(Alpha_pow,2),size(HL_pow,2)); % Initialize the fake_corr_matrix
    %fake_pval_matrix = zeros(61,61); % Initialize the fake_pval_matrix

    for a = 1: size(Alpha_pow,2)
        alpha_pow = Alpha_pow(:,a); % 提取出某个时间窗的 alpha_pow
        temp_alpha_pow = alpha_pow(randperm(nSubj),:); % 打乱 alpha_pow 28个被试的顺序       
        for b = 1: size(HL_pow,2)
            R = corr(temp_alpha_pow,HL_pow(:,b),'type','Spearman');  % 斯皮尔曼相关系数不需要假设数据呈现线性关系，对于一些非线性关系也能较好地反映出相关性
            fake_corr_matrix(a,b) = R;
            %fake_pval_matrix(a,b) = P;
        end    
    end
    
    fake_corr_matrix(abs(fake_corr_matrix) < r_thresh) = 0;  % 将小于阈限的值标记为0
    % 计算所有clusters，找出最大的cluster并记录下元素数量 
    clustinfo = bwconncomp(fake_corr_matrix, conn); 
    clust_sizes = cellfun(@numel, clustinfo.PixelIdxList); 
    if ~isempty(clust_sizes)
        max_clust_info(permi) = max(clust_sizes);
    else
        max_clust_info(permi) = 0;
    end
    
end
fprintf('..Done!\n');

%% 进行cluster-level 校正
clustinfo = bwconncomp(thresh_matrix,conn);
clust_sizes = cellfun(@numel, clustinfo.PixelIdxList); 
% 计算null分布的阈值 
clust_threshold = prctile(max_clust_info, 100-cluster_pval*100);
% 将小于阈值的cluster对应的元素值设为0 
for i_clust = 1:length(clust_sizes) 
    if clust_sizes(i_clust) < clust_threshold 
        thresh_matrix(clustinfo.PixelIdxList{i_clust}) = 0; 
    end
end

%% Save the data to plot
% check output dir
savedir =  'G:\01\eegcode\plot\dat2plot';
if ~exist(savedir,'dir')
    mkdir(savedir)
end

filename =  'pow_corr_leaveout8Hz.mat';
outputfilename = [savedir filesep filename];
save(outputfilename,'corr_matrix','thresh_matrix')


%% Now ready to plot

figure()
set(gcf,'Position',[20 200 360 270]);
imagesc(0:20:1200, 0:20:1200, corr_matrix'); axis xy; colormap(jet);caxis([-1,1]); hold on  % 横坐标为alpha，纵坐标为high-low
if any(logical(thresh_matrix),'all')
    contour(0:20:1200, 0:20:1200,logical(thresh_matrix'),1,'linecolor','k','LineWidth',1);
end
set(gca,'xtick',[200:200:1250]);
set(gca,'ytick',[200:200:1250]);
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
set(gca,'Box','off')
set(gca,'LineWidth',1.2,'xcolor','k','ycolor','k')
set(gca,'tickdir','out') %把刻度线朝外
cb = colorbar;
cb.Location = 'eastoutside';
set(cb,'tickdir','out')  % 朝外
set(cb,'Box','off') 
set(cb,'LineWidth',1.2,'xcolor','k','ycolor','k')
set(cb,'Ticks',[-1 0 1]); %色标值范围及显示间隔 
set(cb,'TickLabels',[]); 
set(cb,'FontSize', 10);
ax = gca;
axpos = ax.Position;
cb.Position(3) = 0.5*cb.Position(3);
ax.Position = axpos;

%% Save the plot
savedir = 'G:\01\figure\summarize_2\fig';
if ~exist(savedir,'dir')
    mkdir(savedir)
end
% 设置保存参数
resolution = 600;  % 分辨率（dpi）
figure(1);  print( [savedir filesep 'pow_corr_leaveout8Hz.tif'], '-dtiff', ['-r' num2str(resolution)]); 
