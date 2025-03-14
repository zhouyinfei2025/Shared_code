%% figS1_new3
% 基线条件不区分高、低概率 (baseline)
% high, low in the learning condition

%% It's always good to start with a clean sheet
clear,clc,close all

%% Average across flickering frequencies in the baseline condition

% Initialize the output matrix
E1_flicker = deal(zeros(28,1876)); 

% Get all the data file names
readdir =  'E:\01\eegdata\proc_data\TF\E1'; % baseline condition
sublist=dir(readdir);
sublist={sublist.name};
sublist=sublist(3:30);

% Analyse
 for subno = 1:length(sublist)
    %% Load data
    dname = sublist{subno};
    fprintf('Loading subject %s to plot ...\n',dname);
    load([readdir filesep dname],'tf_pow','dim')
    tfdata = tf_pow;
   
    %% Parameters for time_domain analysis    
    chan = {'O1','PO7','PO3','O2','PO8','PO4'};
    chan2plot=[];
    for ch=1:length(chan)
        chan2plot(ch) = find(strcmpi({dim.chans.labels},chan(ch)));
    end
    
    dat2plot = squeeze(mean(mean(tfdata(:,chan2plot,1:4,:),2),3)); % channel average and freq average
      
    E1_flicker(subno,:) = dat2plot; 

 end   
 
 %% high, low in the learning condition
 
 % Initialize the output matrix
E1_H = deal(zeros(28,1876)); % high-prob 
E1_L = deal(zeros(28,1876)); % low-prob


% Get all the data file names
readdir =  'E:\01\eegdata\proc_data\TF\E2'; % learning condition
sublist=dir(readdir);
sublist={sublist.name};
sublist=sublist(3:30);

% Analyse
%high-prob dist-loc
 for subno = 1:length(sublist)
    %% Load data
    dname = sublist{subno};
    fprintf('Loading subject %s to plot ...\n',dname);
    load([readdir filesep dname],'tf_pow','dim')
    tfdata = tf_pow;
   
    %% Parameters for time_domain analysis
    % flashed freqs to the high-prob dist_loc of subject 1-28 (corresponding to subno 1-28)
    freq = [5.45,7.5,2.4,4.29,5.45,7.5,2.4,4.29,5.45,7.5,2.4,4.29,5.45,7.5,2.4,4.29,5.45,7.5,2.4,4.29,5.45,7.5,2.4,4.29,5.45,7.5,2.4,4.29];
    f2plot=dsearchn(dim.freqs',freq(subno)')';
    
    chan = {'O1','PO7','PO3','O2','PO8','PO4'};
    chan2plot=[];
    for ch=1:length(chan)
        chan2plot(ch) = find(strcmpi({dim.chans.labels},chan(ch)));
    end
    
    dat2plot = squeeze(mean(tfdata(:,chan2plot,f2plot,:),2)); % channel average
      
    E1_H(subno,:) = dat2plot; 


end   

% low-prob dist-loc
for subno = 1:length(sublist)
    %% Load data
    dname = sublist{subno};
    fprintf('Loading subject %s to plot ...\n',dname);
    load([readdir filesep dname],'tf_pow','dim')
    tfdata = tf_pow; 
    
    %% Parameters for time_domain analysis
    
    % flashed freqs to the low-prob dist_loc of subject 1-28 (corresponding to subno 1-28)
    freq = [2.4,4.29,7.5;2.4,4.29,5.45;4.29,5.45,7.5;2.4,5.45,7.5;2.4,4.29,7.5;2.4,4.29,5.45;4.29,5.45,7.5;2.4,5.45,7.5;2.4,4.29,7.5;2.4,4.29,5.45;4.29,5.45,7.5;2.4,5.45,7.5;2.4,4.29,7.5;2.4,4.29,5.45;4.29,5.45,7.5;2.4,5.45,7.5;2.4,4.29,7.5;2.4,4.29,5.45;4.29,5.45,7.5;2.4,5.45,7.5;2.4,4.29,7.5;2.4,4.29,5.45;4.29,5.45,7.5;2.4,5.45,7.5;2.4,4.29,7.5;2.4,4.29,5.45;4.29,5.45,7.5;2.4,5.45,7.5];
    f2plot=dsearchn(dim.freqs',freq(subno,:)')';
    
    chan = {'O1','PO7','PO3','O2','PO8','PO4'};
    chan2plot=[];
    for ch=1:length(chan)
        chan2plot(ch) = find(strcmpi({dim.chans.labels},chan(ch)));
    end
    
    dat2plot = squeeze(mean(mean(tfdata(:,chan2plot,f2plot,:),2),3)); % channel average and freq average

    E1_L(subno,:) = dat2plot; 

end   
 
 %% permutation test
% parameters
pval = 0.05;
nperm = 1000;
test_statistic = 'sum';
ntests = 1;
tail = 2;
conn = 8;
datset = {E1_flicker,E1_H, E1_L, E1_H - E1_L, E1_flicker-E1_H, E1_flicker-E1_L};
time = [0 3250]; % whole_epoch[0 3250], prestimulus[0 1250],poststimulus[1250 3250]
t2plot=dsearchn(dim.times',time')';
sig_time = deal(zeros(6,1626));  % 保存显著的时间点

for pi = 1:6
  

%     X1 = E1_dat{pi};
%     X2 = E2_dat{pi};
%     
%     X = X2 - X1;

    X = datset{pi}(:,t2plot(1):t2plot(2));
    
    nSubjects = size(X,1);
    voxel_pval   = pval;
    cluster_pval = pval;
    % initialize null hypothesis matrices
    max_clust_info   = zeros(nperm,1);

    %% real t-values

    [~,p,~,tmp] = ttest(X);
    tmap = squeeze(tmp.tstat);
    p = squeeze(p);
    realmean = squeeze(mean(X));

    % uncorrected pixel-level threshold
    threshmean = realmean;
    tmapthresh = tmap;
    if tail == 2
        tmapthresh(abs(tmap)<tinv(1-voxel_pval/2,nSubjects-1))=0;
        threshmean(abs(tmap)<tinv(1-voxel_pval/2,nSubjects-1))=0;
    elseif strcmp(tail,'left')
        tmapthresh(tmap>-1.*tinv(1-voxel_pval,nSubjects-1))=0;
        threshmean(tmap>-1.*tinv(1-voxel_pval,nSubjects-1))=0;
    elseif strcmp(tail,'right')
        tmapthresh(tmap<tinv(1-voxel_pval,nSubjects-1))=0;
        threshmean(tmap<tinv(1-voxel_pval,nSubjects-1))=0;
    end
    %%
    fprintf('Performing %i permutations:\n',nperm);

    for permi=1:nperm

        if mod(permi,100)==0, fprintf('..%i\n',permi); end

        clabels = logical(sign(randn(nSubjects,1))+1);
        tlabels = logical(sign(randn(nSubjects,1))+1);
        flabels = logical(sign(randn(nSubjects,1))+1);

        temp_permute = X; % X is the data structure and is assumed to be a difference score
        temp_permute(clabels,:,:)=temp_permute(clabels,:,:)*-1;

        %% permuted t-maps
        [~,~,~,tmp] = ttest(squeeze(temp_permute));

        faketmap = squeeze(tmp.tstat);
        faketmap(abs(faketmap)<tinv(1-voxel_pval/2,nSubjects-1))=0;
        if tail == 2
            faketmap(abs(faketmap)<tinv(1-voxel_pval/2,nSubjects-1))=0;
        elseif strcmp(tail,'left')
            faketmap(faketmap>-1.*tinv(1-voxel_pval,nSubjects-1))=0;
        elseif strcmp(tail,'right')
            faketmap(faketmap<tinv(1-voxel_pval,nSubjects-1))=0;
        end

        % get number of elements in largest supra-threshold cluster
        clustinfo = bwconncomp(faketmap,conn);
        if strcmp(test_statistic,'count')
            max_clust_info(permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % the zero accounts for empty maps; % using cellfun here eliminates the need for a slower loop over cells
        elseif strcmp(test_statistic,'sum')
            tmp_clust_sum = zeros(1,clustinfo.NumObjects);
            for ii=1:clustinfo.NumObjects
                tmp_clust_sum(ii) = sum(abs(faketmap(clustinfo.PixelIdxList{ii})));
            end
            if  clustinfo.NumObjects>0, max_clust_info(permi) = max(tmp_clust_sum); end
        else
            error('Absent or incorrect test statistic input!');
        end

    end
    fprintf('..Done!\n');

    %% apply cluster-level corrected threshold

    % find islands and remove those smaller than cluster size threshold
    clustinfo = bwconncomp(tmapthresh,conn);
    if strcmp(test_statistic,'count')
        clust_info = cellfun(@numel,clustinfo.PixelIdxList); % the zero accounts for empty maps; % using cellfun here eliminates the need for a slower loop over cells
    elseif strcmp(test_statistic,'sum')
        clust_info = zeros(1,clustinfo.NumObjects);
        for ii=1:clustinfo.NumObjects
            clust_info(ii) = sum(abs(tmapthresh(clustinfo.PixelIdxList{ii})));
        end
    end
    clust_threshold = prctile(max_clust_info,100-cluster_pval*100);

    % identify clusters to remove
    whichclusters2remove = find(clust_info<clust_threshold);

    % compute p-n value for all clusters
    clust_pvals = zeros(1,length(clust_info));
    clust_act = clust_pvals;
    for cp=1:length(clust_info)
        clust_pvals(cp) = length(find(max_clust_info>clust_info(cp)))/nperm;
        clust_act(cp) = sum(tmapthresh(clustinfo.PixelIdxList{cp}));
    end

    % remove clusters
    for i=1:length(whichclusters2remove)
        tmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
    end
    
   %% 保存显著的时间点
   tftime = 0:2:3250;
   sig_time(pi,:) = tmapthresh;
   conn = 8;
   plotclust = bwconncomp(tmapthresh,conn);
   for blob=1:plotclust.NumObjects
       a = tftime(plotclust.PixelIdxList{blob}(1));
       b = tftime(plotclust.PixelIdxList{blob}(end));
       disp([a b])
   end
   disp -----------------------------------------------------------------------------
end

%% Now ready to plot
tftime = -500:2:3250;
datset = {E1_flicker, E1_H, E1_L};
color_list = {'#493131', '#B262F5','#F58AA4'};
figure()
set(gcf,'Position',[20 200 360 270]);
for pi = 1:3
    X = datset{pi};        
    set(legend,'AutoUpdate','off')
    plot(tftime,squeeze(mean(X)),'Color',color_list{pi},'LineWidth',1),hold on 
    nSubjects = size(X,1);
    L1= shadedErrorBar(tftime,squeeze(mean(X)),std(X)/sqrt(nSubjects-1),'k','1');hold on
    L1.patch.FaceColor = color_list{pi};
    L1.patch.FaceAlpha = 0.15;  
    L1.mainLine.Visible = 'off';
    L1.edge(1).Visible = 'off';
    L1.edge(2).Visible = 'off';  
end
%  作图细节
set(gca,'Box','off') % 去掉边框线
set(gca,'XAxisLocation','origin');    %将x轴的位置设置在y=0处
set(gca,'YAxisLocation','origin');     %将y轴的位置设置在x=0处
set(gca,'xlim',[-500 3300]);  % x轴的范围
set(gca,'xTick',[-400,-200,0,200,400,600,800,1000,1200,1250,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200]);  
set(gca,'xticklabel',[])     %去掉x轴上的数字
set(gca,'TickDir','out') % 刻度朝外
set(gca,'ylim',[-3.6 3.5]);  % y轴的范围，考虑到向上留出图注的位置，向下留出差异线的位置
set(gca,'yTick',[-3 -2 -1 0 1 2 3],'Fontsize',7,'Fontname','Arial'); % y轴的刻度，取整数不要留小数点
set(gca,'yticklabel',[])   %set(gca,'ydir','reverse');
set(gca,'LineWidth',1,'xcolor','k','ycolor','k')

% 差异线
L1b = line([438 900],[2.8 2.8],'color',color_list{1},'linewidth',1); 
L2 = line([408 1394],[2.65 2.65],'color',color_list{2},'linewidth',1); 
L3 = line([426 2398],[2.5 2.5],'color',color_list{3},'linewidth',1); 

% baseline vs. high
L12a = line([1036 1886],[1.8 1.8],'color',color_list{1},'linewidth',1); 
L12b = line([1036 1886],[1.83 1.83],'color',color_list{2},'linewidth',1); 

% baseline vs. low
L13a = line([428 2120],[1.6 1.6],'color',color_list{1},'linewidth',1); 
L13b = line([428 2120],[1.53 1.53],'color',color_list{3},'linewidth',1); 

%% Save the plot
savedir = 'E:\01\figure\summarize_2\fig';
if ~exist(savedir,'dir')
    mkdir(savedir)
end
% 设置保存参数
resolution = 600;  % 分辨率（dpi）
figure(1); print([savedir filesep 'figS1(3).tif'], '-dtiff', ['-r' num2str(resolution)]); 