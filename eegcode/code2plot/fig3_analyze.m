%% Time-domain analysis
% 将学习条件拆分为前五个block ‘E2’ 和 后五个block 'E3'

%% It's always good to start with a clean sheet
clear,clc,close all

%% Initialize the output matrix
E1_H = deal(zeros(28,1876)); % high-prob 
E2_H = deal(zeros(28,1876));
E3_H = deal(zeros(28,1876));

E1_L = deal(zeros(28,1876)); % low-prob
E2_L = deal(zeros(28,1876));
E3_L = deal(zeros(28,1876));

E1_A = deal(zeros(28,1876)); % alpha
E2_A = deal(zeros(28,1876));
E3_A = deal(zeros(28,1876));

%% Get all the data file names
E1 = 'E:\01\eegcode\plot_2\TF\E1'; % size(tf_pow): 1 59 13 1876
E2 = 'E:\01\eegcode\plot_2\TF\E2';
E3 = 'E:\01\eegcode\plot_2\TF\E3';

readdir = E3; %%%%%%%%%%%% E1:baseline condition; E2: the former 5 blocks of  learning condition; E3: the latter 5 blocks of learning condition
sublist=dir(readdir);
sublist={sublist.name};
sublist=sublist(3:30);

%% Analyse
%high-prob dist-loc
 for subno = 1:length(sublist)
    %% Load data
    dname = sublist{subno};
    fprintf('Loading subject %s to plot ...\n',dname);
    load([readdir filesep dname])
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

      
    %% Save the data
    if strcmp(readdir,E1) 
        E1_H(subno,:) = dat2plot; 
    elseif strcmp(readdir,E2)
        E2_H(subno,:) = dat2plot;
    else
        E3_H(subno,:) = dat2plot;
    end
end   

%% low-prob dist-loc
for subno = 1:length(sublist)
    %% Load data
    dname = sublist{subno};
    fprintf('Loading subject %s to plot ...\n',dname);
    load([readdir filesep dname])
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
      
    %% Save the data
    if strcmp(readdir,E1) 
        E1_L(subno,:) = dat2plot; 
    elseif strcmp(readdir,E2)
        E2_L(subno,:) = dat2plot;
    else
        E3_L(subno,:) = dat2plot;
    end
end   

%% alpha band
for subno = 1:length(sublist)
    %% Load data
    dname = sublist{subno};
    fprintf('Loading subject %s to plot ...\n',dname);
    load([readdir filesep dname])
    tfdata = tf_pow; 
    
    %% Parameters for time_domain analysis
    freq = [8 12]; 
    f2plot=dsearchn(dim.freqs',freq')';
    
    chan = {'O1','PO7','PO3','O2','PO8','PO4'};
    chan2plot=[];
    for ch=1:length(chan)
        chan2plot(ch) = find(strcmpi({dim.chans.labels},chan(ch)));
    end
    
    dat2plot = squeeze(mean(mean(tfdata(:,chan2plot,f2plot(1):f2plot(2),:),2),3)); % channel average,and then freq average
    
    %% Save the data
    if strcmp(readdir,E1) 
        E1_A(subno,:) = dat2plot; 
    elseif strcmp(readdir,E2)
        E2_A(subno,:) = dat2plot;
    else
        E3_A(subno,:) = dat2plot;
    end
end   

%% Save the output matrix
savedir = 'E:\01\eegcode\summarize_plot\plot_2\dat2plot'; 
if ~exist(savedir,'dir')
    mkdir(savedir)
end
filename =  'time_domain_output.mat';
outputfilename = [savedir filesep filename];
save(outputfilename,'E1_H','E2_H','E3_H','E1_L','E2_L','E3_L','E1_A','E2_A','E3_A')


%% %%%%%%%%%%%%%%%%%%% Now ready to do the permutation test %%%%%%%%%%%%%%%%%%%%%%%%%
% Load data to plot
readdir = 'E:\01\eegcode\summarize_plot\plot_2\dat2plot'; 
filename =  'time_domain_output.mat';
outputfilename = [readdir filesep filename];
load(outputfilename)

%% 
dif_A_21 = E2_A-E1_A;  % 学习条件前5个block减去基线条件的alpha
dif_A_31 = E3_A-E1_A;  % 学习条件后5个block减去基线条件的alpha
dif_A = dif_A_31-dif_A_21; % 前、后5个block的alpha差异

dif_HL_21 = (E2_H - E1_H)-(E2_L - E1_L); % dif_H - dif_L, 学习条件前5个block,high minus low difference
dif_HL_31 = (E3_H - E1_H)-(E3_L - E1_L); 
dif_HL = dif_HL_31 - dif_HL_21; % 前、后5个block的dif_HL 差异

datset = {dif_A_21,dif_A_31,dif_A,dif_HL_21,dif_HL_31,dif_HL};
dlabel = {'Alpha_E2','Alpha_E3','Alpha_E2_vs_E3'...
    'HL_E2','HL_E3','HL_E2_vs_E3'}; 
           
pval = 0.05;
nperm = 1000;
test_statistic = 'sum';
ntests = 1;
tail = 2;
conn = 8;
%%
sig_time = deal(zeros(6,1626));  % 保存显著的时间点
for pi = 1:6
  
    % 要注意的是，数据保存了[-500 3250]的值，但只在 [0 3250] 窗口内做检验
    time = [0 3250]; 
    times = -500:2:3250;
    t2plot=dsearchn(times',time')';
    
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
       disp(a)
       disp(b)
   end
   disp --------------------------------------------------------
end


%% Save the output matrix
savedir = 'E:\01\figure\summarize_2\eegcode\plot_2\dat2plot';
if ~exist(savedir,'dir')
    mkdir(savedir)
end
filename =  'sigtime.mat';
outputfilename = [savedir filesep filename];
save(outputfilename,'sig_time','dlabel','datset','tftime')

