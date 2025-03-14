%% Time-domain analysis

%% It's always good to start with a clean sheet
clear,clc,close all

%% Initialize the output matrix
E1_H = deal(zeros(28,1626)); % high-prob 
E2_H = deal(zeros(28,1626));
E1_L = deal(zeros(28,1626)); % low-prob
E2_L = deal(zeros(28,1626));
E1_A = deal(zeros(28,1626)); % alpha
E2_A = deal(zeros(28,1626));

%% Get all the data file names
E1 = 'G:\01\eegdata\proc_data\TF\E1'; % size(tf_pow): 1 59 13 1876
E2 = 'G:\01\eegdata\proc_data\TF\E2';

readdir = E2; %%%%%%%%%%%% E1:baseline condition; E2: learning condition %%%%%%%%%%%%%%%%
sublist=dir(readdir);
sublist={sublist.name};
sublist=sublist(3:30);

%% Analyse
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
    
    time = [0 3250]; % whole_epoch[0 3250], prestimulus[0 1250],poststimulus[1250 3250]
    t2plot=dsearchn(dim.times',time')';
    
    chan = {'O1','PO7','PO3','O2','PO8','PO4'};
    chan2plot=[];
    for ch=1:length(chan)
        chan2plot(ch) = find(strcmpi({dim.chans.labels},chan(ch)));
    end
    
    dat2plot = squeeze(mean(tfdata(:,chan2plot,f2plot,t2plot(1):t2plot(2)),2)); % channel average
      
    %% Save the data
    if strcmp(readdir,E1) 
        E1_H(subno,:) = dat2plot; 
    else
        E2_H(subno,:) = dat2plot;
    end
end   

%% low-prob dist-loc
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
    
    time = [0 3250]; % whole_epoch[0 3250], prestimulus[0 1250],poststimulus[1250 3250]
    t2plot=dsearchn(dim.times',time')';
    
    chan = {'O1','PO7','PO3','O2','PO8','PO4'};
    chan2plot=[];
    for ch=1:length(chan)
        chan2plot(ch) = find(strcmpi({dim.chans.labels},chan(ch)));
    end
    
    dat2plot = squeeze(mean(mean(tfdata(:,chan2plot,f2plot,t2plot(1):t2plot(2)),2),3)); % channel average and freq average
      
    %% Save the data
    if strcmp(readdir,E1) 
        E1_L(subno,:) = dat2plot; 
    else
        E2_L(subno,:) = dat2plot;
    end
end   

%% alpha band
for subno = 1:length(sublist)
    %% Load data
    dname = sublist{subno};
    fprintf('Loading subject %s to plot ...\n',dname);
    load([readdir filesep dname],'tf_pow','dim')
    tfdata = tf_pow; 
    
    %% Parameters for time_domain analysis
    freq = [9 12]; 
    f2plot=dsearchn(dim.freqs',freq')';
    
    time = [0 3250]; % whole_epoch[0 3250], prestimulus[0 1250],poststimulus[1250 3250]
    t2plot=dsearchn(dim.times',time')';
    chan = {'O1','PO7','PO3','O2','PO8','PO4'};
    chan2plot=[];
    for ch=1:length(chan)
        chan2plot(ch) = find(strcmpi({dim.chans.labels},chan(ch)));
    end
    
    dat2plot = squeeze(mean(mean(tfdata(:,chan2plot,f2plot(1):f2plot(2),t2plot(1):t2plot(2)),2),3)); % channel average,and then freq average
      
    %% Save the data
    if strcmp(readdir,E1) 
        E1_A(subno,:) = dat2plot; 
    else
        E2_A(subno,:) = dat2plot;
    end
end   

%% alpha, high, low的条件间差异
dif_H = E2_H - E1_H;
dif_L = E2_L - E1_L;
dif_A = E2_A-E1_A;

%% Save the output matrix
savedir = 'E:\01\eegcode\plot\dat2plot'; 
if ~exist(savedir,'dir')
    mkdir(savedir)
end
filename =  'time_domain_output.mat';
outputfilename = [savedir filesep filename];
save(outputfilename,'E1_H','E2_H','E1_L','E2_L','E1_A','E2_A', 'dif_H', 'dif_L', 'dif_A')


%% %%%%%%%%%%%%%%%%%%% Now ready to do the permutation test %%%%%%%%%%%%%%%%%%%%%%%%%
% Load data to plot
readdir = 'G:\01\eegcode\plot\dat2plot'; 
filename =  'time_domain_output.mat';
outputfilename = [readdir filesep filename];
load(outputfilename)


%% alpha, high, low，control frequencies 的条件间差异
% alpha, high, low 条件间差异
% high, low 条件间差异的比较
E1_dat = {E1_A, E1_H, E1_L, dif_H};
E2_dat = {E2_A, E2_H, E2_L, dif_L};
dlabel = {'Alpha','High-prob','Low-prob','High vs. Low'}; 
           
pval = 0.05;
nperm = 1000;
test_statistic = 'sum';
ntests = 1;
tail = 2;
conn = 8;
%%
sig_time = deal(zeros(4,1626));  % 保存显著的时间点
for pi = 1:4
  

    X1 = E1_dat{pi};
    X2 = E2_dat{pi};
   

 
    X = X2 - X1;

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
   disp -----------------------------------------------------------------------------
end


%% Save the output matrix
savedir = 'G:\01\eegcode\plot\dat2plot'; 
if ~exist(savedir,'dir')
    mkdir(savedir)
end
filename =  'sig_time_output_2.mat';
outputfilename = [savedir filesep filename];
save(outputfilename,'sig_time','dlabel','E1_dat', 'E2_dat','tftime')
