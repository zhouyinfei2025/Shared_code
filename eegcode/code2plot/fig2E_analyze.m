%% Time_domain_analyis
% 根据 ΔRT [low-high] 将被试分为 superior learners 和 inferior learners
% 根据 50% 的标准，每组14名被试
% 比较两组的 time domain results

% superior learners：subject 1,3,4,9,12,14,15,16,17,18,19,20,24,25  
% inferior learners: subject 2,5,6,7,8,10,11,13,21,22,23,26,27,28   
% subject和subno的顺序是统一的，避免给高、低概率位置匹配频率时出错
% 将tf中的E1分为'E1_superior'和'E1_inferior'，E2分为'E2_superior'和'E2_inferior'

%% %%%%%%%%%%%%%%%%%%%%%%%%% Superior learners %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subject 1,3,4,9,12,14,15,16,17,18,19,20,24,25 (n=14)

%% It's always good to start with a clean sheet
clear,clc,close all

%% Initialize the output matrix
E1_H_s = deal(zeros(14,1876)); 
E2_H_s = deal(zeros(14,1876)); 
E1_L_s = deal(zeros(14,1876));
E2_L_s = deal(zeros(14,1876));
E1_A_s = deal(zeros(14,1876));
E2_A_s = deal(zeros(14,1876));
E1_C_s = deal(zeros(14,1876));
E2_C_s = deal(zeros(14,1876));

%% Get all the data file names
E1 = 'E:\01\eegcode\plot_group\TF\E1\E1_superior'; % size(tf_pow): 1 59 13 1876
E2 = 'E:\01\eegcode\plot_group\TF\E2\E2_superior';

readdir = E2; %%%%%%%%%%%% E1:baseline condition; E2: learning condition %%%%%%%%%%%%%%%%
sublist=dir(readdir);
sublist={sublist.name};
sublist=sublist(3:16); % 14 subjects

%% Analyse
%high-prob dist-loc
 for subno = 1:length(sublist) % subno和subject的顺序是一致的
    %% Load data
    dname = sublist{subno};
    fprintf('Loading subject %s to plot ...\n',dname);
    load([readdir filesep dname],'tf_pow','dim')
    tfdata = tf_pow;
   
    %% Parameters for time_domain analysis
    % flashed freqs to the high-prob dist_loc of subject 1,3,4,9,12,14,15,16,17,18,19,20,24,25
    freq = [5.45,2.4,4.29,5.45,4.29,7.5,2.4,4.29,5.45,7.5,2.4,4.29,4.29,5.45];
    f2plot=dsearchn(dim.freqs',freq(subno)')';
    
    chan = {'O1','PO7','PO3','O2','PO8','PO4'};
    chan2plot=[];
    for ch=1:length(chan)
        chan2plot(ch) = find(strcmpi({dim.chans.labels},chan(ch)));
    end
    
    dat2plot = squeeze(mean(tfdata(:,chan2plot,f2plot,:),2)); % channel average
      
    %% Save the data
    if strcmp(readdir,E1) 
        E1_H_s(subno,:) = dat2plot; 
    else
        E2_H_s(subno,:) = dat2plot;
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
    
    % flashed freqs to the low-prob dist_loc of of subject 1,3,4,9,12,14,15,16,17,18,19,20,24,25
    freq = [2.4,4.29,7.5;4.29,5.45,7.5;2.4,5.45,7.5;2.4,4.29,7.5;2.4,5.45,7.5;2.4,4.29,5.45;4.29,5.45,7.5;2.4,5.45,7.5;2.4,4.29,7.5;2.4,4.29,5.45;4.29,5.45,7.5;2.4,5.45,7.5;2.4,5.45,7.5;2.4,4.29,7.5];
    f2plot=dsearchn(dim.freqs',freq(subno,:)')';
    
    chan = {'O1','PO7','PO3','O2','PO8','PO4'};
    chan2plot=[];
    for ch=1:length(chan)
        chan2plot(ch) = find(strcmpi({dim.chans.labels},chan(ch)));
    end
    
    dat2plot = squeeze(mean(mean(tfdata(:,chan2plot,f2plot,:),2),3)); % channel average and freq average
      
    %% Save the data
    if strcmp(readdir,E1) 
        E1_L_s(subno,:) = dat2plot; 
    else
        E2_L_s(subno,:) = dat2plot;
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
        E1_A_s(subno,:) = dat2plot; 
    else
        E2_A_s(subno,:) = dat2plot;
    end
end   

%% alpha, high, low，control frequencies 的条件间差异 
dif_H_s = E2_H_s - E1_H_s; % s: superior
dif_L_s = E2_L_s - E1_L_s;
dif_A_s = E2_A_s -E1_A_s;

%% Save the output matrix
savedir = 'E:\01\eegcode\summarize_plot\plot_group\dat2plot'; 
if ~exist(savedir,'dir')
    mkdir(savedir)
end
filename =  'superior_learners.mat';
outputfilename = [savedir filesep filename];
save(outputfilename,'dif_H_s','dif_L_s','dif_A_s')


%% %%%%%%%%%%%%%%%%%%%%%%%% inferior learners %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subject 2,5,6,7,8,10,11,13,21,22,23,26,27,28 (n=14)


%% It's always good to start with a clean sheet
clear,clc,close all

%% Initialize the output matrix
E1_H_i = deal(zeros(14,1876)); 
E2_H_i = deal(zeros(14,1876)); 
E1_L_i = deal(zeros(14,1876));
E2_L_i = deal(zeros(14,1876));
E1_A_i = deal(zeros(14,1876));
E2_A_i = deal(zeros(14,1876));
E1_C_i = deal(zeros(14,1876));
E2_C_i = deal(zeros(14,1876));

%% Get all the data file names
E1 = 'E:\01\eegcode\plot_group\TF\E1\E1_inferior'; % size(tf_pow): 1 59 13 1876
E2 = 'E:\01\eegcode\plot_group\TF\E2\E2_inferior';

readdir = E2; %%%%%%%%%%%% E1:baseline condition; E2: learning condition %%%%%%%%%%%%%%%%
sublist=dir(readdir);
sublist={sublist.name};
sublist=sublist(3:16); % 14 subjects

%% Analyse
%high-prob dist-loc
 for subno = 1:length(sublist) % subno和subject的顺序是一致的
    %% Load data
    dname = sublist{subno};
    fprintf('Loading subject %s to plot ...\n',dname);
    load([readdir filesep dname],'tf_pow','dim')
    tfdata = tf_pow;
   
    %% Parameters for time_domain analysis
    % flashed freqs to the high-prob dist_loc of subject 2,5,6,7,8,10,11,13,21,22,23,26,27,28 
    freq = [7.5,5.45,7.5,2.4,4.29,7.5,2.4,5.45,5.45,7.5,2.4,7.5,2.4,4.29];
    f2plot=dsearchn(dim.freqs',freq(subno)')';
    
    chan = {'O1','PO7','PO3','O2','PO8','PO4'};
    chan2plot=[];
    for ch=1:length(chan)
        chan2plot(ch) = find(strcmpi({dim.chans.labels},chan(ch)));
    end
    
    dat2plot = squeeze(mean(tfdata(:,chan2plot,f2plot,:),2)); % channel average      

    %% Save the data
    if strcmp(readdir,E1) 
        E1_H_i(subno,:) = dat2plot; 
    else
        E2_H_i(subno,:) = dat2plot;
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
    
    % flashed freqs to the low-prob dist_loc of of subject 2,5,6,7,8,10,11,13,21,22,23,26,27,28 
    freq = [2.4,4.29,5.45;2.4,4.29,7.5;2.4,4.29,5.45;4.29,5.45,7.5;2.4,5.45,7.5;2.4,4.29,5.45;4.29,5.45,7.5;2.4,4.29,7.5;2.4,4.29,7.5;2.4,4.29,5.45;4.29,5.45,7.5;2.4,4.29,5.45;4.29,5.45,7.5;2.4,5.45,7.5];
    f2plot=dsearchn(dim.freqs',freq(subno,:)')';
    
    chan = {'O1','PO7','PO3','O2','PO8','PO4'};
    chan2plot=[];
    for ch=1:length(chan)
        chan2plot(ch) = find(strcmpi({dim.chans.labels},chan(ch)));
    end
    
    dat2plot = squeeze(mean(mean(tfdata(:,chan2plot,f2plot,:),2),3)); % channel average and freq average

      
    %% Save the data
    if strcmp(readdir,E1) 
        E1_L_i(subno,:) = dat2plot; 
    else
        E2_L_i(subno,:) = dat2plot;
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
        E1_A_i(subno,:) = dat2plot; 
    else
        E2_A_i(subno,:) = dat2plot;
    end
end   

%% alpha, high, low，control frequencies 的条件间差异
dif_H_i = E2_H_i - E1_H_i;  % i: inferior
dif_L_i = E2_L_i - E1_L_i;
dif_A_i = E2_A_i-E1_A_i;

%% Save the output matrix
savedir = 'E:\01\eegcode\summarize_plot\plot_group\dat2plot'; 
if ~exist(savedir,'dir')
    mkdir(savedir)
end
filename =  'inferior_learners.mat';
outputfilename = [savedir filesep filename];
save(outputfilename,'dif_H_i','dif_L_i','dif_A_i')

%%  %%%%%%%%%%%%%%%%%% permutation test %%%%%%%%%%%%%%%%%%%%%%%%%%
% superior 组的 high vs. low
% inferior 组的 high vs. low
% superior 和 inferior 组的 alpha

%% It's always good to start with a clean sheet
clear,clc,close all

%% Load data to analyse
readdir = 'E:\01\eegcode\plot_group\dat2plot';
filename_1 =  'superior_learners.mat';
filename_2 =  'inferior_learners.mat';
outputfilename_1 = [readdir filesep filename_1];
outputfilename_2 = [readdir filesep filename_2];
load(outputfilename_1)
load(outputfilename_2)

%% Ready to analyse

E1_dat = {E1_H_s; E1_L_s; dif_H_s; E1_H_i; E1_L_i;dif_H_i; E1_A_s; E1_A_i; dif_A_s};
E2_dat = {E2_H_s; E2_L_s; dif_L_s; E2_H_i; E2_L_i;dif_L_i; E2_A_s; E2_A_i; dif_A_i};
dlabel = {'superior_high','superior low','superior_high_vs_low'... % superior 组 high, low 各自的差异线以及两者的差异
    'inferior_high','inferior_low','inferior_high_vs_low'... % inferior 组 high, low 各自的差异线以及两者的差异
    'superior_alpha','inferior_alpha','alpha_superior_vs_inferior'}; % superior 和 inferior 组 alpha 的差异以及alpha在两组的差异
           
pval = 0.05; 
nperm = 1000;
test_statistic = 'sum';
ntests = 1;
tail = 2;
conn = 8;
%%
sig_time = deal(zeros(9,1626));  % 保存显著的时间点
for pi = 1:9
  

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
disp -------------------------------------------------
end

%% Save the output matrix
savedir = 'E:\01\eegcode\plot_group\dat2plot'; 
if ~exist(savedir,'dir')
    mkdir(savedir)
end
filename =  'sig_time_output.mat';
outputfilename = [savedir filesep filename];
save(outputfilename,'sig_time','dlabel','E1_dat', 'E2_dat','tftime')