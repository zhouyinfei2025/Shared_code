%% Time domain average of 8-12Hz for all chans seperately

%% It's always good to start with a clean sheet
clear,clc,close all

%% Get all the data file names
E1 = 'E:\01\eegdata\proc_data\TF\E1';  % size(tf_pow): 1 59 13 1876
E2 = 'E:\01\eegdata\proc_data\TF\E2';

readdir = E2; %%%%%%%%%%%% E1:baseline condition; E2: learning condition %%%%%%%%%%%%%%%%
sublist=dir(readdir);
sublist={sublist.name};
sublist=sublist(3:30);

%% Intialize dat2plot for all subjects
dat2plot_E1 = deal(zeros(28,59,1626)); %28个被试的数据
dat2plot_E2 = deal(zeros(28,59,1626));

%% Analyse
for subno = 1:length(sublist)
    %% Load data
    dname = sublist{subno};
    fprintf('Loading subject %s to plot ...\n',dname);
    load([readdir filesep dname],'tf_pow','dim')
    tfdata = tf_pow; 
   
    %% Parameters for time_domain analysis
    freq = [8 12]; 
    f2plot=dsearchn(dim.freqs',freq')';
    
    time = [0 3250]; % whole_epoch[0 3250], prestimulus[0 1250],poststimulus[1250 3250]
    t2plot=dsearchn(dim.times',time')';
    
    dat2plot = squeeze(mean(tfdata(:,:,f2plot(1):f2plot(2),t2plot(1):t2plot(2)),3)); % freq average
      
    %% Save the data
    if strcmp(readdir,E1) 
        dat2plot_E1(subno,:,:) = dat2plot; 
    else
        dat2plot_E2(subno,:,:) = dat2plot;
    end
end   


%% %%%%%%%%%%%%%%%%%%%%  Now ready to plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
ha = tight_subplot(6, 10, [0.04 0.02], [0.06 0.04], [0.03 0.03]);

% RSM_permutation
pval = 0.05;
nperm = 1000;
test_statistic = 'sum';
ntests = 1;
tail = 2;
conn = 8;

% All chans
chanlist = {dim.chans.labels};
    
% loop around chans
for chani = 1:57
  

    X1 = squeeze(dat2plot_E1(:,chani,:));
    X2 = squeeze(dat2plot_E2(:,chani,:));
   

 
    X = X1 - X2;

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
    %% Plot
    index = find(tmapthresh>0);
    %sig_time = dim.times(index(1));
    tftime = 0:2:3250;
    axes(ha(chani))
    set(legend,'AutoUpdate','off')
    plot(tftime,squeeze(mean(X1)),'Color','k');hold on
    plot(tftime,squeeze(mean(X2)),'Color','r');
    %shadedErrorBar(x,y,errBar,lineProps,transparent)
    shadedErrorBar(tftime,squeeze(mean(X1)),std(X1)/sqrt(nSubjects-1),'k','0.2');hold on
    shadedErrorBar(tftime,squeeze(mean(X2)),std(X2)/sqrt(nSubjects-1),'r','0.2');hold on
    yl=get(gca,'ylim');
    line([min(xlim) max(xlim)],[0 0],'LineStyle','--','Color','black','LineWidth',0.5)
    line([1250 1250],[min(ylim) max(ylim)],'LineStyle','--','Color','black');
    set(gca,'xlim',[0 3250]);
    set(gca,'Box','off')
    title(chanlist(chani))
    %xlabel('Time (ms)'), 
    %ylabel('Power（dB）')
    hold on
    plotclust = bwconncomp(tmapthresh,conn);
    for blob=1:plotclust.NumObjects
        plot([tftime(plotclust.PixelIdxList{blob}(1)) tftime(plotclust.PixelIdxList{blob}(end))],[min(yl)+sum(abs(yl))/20 min(yl)+sum(abs(yl))/20],'color','black','linewidth',1.5);
    end
end
set(ha(1:50),'XTickLabel',''); 

%% Save the plot
savedir = 'E:\01\figure\topoplot';
if ~exist(savedir,'dir')
    mkdir(savedir)
end
saveas(figure(1),[savedir filesep 'alpha.png']) 
