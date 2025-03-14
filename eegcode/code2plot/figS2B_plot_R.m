%% plot_CTFs_E2
% 区分高、低概率干扰位置


%% It's always good to start with a clean sheet
clear, close all, warning('off','all'),clc


%% Get data

readdir = 'E:\01\eegdata\proc_data\TF_cond_new\E2\distractor_iem_allchans';
sublist = dir(fullfile(readdir,'*.mat'));
sublist={sublist.name};
nSubs = length(sublist);

for id = 1:nSubs
    dname = sublist{id};
    load([readdir filesep dname],'em')
    
    high_index = 1;
    low_index = [2,3,4];
    ctf_high(id,:,:) = squeeze(mean(mean(mean(em.C2_total_shifted(1,:,:,:,high_index,:),5),4),2));
    ctf_low(id,:,:)= squeeze(mean(mean(mean(em.C2_total_shifted(1,:,:,:,low_index,:),5),4),2));%倒数第二个维度是bin
    
end

%% Alpha CTF Slope

% settings 
nSamps = length(em.time); % # of sample points

% preallocation slope matrix
slopes_high = nan(nSubs,nSamps);
slopes_low = nan(nSubs,nSamps);

% calculate slope values for each subject across time
for id = 1:nSubs
    for samp = 1:nSamps
        dat = squeeze(ctf_high(id,samp,:)); % dat 存的是以4个位置为"中心"（原点）分别的CTF
        x = 1:3;
        d = [dat(1),mean([dat(2),dat(4)]),dat(3)];
        fit = polyfit(x,d,1);
        slopes_high(id,samp)= fit(1); 
    end
end

for id = 1:nSubs
    for samp = 1:nSamps
        dat = squeeze(ctf_low(id,samp,:)); % dat 存的是以4个位置为"中心"（原点）分别的CTF
        x = 1:3;
        d = [dat(1),mean([dat(2),dat(4)]),dat(3)];
        fit = polyfit(x,d,1);
        slopes_low(id,samp)= fit(1); 
    end
end

%% Permutation test

% parameters
pval = 0.05;
nperm = 1000;
test_statistic = 'sum';
ntests = 1;
tail = 2;
conn = 8;

datset = {slopes_high; slopes_low; slopes_high-slopes_low};


for pi = 1:size(datset,1)
  

%     X1 = E1_dat{pi};
%     X2 = E2_dat{pi};
%    
%     X = X2 - X1;

    X = datset{pi};

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
   tftime = -500:2:3250;
   conn = 8;
   plotclust = bwconncomp(tmapthresh,conn);
   for blob=1:plotclust.NumObjects
       a = tftime(plotclust.PixelIdxList{blob}(1));
       b = tftime(plotclust.PixelIdxList{blob}(end));
       disp([a b])
   end
   disp -----------------------------------------------------------------------------
end

%% Plot
tftime = -500:2:3250;
datset = {slopes_high,slopes_low};
color_list = {'#B262F5','#F58AA4'};

figure()
set(gcf,'Position',[20 200 360 270]);
for pi = 1:size(datset,2)
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
%set(gca,'xTick',[-500,-200,0,200,400,600,800,1000,1200,1250,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200]); 
set(gca,'xTick',[-500,0,500,1000,1250,1500,2000,2500,3000]);  
set(gca,'xticklabel',[])     %去掉x轴上的数字
set(gca,'TickDir','out') % 刻度朝外
%set(gca,'ylim',[-0.3 0.3]);  % y轴的范围，考虑到向上留出图注的位置，向下留出差异线的位置
set(gca,'ylim',[-0.1 0.2]);
set(gca,'yTick',[-0.1 0 0.1 0.2],'Fontsize',7,'Fontname','Arial'); % y轴的刻度，取整数不要留小数点
set(gca,'yticklabel',[])   %set(gca,'ydir','reverse');

% yl=get(gca,'ylim');
% line([1250 1250],[min(yl) max(yl)],'LineStyle',':','Color','black');

set(gca,'LineWidth',1.1,'xcolor','k','ycolor','k')
% xlabel('Time (ms)', 'FontName', 'Arial','FontSize',12);
% ylabel('CTF slope', 'FontName', 'Arial','FontSize',12);
%title('Learning ', 'FontName', 'Arial','FontSize',14)


%% Save the output matrix
savedir = 'E:\01\eegcode(new)\IEM_new\dat2plot'; 
if ~exist(savedir,'dir'), mkdir(savedir); end
filename =  'distractor_iem_E2_allchans.mat';
outputfilename = [savedir filesep filename];
save(outputfilename,'ctf_high','ctf_low','slopes_high','slopes_low')

%% Save the plot
savedir = 'E:\01\eegcode(new)\IEM_new\dat2plot\tif';
if ~exist(savedir,'dir'),mkdir(savedir);end
% 设置保存参数
resolution = 600;  % 分辨率（dpi）
figure(1); print([savedir filesep 'learning_allchans.tif'], '-dtiff', ['-r' num2str(resolution)]);
