%% distractor_iem_allchans
% 探究对干扰位置的空间选择性
% 基线条件平均四个位置
% 学习条件区分高、低概率位置
% 2025/2/22

%% It's always good to start with a clean sheet
clear, close all, warning('off','all'),clc

%% Get all the data file names
%[~, filepath]=uigetfile('*.mat'); 
readdir = 'E:\01\eegdata\proc_data\TF_cond_new\E1'; % E1: baseline condition; E2: learning condition
sublist = dir(fullfile(readdir,'*.mat'));
sublist={sublist.name};

%% Analysis
for subno = 1:length(sublist)  
    %% Load data
    clear train_data_all
    dname = sublist{subno};
    fprintf('Loading subject %s for analysis ...\n',dname);
    load([readdir filesep dname],'tf_pow','dim')
    
    %% Secelt chans for analysis
%     chan = {'O1','PO7','PO3','O2','PO8','PO4'};
%     chan2plot=[];
%     for ch=1:length(chan)
%         chan2plot(ch) = find(strcmpi({dim.chans.labels},chan(ch)));
%     end
%     
    for i = 1:size(z_pow,2)
        tf_pow{i} = z_pow{i}(1:57,:,:);
    end
    
    %% Determine random seed
    rng('default')
    rng shuffle;
    em.randSeed = rng; 
    
    %% Prepare for analysis
    for i = 1:4  % 得到每个位置的数据
        eval([ 'location_' num2str(i) '_data_fre =tf_pow{' num2str(i) '};']);
    end
    
    trian_events=[];
    train_data_all=[];
    
    for i = 1:4  % 各位置数据合并，放在train_data_all
        eval(['train_data_all = cat(3,train_data_all,location_' num2str(i ) '_data_fre);'])
    end
    
    for i = 1:4 %打上marker
        eval([ 'location_' num2str(i) '_marker = zeros(1,size(location_' num2str(i) '_data_fre,3)) + ' num2str(i) ';']);
    end
    
    for i = 1:4 % 生成数据 each_location_trial 汇总各个位置的trial数
        eval([ 'each_location_trial(' num2str(i) ') = size(location_' num2str(i) '_data_fre,3);']);
    end
    
    for i = 1:4 % 将marker合在一起
        eval(['trian_events = cat(2,trian_events,location_' num2str(i) '_marker);'])
    end
    
    %% Set parameters
    em.nChans = 4; % # of channels 通道   
    em.nBins = em.nChans; % # of stimulus bins 等于通道
    em.nIter = 10; % # of iterations 迭代次数
    em.nBlocks = 3; % # of blocks for cross-validation 用于交叉验证的组块
    em.frequencies = [8 12]; % frequency bands to analyze
    em.bands = {'Alpha'};
    em.nElectrodes = size(location_1_data_fre,1);
    em.time = dim.times; % time points of interest
    em.window = 2; % 时间窗口
    em.Fs = 500; %1s有250个周期，4ms一个周期
    em.posBin = trian_events;
    em.nTrials = length(trian_events);
    
    %% For brevity in analysis
    posBin = em.posBin;
    nChans = em.nChans;
    nBins = em.nBins;
    nIter = em.nIter;
    nBlocks = em.nBlocks;
    freqs = em.frequencies;
    times = em.time;
    nFreqs = size(em.frequencies,1); %读取em.frequency矩阵的第一列数量
    nElectrodes = em.nElectrodes; 
    nSamps = length(em.time);
    Fs = em.Fs;
    nTrials = em.nTrials;
    nTimes = length(times);
    
    %% Specify basis set 形成通道假设模型，数据按照一定规律对齐
    em.sinPower = 7; 
    em.x = linspace(0, 2*pi-2*pi/nBins, nBins);
    em.cCenters = linspace(0, 2*pi-2*pi/nChans, nChans);
    em.cCenters = rad2deg(em.cCenters);
    pred = sin(0.5*em.x).^em.sinPower; % hypothetical channel responses
    pred = wshift('1D',pred,2); % shift the initial basis function; 4-2, 6-3, 8-4, 10-5, 12-6 ...
    basisSet = nan(nChans,nBins);
    for c = 1:nChans
        basisSet(c,:) = wshift('1D',pred,-c+1); % generate circularly shifted basis functions，使得对角线为1
    end
    em.basisSet = basisSet; % save basis set to data structure
    
    %% Preallocate Matrices
    tf_evoked = nan(nFreqs,nIter,nSamps,nBlocks,nChans); tf_total = tf_evoked;
    C2_evoked = nan(nFreqs,nIter,nSamps,nBlocks,nBins,nChans); C2_total = C2_evoked; C2_total_shifted = C2_evoked;
    em.blocks = nan(nTrials,nIter);  % create em.block to save block assignments
    
    %% Ready to analyze
    % Loop around each frequency
    for f = 1:nFreqs
        tic % 开始计时
        fdata_total = train_data_all;
        % Loop around each iteration
        for iter = 1:nIter
            
            %--------------------------------------------------------------------------
            % Assign trials to blocks (such that trials per position are
            % equated within blocks)
            %--------------------------------------------------------------------------
            
            % preallocate arrays
            blocks = nan(size(posBin));
            shuffBlocks = nan(size(posBin));
            
            % count number of trials within each position bin
            clear binCnt
            for bin = 1:nBins
                binCnt(bin) = sum(posBin == bin);
            end
            
            minCnt = min(binCnt); % # of trials for position bin with fewest trials
            nPerBin = floor(minCnt/nBlocks); % max # of trials such that the # of trials for each bin can be equated within each block
            
            % shuffle trials
            shuffInd = randperm(nTrials)'; % create shuffle index
            shuffBin = posBin(shuffInd); % shuffle trial order
            
            % take the 1st nPerBin x nBlocks trials for each position bin.
            for bin = 1:nBins
                idx = find(shuffBin == bin); % get index for trials belonging to the current bin
                idx = idx(1:nPerBin*nBlocks); % drop excess trials
                x = repmat(1:nBlocks',nPerBin,1); shuffBlocks(idx) = x; % assign randomly order trials to blocks
            end
            
            % unshuffle block assignment
            blocks(shuffInd) = shuffBlocks;
            
            % save block assignment
            em.blocks(:,iter) = blocks; % block assignment
            em.nTrialsPerBlock = length(blocks(blocks == 1)); % # of trials per block
            
            %-------------------------------------------------------------------------
            
            % Average data for each position bin across blocks
            posBins = 1:nBins;
            blockDat_total = nan(nBins*nBlocks,nElectrodes,nSamps);  % averaged total data
            labels = nan(nBins*nBlocks,1);                           % bin labels for averaged data
            blockNum = nan(nBins*nBlocks,1);                         % block numbers for averaged data
            c = nan(nBins*nBlocks,nChans);                           % predicted channel responses for averaged data
            bCnt = 1;
            for ii = 1:nBins
                for iii = 1:nBlocks
                    %evoked_data = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,:),1))).^2;
    %               blockDat_total(bCnt,:,:) = squeeze(mean(fdata_total(:,:,posBin==posBins(1) & blocks==1),3)); %%跑一个才跑这个
    %               x = posBin==posBins(1) & blocks==1;
                 
                    blockDat_total(bCnt,:,:) = squeeze(mean(fdata_total(:,:,posBin==posBins(ii) & blocks==iii),3)); %找trial,trial数在第三位
                    labels(bCnt) = ii;
                    blockNum(bCnt) = iii;
                    c(bCnt,:) = basisSet(ii,:);
                    bCnt = bCnt+1;
                end
            end
            
            for t = 1:nSamps
                % grab data for timepoint t
                toi = ismember(times,times(t)-em.window/2:times(t)+em.window/2); 
                dt = squeeze(mean(blockDat_total(:,:,toi),3));  % total data
                
                % Do forward model
                for i=1:nBlocks % loop through blocks, holding each out as the test set
                    
                    trnl = labels(blockNum~=i); % training labels
                    tstl = labels(blockNum==i); % test labels
                    
                    %-----------------------------------------------------%
                    % Analysis on Evoked Power                            %
                    %-----------------------------------------------------%
    %                 B1 = de(blockNum~=i,:);    % training data
    %                 B2 = de(blockNum==i,:);    % test data
    %                 C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
    %                 W = C1\B1;          % estimate weight matrix
    %                 C2 = (W'\B2')';     % estimate channel responses
    %                 
    %                 C2_evoked(f,iter,t,i,:,:) = C2; % save the unshifted channel responses
                    
                    % shift eegs to common center
    %                 n2shift = ceil(size(C2,2)/2);
    %                 for ii=1:size(C2,1)
    %                     [~, shiftInd] = min(abs(posBins-tstl(ii)));
    %                     C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
    %                 end
    %                 
    %                 tf_evoked(f,iter,t,i,:) = mean(C2,1); % average shifted channel responses
                    
                    %-----------------------------------------------------%
                    % Analysis on Total Power                             %
                    %-----------------------------------------------------%
                    B1 = dt(blockNum~=i,:);    % training data
                    B2 = dt(blockNum==i,:);    % test data
                    C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
                    W = C1\B1;          % estimate weight matrix
                    C2 = (W'\B2')';     % estimate channel responses
                    
                    C2_total(f,iter,t,i,:,:) = C2;
                    
                    % shift eegs to common center
                    n2shift = ceil(size(C2,2)/2);
                    for ii=1:size(C2,1)
                        [~, shiftInd] = min(abs(posBins-tstl(ii)));
                        C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
                    end
                    
                    tf_total(f,iter,t,i,:) = mean(C2,1); % average shifted channel responses
                    C2_total_shifted(f,iter,t,i,:,:) = C2; %如果需要挑出高概率位置，则需要保存这个shifted但没有平均Bin的，Bin是倒数第二个维度
                    %-----------------------------------------------------%
                    
                end               
            end           
        end
        toc % 停止计时
    end
    
    %tf_tatol_mean = squeeze(mean(mean(tf_total,2),4));
    
    %% Save
    save_folder_path = fullfile(readdir, 'distractor_iem_allchans');    
    if ~exist(save_folder_path, 'dir'), mkdir(save_folder_path); end    
    fName = [save_folder_path filesep sublist{subno}(1:end-20) '_IEM.mat' ];
    %em.C2.evoked = C2_evoked;
    em.C2.total = C2_total;
    %em.tfs.evoked = tf_evoked;
    em.tfs.total = tf_total; 
    em.C2_total_shifted = C2_total_shifted;
    em.nBlocks = nBlocks;
    save(fName,'fdata_total','em','-v7.3');
end