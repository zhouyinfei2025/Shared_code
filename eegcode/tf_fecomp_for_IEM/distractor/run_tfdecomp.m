%% run_tfdecomp
% 根据干扰出现位置分类，平均8-12Hz，保留试次
% 2025/2/23


%% It's always good to start with a clean sheet
clear, close all, warning('off','all'),clc

%% triggers and events

fix_loc = [131,141,151];   % 干扰出现在高概率干扰位置
right    = [122,142,152];   %  干扰出现在右侧
left       = [123,133,153];   % 干扰出现在左侧 
opp_loc  = [124,134,144];   % 干扰出现在高概率位置对侧

Cond = {'fix_loc';'right';'opp_loc';'left'};
All_cond = {'Fix_loc';'Right';'Opp_loc';'Left'};

%% Get all the data file names
%[~, filepath]=uigetfile('*.mat'); 
filepath = 'E:\01\eegdata\proc_data';
sublist = dir(fullfile(filepath,'*.mat'));
sublist={sublist.name};

%% Analysis
for subno = 1:length(sublist)
    %% Load data
    dname = sublist{subno};
    fprintf('Loading subject %s for analysis ...\n',dname(1:4));
    load([filepath filesep dname])
    
    %% Delete error trials
    fprintf('Loading subject %s for error trial delete...\n',dname(1:4));
    %EEG = pop_selectevent(EEG,'latency','-100.0<=2000.01','deleteevents','on');
    numtrl=EEG.trials;
    rejected_trials = zeros(1,numtrl);
    for m=1:numtrl   
        for n = 1:numel(EEG.epoch(m).eventtype)
           if cell2mat(EEG.epoch(m).eventtype(n)) == 111 || cell2mat(EEG.epoch(m).eventtype(n)) == 222
               rejected_trials(m) = 1;
           end
        end
    end
    EEG=pop_select(EEG,'notrial',find(rejected_trials));
    
    %% Read all the events and load data for different condition
    clear events;
    events = zeros(1,EEG.trials);
    for cEv = 1:numel(EEG.epoch)
         events(cEv) = cell2mat(EEG.epoch(cEv).eventtype(cell2mat(EEG.epoch(cEv).eventlatency)==0));
    end
    eegdat = cell(1,length(Cond));
    for j =1:length(Cond)
        % Get data
        trialtype = ismember(events,eval(cell2mat(Cond(j))));
        eegdat{j} = EEG.data(:,:,trialtype==1);
    end
    
    %% Parameters for tfdecomp
    cfg = [];

    % -- Path/filenames for saving:
    cfg.writdir = [filepath filesep 'TF_cond_new']; 
    cfg.filename = [dname(1:4) '_tfdecomp_output.mat'];

    cfg.srate = EEG.srate; 
    cfg.eegtime = EEG.times; 
    cfg.channels = 1:59;
    cfg.chanlocs = EEG.chanlocs;
    cfg.frequencies = [2.4,4.29,5.45,7.5,8,9,10,11,12,20,25,30];
    cfg.cycles = [3 12]; % min max number of cycles used for min max frequency
    cfg.times2save = -500:2:3250; % epoch of interest,每隔2ms采一个点，保持500Hz的采样率
%    cfg.basetime = [-500 -300]; % pre-stim baseline
%    cfg.baselinetype = 'conavg'; % 'conavg' or 'conspec'
    cfg.erpsubract = false; % if true, non-phase-locked (i.e. "induced") power will be computed by subtracting the ERP from each single trial
    cfg.matchtrialn = false; % if true, conditions will be equated in terms of trial count, so SNR(signal noise ratio) is comparable across conditions
    cfg.singletrial = true;
    
    % -- other administrative stuff:
    cfg.report_progress = true;
    cfg.save_output = true;
    cfg.overwrite = false;
    
    %% Call the tfdecomp function
    [tf_pow, z_pow, dim] = tfdecomp(cfg,eegdat);
        
end
    
