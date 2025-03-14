%% Time-frequency analysis


%% It's always good to start with a clean sheet
clear, close all, warning('off','all'),clc

%% Add toolbox
%eeglab_dir = 'E:\package\eeglab2020_0';
%addpath(eeglab_dir)
%eeglab
%% triggers and events
triggers = {'L_D_V_T_1' {'134'};
            'L_D_V_T_2' {'144'};
            'H_D_V_T_1' {'131'};
            'H_D_V_T_2' {'141'};
            
            'N_D_L_T' {'150'};
            'N_D_H_T' {'120'};
            
            'H_D_L_T' {'151'};
            
            'V_D_V_T_1' {'133'};
            'V_D_V_T_2' {'142'};
            'V_D_L_T_1' {'152'};
            'V_D_L_T_2' {'153'};
            'N_D_V_T_1' {'130'};
            'N_D_V_T_2' {'140'}
            'L_D_H_T_1' {'122'}
            'L_D_H_T_2' {'123'}
            'L_D_H_T_3' {'124'}};

% Not differentiate condition
noconditions = [triggers{:,2}];
Cond = {'noconditions'};
All_cond = {'Noconditions'};

%% Get all the data file names
[~, filepath]=uigetfile('*.mat'); 
sublist = dir(fullfile(filepath,'*.mat'));
sublist={sublist.name};


%% Analysis for baseline condition
for subno = 1:length(sublist)
    %% Load data
    dname = sublist{subno};
    fprintf('Loading subject %s for analysis ...\n',dname(1:4));
    load([filepath filesep dname])
    EEG = EEG_pre;
    
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
    
    %% Read all the events and load data 
    eegdat = cell(1,length(Cond));
    eegdat{1} = EEG.data(:,:,:);
    
    %% Parameters for tfdecomp
    cfg = [];

    % -- Path/filenames for saving:
    cfg.writdir = [filepath 'TF_E2']; 
    cfg.filename = [dname(1:4) '_tfdecomp_output.mat'];

    cfg.srate = EEG.srate; 
    cfg.eegtime = EEG.times; 
    cfg.channels = 1:59;
    cfg.chanlocs = EEG.chanlocs;
    cfg.frequencies = [2.4,4.29,5.45,7.5,8,9,10,11,12,20,25,30];
    cfg.cycles = [3 12]; % min max number of cycles used for min max frequency
    cfg.times2save = -500:2:3250; % epoch of interest
    cfg.basetime = [-500 -300]; % pre-stim baseline
    cfg.baselinetype = 'conavg'; % 'conavg' or 'conspec'
    cfg.erpsubract = false; % if true, non-phase-locked (i.e. "induced") power will be computed by subtracting the ERP from each single trial
    cfg.matchtrialn = false; % if true, conditions will be equated in terms of trial count, so SNR(signal noise ratio) is comparable across conditions
    
    
    % -- other administrative stuff:
    cfg.report_progress = true;
    cfg.save_output = true;
    cfg.overwrite = false;
    
    %% Call the tfdecomp function
    [tf_pow,tf_phase,dim] = tfdecomp(cfg,eegdat);
        
end
