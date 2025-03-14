%% It's always good to start with a clean sheet
clear, close all, warning('off','all'),clc

%% part I: loading raw data and filter
cfg = []; % structure array
cfg.eeglab_path = 'D:\package\eeglab2020_0';          % path to eeglab package
cfg.ft_path         = 'D:\package\fieldtrip-20211020';   % path to fieldtrip package
cfg.layout          = 'D:\package\new-64.ced';               % path to channel layout
cfg.parent_folder='E:\01';                                             % Folder for Project
cfg.readdir        = 'E:\01\eegdata\raw_data';                % path to raw files
cfg.behavdir     = 'E:\01\behav_data';                          % path to behavioral data
cfg.writdir         = 'E:\01\eegdata\proc_data';              % path to output files
cfg.codedir       = 'E:\01\eegcode\preproc';                 % path to codes

% Define a project name  
%cfg.projectname = '2001'; 

cfg.nchan       = 61;                  % the number of channels/electrodes
cfg.ref         = 'mastoids';          % reference the data to mastoids
cfg.reref       = {'TP9','TP10'};      % channel labels for re-reference
cfg.veog        = {'T8','FT10'};       % channel labels that correspond to VEGO
cfg.heog        = {'FT9','T7'};        % channel labels that correspond to HEGO
cfg.resrate     = 512;                 % new sampling rate 
cfg.highcut     = 0.1;                 % high-pass cut-off of final analysis
cfg.icacut      = 1.5;                 % higher cut-off for ICA 
cfg.fltord      = 2;                   % filter order
cfg.read_behav  = false;        % 是否要将行为数据正确率和RT添加到EEG数据中，需要填true,不需要填false,行为数据格式见behav这个文件夹下       
cfg.single_run  = false;           % 是否需要一次只跑一个被试，如果需要，就填true,会出现窗口要求你选择被试数据
cfg.reject_heog = false;          % 剔除眨眼后，是否还需要手动来剔除可能的水平眼动trial,是的选true

%% part II: epoching and trial rejection marking
triggers = {'NO_DIS_1' {'120'};
            'NO_DIS_2' {'130'};
            'NO_DIS_3' {'140'};
            'NO_DIS_4' {'150'};
            
            'HIGH_DIS_1' {'131'};
            'HIGH_DIS_1' {'141'};
            'HIGH_DIS_1' {'151'};
            
            'LOW_DIS_1' {'122'}
            'LOW_DIS_2' {'123'}
            'LOW_DIS_3' {'124'}
            'LOW_DIS_4' {'133'}
            'LOW_DIS_5' {'134'}
            'LOW_DIS_6' {'142'}
            'LOW_DIS_7' {'144'}
            'LOW_DIS_8' {'152'}
            'LOW_DIS_9' {'153'};};

% [-0.5 3.25] time window of interest, including the baseline period
cfg.epochtime   = [-2 4.75];  % a rather wide epoch window (1.5 second longer at each side) for time frequency analysis
cfg.artiftime   = [-0.3 3.55];  % a rather wide window (300ms longer at each side) for rejecting the artifacts only
cfg.artcutoff   = 12;           % Standard for automatically cutting off the artifacts with fieldtrip
cfg.triggers    = triggers;
clear triggers
cfg.trigger_subtract = [];    % depending on physical lab settings, sometimes weird high numbers get added to the trigger values specified in your experiment script
cfg.inspect_chans = true;     % pauses function and shows figures with topomaps of variance, to highlight bad channels

%% part III: ICA
cfg.chanfilename = 'chans2interp.txt';

%% part IV: final cleaning
cfg.icafilename = 'ICs2remove.txt';
cfg.inspect_ica = true;

%% now run it
zyf_eegpreproc_bv(cfg);