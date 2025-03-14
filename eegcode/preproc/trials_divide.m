%% 将学习条件预处理之后的数据分为前5个block和后五个block
% EEG_pre 包含前5个block的数据
% EEG_post 包含后5个block的数据

%% It's always good to start with a clean sheet
clear,clc

%% Get all the data file names
[~, filepath]=uigetfile('*.mat'); 
sublist = dir(fullfile(filepath,'*.mat'));
sublist={sublist.name};
sublist=sublist(2:2:56); % learning condition

%%
for subno = 1:28
    %% Load data
    clear EEG EEG_pre EEG_post
    dname = sublist{subno};
    fprintf('Loading subject %s for analysis ...\n',dname);
    load([filepath filesep dname])

    %% 找到前5个block和后5个block的分界试次，第600个试次
    for i = 1:EEG.trials
        if EEG.epoch(i).trialnum == 598 % 找不到600，就找599，598...以此类推。
            disp(i)
            break
        end
    end

    %%
    EEG_pre = EEG;
    EEG_pre.trials = i;
    EEG_pre.data = EEG.data(:,:,1:i);
    EEG_pre.epoch = EEG.epoch(1:i);

    EEG_post = EEG;
    EEG_post.trials = EEG.trials - i;
    EEG_post.data = EEG.data(:,:,i+1:end);
    EEG_post.epoch = EEG.epoch(i+1:end);
    
    %% Save the data    
    save([filepath dname(1:4) '_' 'former_5blocks_cleaned.mat'],'EEG_pre','-v7.3');
    save([filepath dname(1:4) '_' 'latter_5blocks_cleaned.mat'],'EEG_post','-v7.3');
    
end
