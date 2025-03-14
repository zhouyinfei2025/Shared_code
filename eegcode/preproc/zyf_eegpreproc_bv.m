function zyf_eegpreproc_bv(cfg)
%% Change dir to files and make list of subject filenames
cd(cfg.readdir);        % 切换当前工作目录
sublist=dir('*.vhdr');  % 找到当前工作目录下所有.vhdr文件
sublist={sublist.name}; % 要注意的是sublist的顺序(subno)可能和被试顺序不一致，取决于文件命名方式

if cfg.single_run      
    [filename, filepath]=uigetfile('*.vhdr'); % 跳出文件对话框，选择要打开的文件
    % full_filepath = [filepath filename];
    sublist = {filename};
end
%% Loop around the subjects
for subno=1:length(sublist)
    %% Make the folder for output
    wdir = [cfg.writdir filesep sublist{subno}(1:4)];    
    if ~exist(wdir,'dir')
        mkdir(wdir)
    end
    % following variables
    outfilename1 = [ wdir filesep sublist{subno}(1:4) '_reref_hpfilt.mat' ];
    outfilename1b= [wdir filesep sublist{subno}(1:4) '_reref_hpfilt_4ica.mat' ];
    outfilename2 = [ wdir filesep sublist{subno}(1:4) '_rejectedtrials.mat' ];
    outfilename3 = [ wdir filesep sublist{subno}(1:4)  '_icaweights.mat' ];
    outfilename4 = [ wdir filesep sublist{subno}(1:4)  '_cleaned.mat' ];
    outfilename_exclude_heog = [ wdir filesep sublist{subno}(1:4) '_cleaned_exclude_heog.mat' ];
    outfilename_ica_reject_trial = [ wdir filesep sublist{subno}(1:4) '_ica_reject_trials.mat' ];
    % do all steps consecutively for one subject, or go to the step where you left off
    reload = true;
    go_on = true;
    
    while go_on
        
       %% Go to preprocessing step
        if ~exist(outfilename1,'file')  % go to preprocessing step 1
           %% Load data
            fprintf('Step1: Loading subject %s for re-referencing and high-pass filtering...\n',sublist{subno}(1:4));
            if cfg.single_run
                EEG = pop_loadbv(filepath,filename);
            else
                data_dir = [cfg.readdir '\'];
                EEG = pop_loadbv(data_dir,sublist{subno});
            end
            % erp: a structure array 类似EEG的一个结构体数组
            erp.event = struct( 'type', { EEG.event.type }, 'latency', {EEG.event.latency});
            event_number = length(EEG.event);
            erp.eventCodes = cell(1,event_number);
            for i=1:event_number
                erp.eventCodes{i} = EEG.event(i).type;
            end
            erp.eventTimes = round(cell2mat({erp.event.latency})); % Event Times
            erp.eventCodes = nan(1,length(erp.eventCodes));
            for ii = 1 : length(erp.eventTimes) %不要结尾
                s = erp.event(1,ii).type(2:4); % get characters 2-4 (i.e. the digits). 1 is actually stored as blank-blank-1.
                s = str2double(s);             % covert string to double BP专用
                erp.eventCodes(1,ii) = s;      % save the event code structure
            end
            for i=1:event_number
                EEG.event(i).type = erp.eventCodes(i);
            end
            
           %% Re-reference 
            if strcmp(cfg.ref,'average')  % 全脑平均
                EEG = pop_reref(EEG,[],'refstate','averef');
            else                          % 根据选择的两个位点进行平均
                EEG = pop_reref(EEG,[find(strcmpi(cfg.reref(1),{EEG.chanlocs.labels})) find(strcmpi(cfg.reref(2),{EEG.chanlocs.labels}))],'refstate','averef');
            end
            
           %% 建立双极导联 VEOG HEOG
            veogdat = squeeze(EEG.data(strcmpi(cfg.veog(1),{EEG.chanlocs.labels}),:)-EEG.data(strcmpi(cfg.veog(2),{EEG.chanlocs.labels}),:));
            heogdat = squeeze(EEG.data(strcmpi(cfg.heog(1),{EEG.chanlocs.labels}),:)-EEG.data(strcmpi(cfg.heog(2),{EEG.chanlocs.labels}),:));
            noincludeIndex = [cfg.veog cfg.heog] ;
            EEG = pop_select(EEG,'nochannel',noincludeIndex); % 删除眼电电极
            EEG.data(EEG.nbchan+1,:) = veogdat;
            EEG.data(EEG.nbchan+2,:) = heogdat;
            EEG.chanlocs(EEG.nbchan+1).labels = 'VEOG';
            EEG.chanlocs(EEG.nbchan+2).labels = 'HEOG';
            clear veogdat heogdat
            chan_number = EEG.nbchan+2;
            EEG.nbchan = chan_number;
            
           %% resample if asked for
            try
                if EEG.srate>resrate
                    EEG = pop_resample(EEG,resrate);
                end
            catch
                disp(['No resampling done. Sampling rate is ' num2str(EEG.srate)])
                
            end
           %% High-pass filter 
            try
                EEG  = pop_eegfiltnew(EEG, 'locutoff',cfg.highcut); % high-pass cut-off of final analysis
                EEGb = pop_eegfiltnew(EEG, 'locutoff',cfg.icacut);  % higher cut-off for ICA
                % EEG = pop_basicfilter(EEG,1:EEG.nbchan,'Filter','highpass','Design','butter','Cutoff',highcut,'Order',fltord,'RemoveDC','on'); 
                % EEGb = pop_basicfilter(EEG,1:EEG.nbchan,'Filter','highpass','Design','butter','Cutoff',icacut,'Order',fltord,'RemoveDC','on');% high edge cut-off for stable ICA 
            catch
                disp('eegfiltnew not present in eeglab package!')
            end
            
            %% 导入电极位置
            EEG  = pop_chanedit(EEG,'lookup',cfg.layout);
            EEGb  = pop_chanedit(EEGb,'lookup',cfg.layout);
            
           %% 比较marker数目是否正确
            EEGm = EEG;
            m=[EEGm.event.type];
            if ~isempty(cfg.trigger_subtract)
                m=m-(cfg.trigger_subtract);
            end
            
            for i=1:length(m)
                EEGm.event(i).type=m(i);
            end
            
            mtable=tabulate(m);
            mtable(mtable(:,2)==0,:)=[];
            mtable=mtable(:,[1 2]);
            disp(mtable)
            disp('Make sure the triggers and their number of ocurrence in the experiment make sense! If so, type in dbcont and hit enter...')
            keyboard
            
            % 确认无误后，把正确的event放回到EEG.event和EEGb.event中
            EEG.event  = EEGm.event; 
            EEGb.event = EEGm.event;
            
          %% Save data
           disp('Saving data of step I')
            if exist('NUM','var')
                save(outfilename1,'EEG','NUM','TXT','raw')
            else
                save(outfilename1,'EEG') % 把'EEG'这个结构体数组变量保存到outfilename1文件中
            end
            save(outfilename1b,'EEGb')  % 把'EEGb'这个结构体数组变量保存到outfilename1b文件中
            reload = false;
           
        elseif ~exist(outfilename2,'file')   % go to preprocessing step 2 
           %% Load data
            if reload
                fprintf('Step2: Loading subject %s for epoching and bad trial marking...\n',sublist{subno}(1:4));
                load(outfilename1)
                load(outfilename1b)
            else
                fprintf('Epoching and bad trial marking for subject %s...\n',sublist{subno}(1:4));
            end
           %% Epoch the data
            EEGb = pop_epoch( EEGb, [cfg.triggers{:,2}],[cfg.artiftime]);
            
            try
                EEGb = pop_rmbase( EEGb, []);
            catch
                EEGb = pop_rmbase( EEGb, [],[]);
            end
           %% Check for bad channels
            nchan = EEGb.nbchan;
            if cfg.inspect_chans
                figure(1);
                nblocks = 10; % divide data into data blocks of 10
                ntrials = floor(EEGb.trials/nblocks);
                for b=1:nblocks
                    subplot(3,ceil(nblocks/3),b)
                    newdata = reshape(squeeze(EEGb.data(1:nchan,:,1+((b-1)*ntrials):b*ntrials)),nchan,size(EEGb.data,2)*ntrials);
                    zstd = std(zscore(newdata),[],2);
                    topoplot(zstd,EEGb.chanlocs(1:nchan),'electrodes','on','maplimits',[min(zstd) max(zstd)]);
                    colorbar
                end
                colormap hot
                saveas(gcf,[cfg.writdir filesep sublist{subno}(1:4) filesep sublist{subno}(1:4) '_bad_electrode_check.png'])
                figure(2); topoplot([],EEGb.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEGb.chaninfo);
                pop_eegplot( EEGb, 1, 1, 1) % Channel data(scroll),to check flat or noisy channels.
                disp('挑选出坏通道,将可能的坏通道写在chans2interp这个文本里')
                disp('Make sure the you know which channel is a possibly broken one! If so, type in dbcont and hit enter...')
                keyboard
                close 'Scroll channel activities -- eegplot()'
            end
           %% Detection of muscle artifact based on Fieldtrip routines
            ft_defaults;
            cd(cfg.codedir); % 切换到代码文件夹，调用eeglab2ft函数
            FT_EEG = eeglab2ft(EEGb); % eeglab2ft这个函数在师姐给的代码文件包里
            
            cfg1 = [];
            cfg1.channel = {EEGb.chanlocs(1:nchan).labels}';
            dat = ft_selectdata(cfg1,FT_EEG);
            dat = dat.trial;
            dat_filt = cell(size(dat));
            
            cfg1=[];
            % these are standard fieldtrip settings taken from online tutorial
            cfg1.bpfreq = [110 140];
            cfg1.bpfiltord = 6;
            cfg1.bpfilttype = 'but';
            cfg1.hilbert = 'yes';
            cfg1.boxcar = 0.2;
            cfg1.channels = 1:nchan;
            cfg1.cutoff = cfg.artcutoff;
            cfg1.art_time = cfg.artiftime; % window within which artifacts are not allowed; note: also use this window for ICA!
            
            %-% loop over trials
            reverseStr='';
            numtrl=EEGb.trials;
            tidx = dsearchn(FT_EEG.time{1}',cfg1.art_time')';
    
            
            for ei=1:numtrl
                
                % display progress
                msg = sprintf('Filtering trial %i/%i...',  ei,numtrl);
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
                
                % filter in high-freq band
                tmpdat = ft_preproc_bandpassfilter(dat{ei},EEGb.srate, cfg1.bpfreq, cfg1.bpfiltord, cfg1.bpfilttype);
                tmpdat = ft_preproc_hilbert(tmpdat,cfg1.hilbert);
                nsmp = round(cfg1.boxcar*EEGb.srate);
                if ~rem(nsmp,2)
                    % the kernel should have an odd number of samples
                    nsmp = nsmp+1;
                end
                tmpdat = ft_preproc_smooth(tmpdat, nsmp); % better edge behaviour
                dat_filt{ei} = double(tmpdat(:,tidx(1):tidx(2)));
                
                if ei==1
                    sumval = zeros(size(dat_filt{1},1), 1);
                    sumsqr = zeros(size(dat_filt{1},1), 1);
                    numsmp = zeros(size(dat_filt{1},1), 1);
                    numsgn = size(dat_filt{1},1);
                end
                
                % accumulate the sum and the sum-of-squares
                sumval = sumval + sum(dat_filt{ei},2);
                sumsqr = sumsqr + sum(dat_filt{ei}.^2,2);
                numsmp = numsmp + size(dat_filt{ei},2);
            end
            fprintf('\n')
            
            % avg and std
            datavg = sumval./numsmp;
            datstd = sqrt(sumsqr./numsmp - (sumval./numsmp).^2);
            
            zmax = cell(1, numtrl);
            zsum = cell(1, numtrl);
            zindx = cell(1, numtrl);
            
            indvec = ones(1,numtrl);
            for ei = 1:numtrl
                % initialize some matrices
                zmax{ei}  = -inf + zeros(1,size(dat_filt{ei},2));
                zsum{ei}  = zeros(1,size(dat_filt{ei},2));
                zindx{ei} = zeros(1,size(dat_filt{ei},2));
                
                nsmp          = size(dat_filt{ei},2);
                zdata         = (dat_filt{ei} - datavg(:,indvec(ei)*ones(1,nsmp)))./datstd(:,indvec(ei)*ones(1,nsmp));  % convert the filtered data to z-values
                zsum{ei}   = nansum(zdata,1);                   % accumulate the z-values over channels
                [zmax{ei},ind] = max(zdata,[],1);            % find the maximum z-value and remember it
                zindx{ei}      = cfg1.channels(ind);                % also remember the channel number that has the largest z-value
                
                zsum{ei} = zsum{ei} ./ sqrt(numsgn);
            end
            cfg1.threshold = median([zsum{:}]) + abs(min([zsum{:}])-median([zsum{:}])) + cfg1.cutoff;
            
            figure(3)
            plot([zsum{:}])
            hold on
            plot([1 length([zsum{:}])], [cfg1.threshold cfg1.threshold],'m')
            saveas(gcf,[wdir filesep sublist{subno}(1:4) '_zvals_rejectedtrials.png'])
            rejected_trials = zeros(1,numtrl);
            for ei=1:numtrl
                if sum(zsum{ei}>cfg1.threshold)>0
                    rejected_trials(ei) = 1;
                end
            end
            
            rejeegplottmp = trial2eegplot(  rejected_trials, zeros(EEGb.nbchan,EEGb.trials), EEGb.pnts, [1 1 0.5]);
            eegplot(EEGb.data,'srate', EEGb.srate,'limits', [EEGb.xmin EEGb.xmax]*1000 ,'eloc_file', EEGb.chanlocs,'winrej',rejeegplottmp);
            disp('检查trial，点击认为不好的trials')
            fprintf(['At this point, do following checks:\n'...
                '- check if rejected trials make sense, also check zsum and threshold figure\n' ...
                '- turn back to topoplots of possible bad channels, and verify (add them to .txt file with name as specfied in cfg.chanfilename)\n'...
                '- bad trials and bad channels need to be removed before accurate ICA can be done!\n'...
                '>> DO NOT close the eegplot figure if you want to change the marked trials!\n'...
                '>> Click on epochs to add/remove marked trials; newly added trials are blue, already marked trials are yellow\n'...
                '>> If satisfied, just type "dbcont" in the command window. The figure will close automatically.\n\n'])
            keyboard
            
            close(figure(1));
            close(figure(2));
            close(figure(3));
            h=get(gcf);
            tmprej = h.UserData.winrej;
            tmprej = eegplot2trial(tmprej,EEGb.pnts,EEGb.trials);
            close 'Scroll activity -- eegplot()'

            if sum(tmprej)~=sum(rejected_trials)
                fprintf('Manually Added/removed %i trials to the automatically detected trials...\n', sum(tmprej)-sum(rejected_trials))
                rejected_trials = tmprej;
            end
  
            

           %% Save data
            disp('Saving data of step II')
            if exist('wrong','var')
                save(outfilename2,'rejected_trials','cfg1','wrong');
            else
                save(outfilename2,'rejected_trials','cfg1'); % the cfg variable contains the rejection settings, so these can be traced back (e.g. z-val cutoff)
            end
            reload = false;
        elseif ~exist(outfilename3,'file')  % go to preprocessing step 3
           %% Load data
            if reload
                fprintf('Step3: Loading subject %s for ICA...\n',sublist{subno}(1:4));
                load(outfilename1b)
                load(outfilename2)
                EEGb = pop_epoch( EEGb, [cfg.triggers{:,2}],[cfg.artiftime]);
                try
                    EEGb = pop_rmbase( EEGb, []);
                catch
                    EEGb = pop_rmbase( EEGb, [],[]);
                end
                    
            else
                fprintf('Running ICA for subject %s...\n',sublist{subno}(1:4));
            end
            
           %% Remove trials with artifacts detected in previous step
            EEGb = pop_select(EEGb,'notrial',find(rejected_trials));
            
           %% Remove bad channels detected in previous step
            fid=fopen([cfg.parent_folder filesep cfg.chanfilename],'r');
            chans2interp={};
            while ~feof(fid)
                aline = regexp(fgetl(fid),' ','split');
                chans2interp{size(chans2interp,1)+1,1}=aline{1};
                for i=2:length(aline)
                    chans2interp{size(chans2interp,1),i}=aline{i};
                end
            end
            nchan = EEGb.nbchan;
            %如果要排除双极导联，此处改为nchan = EEG.nbchan-2,因为双极导联是最后两个，所以只要-2即可
            %跑ICA不排除双极导联
            chanind=1:nchan;
            subject_prefix = sublist{subno}(1:4);
            chans = chans2interp(strcmpi(chans2interp(:,1),subject_prefix),2:end);
            subject_name = chans2interp(:,1);
            subject_id = find(strcmpi(subject_name,subject_prefix));
            if ~isempty(chans{1})
                bad_chansidx = zeros(1,length(chans));
                for c=1:length(chans)
                    if ~isempty(chans{c}), bad_chansidx(c) = find(strcmpi({EEGb.chanlocs.labels},chans(c))); end
                end
                bad_chansidx(bad_chansidx==0)=[];
                chanind(bad_chansidx)=[];
            end
            
           %% Run ICA
            EEGb=pop_runica(EEGb,'icatype','runica','chanind',chanind,'extended',1,'interrupt','on');
            
            ICAEEG.icawinv = EEGb.icawinv;
            ICAEEG.icasphere = EEGb.icasphere;
            ICAEEG.icaweights = EEGb.icaweights;
            ICAEEG.icachansind = EEGb.icachansind;
            
           %% Save data
            disp('Saving data of step III')
            save(outfilename3,'ICAEEG')
            reload = true;
            
           %% 将EEGb保存为.set文件，可在图形交互界面查看ICA结果
            savename = [sublist{subno}(1:4) '.set'];
            pop_saveset( EEGb, 'filename',savename,'filepath', wdir);
            %eeglab; % 打开图形交互界面
            
            
        elseif ~exist(outfilename4,'file')
            %% Load data 
            if reload
                fprintf('Step4: Loading subject %s for final cleaning...\n',sublist{subno}(1:4));
                load(outfilename1)
                load(outfilename1b)
                load(outfilename2)
                load(outfilename3)
                
                EEGb = pop_epoch( EEGb, [cfg.triggers{:,2}],[cfg.artiftime]);
                EEGb = pop_rmbase( EEGb, [],[]);
                EEGb = pop_select(EEGb,'notrial',find(rejected_trials));
                EEGb.icachansind = ICAEEG.icachansind;
                EEGb.icasphere = ICAEEG.icasphere;
                EEGb.icaweights = ICAEEG.icaweights;
                EEGb.icawinv = ICAEEG.icawinv;
                EEGb = eeg_checkset(EEGb);

            else
                fprintf('Final cleaning steps for subject %s ...\n',sublist{subno}(1:4));
            end
            
           %% Inspect components to remove
            if cfg.inspect_ica
                pop_selectcomps(EEGb,[1:20]); % plot first 20 ICs
                saveas(gcf,[cfg.writdir filesep sublist{subno}(1:4) filesep sublist{subno}(1:4) '_ICs_check.png'])
                pop_eegplot( EEGb, 0, 1, 1);  % plot the component activation,成分卷轴图
                
                inspecting = true;
                while inspecting
                    disp('检查ica成分，判断哪些是眼动，将需要剔除的成分写在ICs2remove这个文本里')
                    prompt = 'Which component do you want to inspect? (type number or "done" when done) --> ';
                    comp2check = input(prompt,'s');
                    if strcmp(comp2check,'done')
                        inspecting = false;
                    else
                        pop_prop( EEGb, 0, str2double(comp2check));
                    end
                end
            end
            close all; % close all the open figure windows
          
           %% Mark bad channels 
            fid=fopen([cfg.parent_folder filesep cfg.chanfilename],'r');
            chans2interp={};
            while ~feof(fid)
                aline=regexp(fgetl(fid),' ','split');
                chans2interp{size(chans2interp,1)+1,1}=aline{1};
                for i=2:length(aline)
                    chans2interp{size(chans2interp,1),i}=aline{i};
                end
            end
            
            nchan = EEGb.nbchan;
            chanind=1:nchan;
            subject_prefix = sublist{subno}(1:4);
            chans = chans2interp(strcmpi(chans2interp(:,1),subject_prefix),2:end);
            subject_name = chans2interp(:,1);
            subject_id = find(strcmpi(subject_name,subject_prefix));
            bad_chansidx = [];
            if ~isempty(chans{1})
                bad_chansidx = zeros(1,length(chans));
                for c=1:length(chans)
                    if ~isempty(chans{c}), bad_chansidx(c) = find(strcmpi({EEG.chanlocs.labels},chans(c))); end
                end
                bad_chansidx(bad_chansidx==0)=[];
                chanind(bad_chansidx)=[];
           end
           %% Mark ICs to remove
            fid=fopen([cfg.parent_folder filesep cfg.icafilename],'r');
            comps2remove={};
            while ~feof(fid)
                aline=regexp(fgetl(fid),'\t','split')
                comps2remove{size(comps2remove,1)+1,1}=aline{1}
                comps2remove{size(comps2remove,1)  ,2}=sscanf(aline{2},'%g');
            end
            
           %% Do the final cleaning:
            %- epoch in wider window
            %- baseline correct
            %- reject trials
            %- add ICA weights
            %- remove EOG
            %- interpolate channels
            %- remove components
            
            EEG = pop_epoch( EEG, [cfg.triggers{:,2}],cfg.epochtime);
         
            try
                EEG = pop_rmbase( EEG, []);
            catch
                EEG = pop_rmbase( EEG, [],[]);
            end
            for ei=1:EEG.trials
                EEG.epoch(ei).trialnum = ei;
            end
 
            if exist('NUM','var')
                for ei=1:EEG.trials
                    EEG.epoch(ei).correct = NUM(ei,1);
                    EEG.epoch(ei).RT = NUM(ei,2);
                end
            end
      
        
            EEG = pop_select(EEG,'notrial',find(rejected_trials));
     
            
            EEG.icachansind = ICAEEG.icachansind;
            EEG.icasphere = ICAEEG.icasphere;
            EEG.icaweights = ICAEEG.icaweights;
            EEG.icawinv = ICAEEG.icawinv;
            EEG = eeg_checkset(EEG);
            
            clear ICAEEG
            nchan = length(EEG.chanlocs);
            subject_name = comps2remove(:,1);
            subject_id = find(strcmpi(subject_name,sublist{subno}(1:4)));
       
            if sum(bad_chansidx)>0
                
                EEG2 = pop_select(EEG,'nochannel',bad_chansidx);
                EEG2 = pop_subcomp( EEG2, comps2remove{subject_id,2}, 0);
                %EEG2 = pop_subcomp(EEG2); % 图形交互界面，手动填入要删除的眼电成分
                %for d = 2:length(ICs2remove)
                   %EEG2 = pop_subcomp( EEG2, str2num(ICs2remove{d}));
                %end
                
                
               %% put IC-cleaned channels back in original EEG structure and interpolate bad ones
                good_chans = 1:EEG.nbchan;
                good_chans(bad_chansidx)=[];
                EEG.data(good_chans,:,:) = EEG2.data;
                EEG = eeg_interp(EEG,bad_chansidx);
                
                clear EEG2
            else
                EEG = pop_subcomp( EEG, comps2remove{subject_id,2}', 0);
                %EEG = pop_subcomp(EEG); % 图形交互界面，手动填入要删除的眼电成分
                %for d = 2:length(ICs2remove)
                   %EEG = pop_subcomp( EEG, str2num(ICs2remove{d}));
                %end
            end
            EEG.reject = EEGb.reject;
            clear EEGb
            
            % 根据水平眼动电极来检查数据
            if cfg.reject_heog
                pop_eegplot(EEG,1,0,0);
                inspecting = true;
                while inspecting
                    disp('根据水平眼动电极来检查数据，如果需要剔除trial，点击认为不好的trial')
                    prompt = 'according to heog which trials want to exclude? (type number or "done" when done) --> ';
                    comp2check = input(prompt,'s');
                    if strcmp(comp2check,'done')
                        inspecting = false;
                    end
                end
                clear tmprej
                h=get(gcf);
                tmprej = h.UserData.winrej;
                if length(tmprej) ~= 0;
                    tmprej = eegplot2trial(tmprej,EEG.pnts,EEG.trials); %% to get the rejecet trial for ica plot
                    heog_reject_trials = find(tmprej);
              
                    save(outfilename_exclude_heog,'heog_reject_trials')
                    EEG = pop_select(EEG,'notrial',heog_reject_trials);
                  
                   
                end
                
            end
           %% Save data
           
            fprintf('Saving data of step IV, done preprocessing for subject %s!',sublist{subno}(1:4));
            save(outfilename4,'EEG','comps2remove','bad_chansidx');
            
            %Savename = [sublist{subno}(1:4) '_cleaned.set'];
            %pop_saveset( EEG, 'filename',Savename,'filepath', wdir);
            
            go_on = false; % this gets you out of the while loop, so the subject loop continues
            
        else
            go_on = false; % this gets you out of the while loop, so the subject loop continues
        end 
    end
end
end

