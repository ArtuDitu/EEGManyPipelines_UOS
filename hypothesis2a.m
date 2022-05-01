%Hypothesis 2
%There are effects of image novelty within the time-range from 300–500 ms ...
%a.on EEG voltage at fronto-central channels.


cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS % path to data
addpath(genpath('/net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines')) 
eeglabpath = fileparts(which('eeglab.m'));

eeglab; % start toolbox


%add triggers
% 2nd digit 0 and 1
trigger_new = {1010 1011 1019 1020 1021 1029 1030 1031 1039 1040 1041 1049 1090 1091 1099 2010 2011 2019 2020 2021 2029 2030 2031 2039 2040 2041 2049 2090 2091 2099};
trigger_old = {1110 1111 1119 1120 1121 1129 1190 1191 1199 2110 2111 2119 2120 2121 2129 2190 2191 2199};
window_epoch = [-.2 .8]; 

%load data

cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS/derivatives % path to data, 
addpath(genpath('/net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS/derivatives')) % add your folder to the path
list_of_files = dir('**/final*'); %save all cleaned data sets in one file

%% STEP 1

for eeg_file = 1:size(list_of_files)
    %load file
    cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS/derivatives
    EEG = pop_loadset(list_of_files(eeg_file).name);
    %epoch data
    EEG_new = pop_epoch( EEG, trigger_new, window_epoch, 'epochinfo', 'yes');
    EEG_old = pop_epoch( EEG, trigger_old, window_epoch, 'epochinfo', 'yes');
    %eegplot(EEG_old.data,'srate',EEG_old.srate,'eloc_file',EEG_old.chanlocs,'events',EEG_old.event)
    
    %remove baseline
    EEG_new = pop_rmbase(EEG_new, [-199 0]);
    EEG_old = pop_rmbase(EEG_old, [-199 0]);
    
    %average data on epoch level
    EEG_new.data=mean(EEG_new.data(:,:,:),3);
    EEG_old.data=mean(EEG_old.data(:,:,:),3);
    
    %concatenate data for all participants
    if eeg_file==1
        EEG_new_all=EEG_new;
        EEG_old_all= EEG_old;
    else
        EEG_new_all.data=cat(3, EEG_new_all.data, EEG_new.data);
        EEG_old_all.data=cat(3, EEG_old_all.data, EEG_old.data);
    end
end

%eegplot(EEG_new_all.data,'srate',EEG_new_all.srate,'eloc_file',EEG_new_all.chanlocs)
%eegplot(EEG_old_all.data,'srate',EEG_old_all.srate,'eloc_file',EEG_old_all.chanlocs)


%% STEP 2
% Keep and average over participants

channels_to_keep = {'F1', 'FC1', 'Fz', 'FCz', 'F2', 'FC2'} % 6 channels (Frontal) / 
EEG_new_all = pop_select(EEG_new_all, 'channel', channels_to_keep); % keep channel 
EEG_old_all = pop_select(EEG_old_all, 'channel', channels_to_keep); % keep channel 

% time window of 300 to 500 ms
EEG_new_all.data = EEG_new_all.data(:,257:359,:)
EEG_old_all.data = EEG_old_all.data(:,257:359,:)

%eegplot(EEG_new_all.data,'srate',EEG_new_all.srate,'eloc_file',EEG_new_all.chanlocs)
%eegplot(EEG_old_all.data,'srate',EEG_old_all.srate,'eloc_file',EEG_old_all.chanlocs)


average_over_participant_new_all= squeeze(mean(mean(EEG_new_all.data, 1), 2));
average_over_participant_old_all= squeeze(mean(mean(EEG_old_all.data, 1), 2));

[h,p, ci, stats] = ttest(average_over_participant_new_all,average_over_participant_old_all)
% h = 1 // p=0.0034 // ci= -0.4127, -0.0894 // stats =   tstat: -3.1639, df: 32, sd: 0.4559
% t(32) = -3.1639, p=0.0034

%%


