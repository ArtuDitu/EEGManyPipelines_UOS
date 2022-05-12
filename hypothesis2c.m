%Matlab version: R2020a
% eeglab version: 2020_0

%Hypothesis 2c
%There are effects of image novelty (i.e., between images shown for the first time/new
%vs. repeated/old images) within the time-range from 300–500 ms on alpha power at posterior channels..
% Posterior channels used: CPz, Pz, POz (32, 31, 30)


cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS % path where I have data, you should change it to the path on your PC
addpath(genpath('/net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines')) % add your folder to the path
eeglabpath = fileparts(which('eeglab.m'));
% load eeglab (it is a toolbox that we use to work with eeg data, you can
% download it online, and add it in the same folder as above)
eeglab; % start toolbox

%add TFCE toolbox 
addpath(genpath('/net/store/nbp/projects/joint_error/EEG_Belt/TFCE/ept_TFCE-matlab-master'));
%add triggers
trigger_new = {1010 1011 1019 1020 1021 1029 1030 1031 1039 1040 1041 1049 1090 1091 1099 2010 2011 2019 2020 2021 2029 2030 2031 2039 2040 2041 2049 2090 2091 2099};
trigger_old = {1110 1111 1119 1120 1121 1129 1190 1191 1199 2110 2111 2119 2120 2121 2129 2190 2191 2199};
window_epoch = [-.2 .8]; 
posterior_channels=[30, 31, 32]; %channesl i want to investigate
pnts = 103; % number of timepoints btw 300 and 500ms

power_new_all=[];
power_old_all=[];
%load data

cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS/derivatives % path where I have data, you should change it to the path on your PC
addpath(genpath('/net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS/derivatives')) % add your folder to the path
list_of_files = dir('**/final*'); %save all cleaned data sets in one file


%% STEP 1
% create two datasets, one which contains only new images and one
% which contains only repeated scenery, epoch [-200,800]
for eeg_file = 1:size(list_of_files)
    %load file
    %cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS/derivatives
    EEG = pop_loadset(list_of_files(eeg_file).name);
    %epoch data
    EEG_new = pop_epoch( EEG, trigger_new, window_epoch, 'epochinfo', 'yes');
    EEG_old = pop_epoch( EEG, trigger_old, window_epoch, 'epochinfo', 'yes');
    %eegplot(EEG_manmade.data,'srate',EEG_manmade.srate,'eloc_file',EEG_manmade.chanlocs,'events',EEG_manmade.event)
    
    %remove baseline
    EEG_new = pop_rmbase(EEG_new, [-199 0]);
    EEG_old = pop_rmbase(EEG_old, [-199 0]);
    
    %remove veog and heog
    EEG_new=pop_select(EEG_new, 'nochannel', {'VEOG', 'HEOG'});
    EEG_old=pop_select(EEG_old, 'nochannel', {'VEOG', 'HEOG'});
    
    %find out indices of posterior channels
    chanlist=[];
    for i=1:70
        chanlist=[chanlist, {EEG_new.chanlocs(i).labels}];
    end
    posterior_channels={'CPz', 'POz', 'Pz'};
    chan_Index=[];
    for j=1:3
        chan_Index = [chan_Index, find(strcmp(chanlist,posterior_channels(j)))];
    end
    
    %select posterior channels
    EEG_new.data = squeeze(mean(EEG_new.data(chan_Index,:,:),1));
    EEG_old.data = squeeze(mean(EEG_old.data(chan_Index,:,:),1));
    
    n=.201*EEG_new.srate;

    window = hann(n);
    
    %calculate channel power for both conditions
    chanpowr_new = (2*abs(fft(EEG_new.data(257:359,:).*window,[],1))/pnts).^2; % EEG.data here is our 300-500 ms vector ad EEG.pnts is amount of data points 
    chanpowr_new = mean(chanpowr_new,2); %average over trials 

    %chanpowr_new_all=cat(2,chanpowr_new_all,chanpowr_new); 
    % select alpha
    power_new =mean(chanpowr_new(8:12)); 
    power_new_all=cat(2,power_new_all,power_new); 
    
    chanpowr_old = (2*abs(fft(EEG_old.data(257:359,:).*window,[],1))/pnts).^2; % EEG.data here is our 300-500 ms vector ad EEG.pnts is amount of data points - careful, the parentheses in these order are crucial
    chanpowr_old = mean(chanpowr_old,2); %average over trials 

    %chanpowr_old_all=cat(2,chanpowr_old_all,chanpowr_old); 
    % select alpha
    power_old =mean(chanpowr_old(8:12));
    power_old_all=cat(2,power_old_all,power_old); 
end



%% STEP 2
% run non parametric
[p,h,stats] = signrank(power_new_all, power_old_all)

%run parametric
%[h,p,ci,stats]=ttest(power_new_all, power_old_all);
