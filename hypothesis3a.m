%Matlab version: R2020a
% eeglab version: 2020_0
% TFCE toolbox: https://github.com/Mensen/ept_TFCE-matlab

%Hypothesis 3a
%There are effects of successful recognition of old images (i.e., a difference between
%old images correctly recognized as old [hits] vs. old images incorrectly judged as new
%[misses]) 
%a.on EEG voltage at any channels, at any time.

cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS % path to data, 
addpath(genpath('/net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines')) % add your folder to the path
eeglabpath = fileparts(which('eeglab.m'));

eeglab; % start toolbox

%add TFCE toolbox 
addpath(genpath('/net/store/nbp/projects/joint_error/EEG_Belt/TFCE/ept_TFCE-matlab-master'));
%add triggers
% 3rd digit 1 and 2
trigger_hit = {1010 1011 1019 1110 1111 1119 2010 2011 2019  2110 2111 2119};
trigger_miss = {1020 1021 1029 1120 1121 1129 2020 2021 2029 2120 2121 2129};
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
    EEG_hit = pop_epoch( EEG, trigger_hit, window_epoch, 'epochinfo', 'yes');
    EEG_miss = pop_epoch( EEG, trigger_miss, window_epoch, 'epochinfo', 'yes');
    %eegplot(EEG_miss.data,'srate',EEG_miss.srate,'eloc_file',EEG_miss.chanlocs,'events',EEG_miss.event)
    
    %remove baseline
    EEG_hit = pop_rmbase(EEG_hit, [-199 0]);
    EEG_miss = pop_rmbase(EEG_miss, [-199 0]);
    
    %average data on epoch level
    EEG_hit.data=mean(EEG_hit.data(:,:,:),3);
    EEG_miss.data=mean(EEG_miss.data(:,:,:),3);
    
    %concatenate data for all participants
    if eeg_file==1
        EEG_hit_all=EEG_hit;
        EEG_miss_all= EEG_miss;
    else
        EEG_hit_all.data=cat(3, EEG_hit_all.data, EEG_hit.data);
        EEG_miss_all.data=cat(3, EEG_miss_all.data, EEG_miss.data);
    end
end



% eegplot(EEG_miss_all.data,'srate',EEG_miss_all.srate,'eloc_file',EEG_miss_all.chanlocs)

%% STEP 2
% change shape of data to be applicabel for TFCE 
EEG_hit_all=double(permute(EEG_hit_all.data, [3 1 2]));
EEG_miss_all=double(permute(EEG_miss_all.data, [3 1 2]));

% remove VEOG and HEOG
EEG_hit_all = EEG_hit_all(:,1:70,:);
EEG_miss_all = EEG_miss_all(:,1:70,:);



%% STEP 3
% run t-test tfce
ept_TFCE(EEG_hit_all, EEG_miss_all, EEG_hit.chanlocs(1:70), 'rSample', 512, 'nPerm', 5000, 'type', 'd', 'saveName', 'tfce_hit_miss_delete.mat');


%% STEP 4 
%prepare data for plotting

% calculate difference between natural and manmade images
diff=EEG_hit_all-EEG_miss_all;
%average over all participants
diff_mean=mean(diff(:,:,:),1);
%reshape data for plotting toolbox
mean_permuted=double(permute(diff_mean, [2 3 1]));


%% STEP 5 
% plot results 
% plotting with Ehinger, BV 2018 "EEGVIS toolbox"
%https://github.com/behinger/eegvis/tree/master/topo_butter
addpath(genpath('/net/store/nbp/projects/joint_error/EEG_Belt/TFCE/eegvis'));
load('tfce_hit_miss.mat')
plot_topobutter(mean_permuted ,EEG_hit.times,EEG_hit.chanlocs,'pvalues',Results.P_Values, 'colormap',{{'div','RdYlBu'},{'seq','YlGnBu'},'seq'},'topoalpha',0.05)

% ept_TFCE_Toolbox
% https://www.frontiersin.org/articles/10.3389/fpsyg.2019.00361/full


%% STEP 6 For Result section

load('tfce_hit_miss.mat')
% or
cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS/derivatives
%add TFCE toolbox 
addpath(genpath('/net/store/nbp/projects/joint_error/EEG_Belt/TFCE/ept_TFCE-matlab-master'));
ept_TFCE_Toolbox

% save csv files to check statistics 
p_values_hypo3a = Results.P_Values(:,:)
csvwrite('p_values_hypo3a.csv',p_values_hypo3a)

