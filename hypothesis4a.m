%Matlab version: R2020a
% eeglab version: 2020_0
% TFCE toolbox: https://github.com/Mensen/ept_TFCE-matlab

%Hypothesis 4a
%There are effects of subsequent memory (i.e., a difference between images that will
%be successfully remembered vs. forgotten on a subsequent repetition) ...
%a.on EEG voltage at any channels, at any time.

cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS % path to data
addpath(genpath('/net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines')) % add your folder to the path
eeglabpath = fileparts(which('eeglab.m'));

eeglab; % start toolbox

%add TFCE toolbox 
addpath(genpath('/net/store/nbp/projects/joint_error/EEG_Belt/TFCE/ept_TFCE-matlab-master'));
%add triggers
% 3rd digit 1 and 2
% 0 and 1 last digit
trigger_forg = {1010 1020 1030 1040 1090 1110 1120 1190 2010 2020 2030 2040 2090 2110 2120 2190};
trigger_rem = {1011 1021 1031 1041 1091 1111 1121 1191 2011 2021 2031 2041 2091 2111 2121 2191};
window_epoch = [-.2 .8]; 

%load data

cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS/derivatives % path to data, 
addpath(genpath('/net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS/derivatives')) 
list_of_files = dir('**/final*'); %save all cleaned data sets in one file

%% STEP 1

for eeg_file = 1:size(list_of_files)
    %load file
    cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS/derivatives
    EEG = pop_loadset(list_of_files(eeg_file).name);
    %epoch data
    EEG_forg = pop_epoch( EEG, trigger_forg, window_epoch, 'epochinfo', 'yes');
    EEG_rem = pop_epoch( EEG, trigger_rem, window_epoch, 'epochinfo', 'yes');
    %eegplot(EEG_rem.data,'srate',EEG_rem.srate,'eloc_file',EEG_rem.chanlocs,'events',EEG_rem.event)
    
    %remove baseline
    EEG_forg = pop_rmbase(EEG_forg, [-199 0]);
    EEG_rem = pop_rmbase(EEG_rem, [-199 0]);
    
    %average data on epoch level
    EEG_forg.data=mean(EEG_forg.data(:,:,:),3);
    EEG_rem.data=mean(EEG_rem.data(:,:,:),3);
    
    %concatenate data for all participants
    if eeg_file==1
        EEG_forg_all=EEG_forg;
        EEG_rem_all= EEG_rem;
    else
        EEG_forg_all.data=cat(3, EEG_forg_all.data, EEG_forg.data);
        EEG_rem_all.data=cat(3, EEG_rem_all.data, EEG_rem.data);
    end
end


%eegplot(EEG_rem_all.data,'srate',EEG_rem_all.srate,'eloc_file',EEG_rem_all.chanlocs)

%% STEP 2
% change shape of data to be applicabel for TFCE 
EEG_forg_all=double(permute(EEG_forg_all.data, [3 1 2]));
EEG_rem_all=double(permute(EEG_rem_all.data, [3 1 2]));

% remove VEOG and HEOG
EEG_forg_all = EEG_forg_all(:,1:70,:);
EEG_rem_all = EEG_rem_all(:,1:70,:);



%% STEP 3
% run t-test tfce
ept_TFCE(EEG_forg_all, EEG_rem_all, EEG_forg.chanlocs(1:70), 'rSample', 512, 'nPerm', 5000, 'type', 'd', 'saveName', 'tfce_forg_rem.mat');


%% STEP 4 
%prepare data for plotting

% calculate difference between natural and manmade images
diff=EEG_forg_all-EEG_rem_all;
%average over all participants
diff_mean=mean(diff(:,:,:),1);
%reshape data for plotting toolbox
mean_permuted=double(permute(diff_mean, [2 3 1]));


%% STEP 5 
% plot results 
% plotting with Ehinger, BV 2018 "EEGVIS toolbox"
%https://github.com/behinger/eegvis/tree/master/topo_butter
addpath(genpath('/net/store/nbp/projects/joint_error/EEG_Belt/TFCE/eegvis'));
load('tfce_forg_rem.mat')
plot_topobutter(mean_permuted ,EEG_forg.times,EEG_forg.chanlocs,'pvalues',Results.P_Values, 'colormap',{{'div','RdYlBu'},{'seq','YlGnBu'},'seq'},'topoalpha',0.05)

% ept_TFCE_Toolbox
% https://www.frontiersin.org/articles/10.3389/fpsyg.2019.00361/full

%% for result

load('tfce_forg_rem.mat')
% or
cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS/derivatives
%add TFCE toolbox 
addpath(genpath('/net/store/nbp/projects/joint_error/EEG_Belt/TFCE/ept_TFCE-matlab-master'));
ept_TFCE_Toolbox


% save csv files to check statistics 
p_values_hypo4a = Results.P_Values(:,:);
csvwrite('p_values_hypo4a.csv',p_values_hypo4a)

