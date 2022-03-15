%Hypothesis 1
%There is an effect of scence category on the amplitude of the N1
%component.

cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS % path where I have data, you should change it to the path on your PC
addpath(genpath('/net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS')) % add your folder to the path
eeglabpath = fileparts(which('eeglab.m'));
% load eeglab (it is a toolbox that we use to work with eeg data, you can
% download it online, and add it in the same folder as above)
eeglab; % start toolbox

%add TFCE toolbox 
addpath(genpath('/net/store/nbp/projects/joint_error/EEG_Belt/TFCE/ept_TFCE-matlab-master'));
%add triggers
trigger_manmade = {1010 1011 1019 1020 1021 1029 1030 1031 1039 1040 1041 1049 1090 1091 1099 1110 1111 1119 1120 1121 1129 1190 1191 1199}; 
trigger_natural={2010 2011 2019 2020 2021 2029 2030 2031 2039 2040 2041 2049 2090 2091 2099 2110 2111 2119 2120 2121 2129 2190 2191 2199};
window_epoch = [-.2 .2];

%load data

cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS/derivatives % path where I have data, you should change it to the path on your PC
addpath(genpath('/net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS/derivatives')) % add your folder to the path
list_of_files = dir('**/final*'); %save all cleaned data sets in one file

%% STEP 1
% create two datasets, one which contains only man made scenery and one
% which contains only natural scenery, epoch [-100,200]
for eeg_file = 1:size(list_of_files)
    %load file
    cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS/derivatives
    EEG = pop_loadset(list_of_files(eeg_file).name);
    %epoch data
    EEG_natural = pop_epoch( EEG, trigger_natural, window_epoch, 'epochinfo', 'yes');
    EEG_manmade = pop_epoch( EEG, trigger_manmade, window_epoch, 'epochinfo', 'yes');
    %eegplot(EEG_manmade.data,'srate',EEG_manmade.srate,'eloc_file',EEG_manmade.chanlocs,'events',EEG_manmade.event)
    
    %remove baseline
    EEG_natural = pop_rmbase(EEG_natural, [-199 0]);
    EEG_manmade = pop_rmbase(EEG_manmade, [-199 0]);
    
    %average data on epoch level
    EEG_natural.data=mean(EEG_natural.data(:,:,:),3);
    EEG_manmade.data=mean(EEG_manmade.data(:,:,:),3);
    
    %concatenate data for all participants
    if eeg_file==1
        EEG_natural_all=EEG_natural;
        EEG_manmade_all= EEG_manmade;
    else
        EEG_natural_all.data=cat(3, EEG_natural_all.data, EEG_natural.data);
        EEG_manmade_all.data=cat(3, EEG_manmade_all.data, EEG_manmade.data);
    end
end

%% STEP 2
% change shape of data to be applicabel for TFCE 
EEG_natural_all=double(permute(EEG_natural_all.data, [3 1 2]));
EEG_manmade_all=double(permute(EEG_manmade_all.data, [3 1 2]));

% remove VEOG and HEOG
EEG_natural_all = EEG_natural_all(:,1:70,:);
EEG_manmade_all = EEG_manmade_all(:,1:70,:);



%% STEP 3
% run t-test tfce
ept_TFCE(EEG_natural_all, EEG_manmade_all, EEG_natural.chanlocs(1:70), 'rSample', 512, 'nPerm', 5000, 'type', 'd', 'saveName', 'tfce_natural_mammade.mat');


%% STEP 4 
%prepare data for plotting

% calculate difference between natural and manmade images
diff=EEG_natural_all-EEG_manmade_all;
%average over all participants
diff_mean=mean(diff(:,:,:),1);
%reshape data for plotting toolbox
mean_permuted=double(permute(diff_mean, [2 3 1]));


%% STEP 5 
% plot results 
% plotting with Ehinger, BV 2018 "EEGVIS toolbox"
%https://github.com/behinger/eegvis/tree/master/topo_butter
addpath(genpath('/net/store/nbp/projects/joint_error/EEG_Belt/TFCE/eegvis'));
load('tfce_natural_mammade.mat')
plot_topobutter(mean_permuted ,EEG_natural.times,EEG_natural.chanlocs,'pvalues',Results.P_Values, 'colormap',{{'div','RdYlBu'},{'seq','YlGnBu'},'seq'},'topoalpha',0.05)

% ept_TFCE_Toolbox
