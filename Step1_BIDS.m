% BIDS extension for step one
%% STEP ONE: Load data + Filtering (high-low), Resampling

% add paths and create a folder
cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS % change directory
addpath(genpath('/net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS')) % add directory with all subfolders to the path
eeglabpath = fileparts(which('eeglab.m')); % create variable storing a path to eeglab toolbox
eeglab;

%two plugins have to be installed before loading the dataset is possible:
%bva.io to load Brain Vision Analyzer EEG data files
%bids_tools.io to manage bids files 
cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS
[STUDY ALLEEG] = pop_importbids('/net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS'); % load data: select data from your disk
CURRENTSTUDY = 1;
EEG = ALLEEG; 
CURRENTSET = [1:length(EEG)];

%download chanlocs_ced.txt and add it to BIDS folder
eloc=readlocs('chanlocs_ced.txt', 'filetype', 'chanedit'); % load channel locations
% parameters
low_pass = 100;
high_pass = .1;
power_line = 50;


for i=1:length(EEG)
    EEG(i).preprocessing = []; % add field to keep track of changes
    EEG(i)=eeg_checkset(EEG(i), 'loaddata'); %load data for each dataset
    EEG(i).chanlocs=eloc; %add chanlocations to each dataset
    %%% DOWNSAMPLING %%%
    %EEG = pop_resample(EEG, 512);
    %EEG.preprocessing = [EEG.preprocessing 'Resampled,'];

    %pop_spectopo(EEG)

    %%% HighPass FILTER %%%
    EEG(i) = pop_eegfiltnew(EEG(i), high_pass, []);   % highpass
    EEG(i).preprocessing = [EEG(i).preprocessing 'Highpass,'];

    %%% LowPass FILTER %%%
    EEG(i) = pop_eegfiltnew(EEG(i), [], low_pass);   % lowpass
    EEG(i).preprocessing = [EEG(i).preprocessing 'Lowpass,'];

    %pop_spectopo(EEG)

    % remove line noise with zapline 
    d_tmp = permute(EEG(i).data, [2,1]); % change dimensions 
    d_tmp = nt_zapline(d_tmp, power_line/EEG(i).srate); % use zapline
    EEG(i).data = permute(d_tmp,[2,1]); % change dimensions back
    EEG(i).preprocessing = [EEG(i).preprocessing 'Zaplined,'];

    %pop_spectopo(EEG(i))

    % referenc data to average
    EEG(i) = pop_reref( EEG(i), []);

    EEG(i) = pop_editset(EEG(i), 'setname', sprintf('Step1_%s',EEG(i).setname));
    EEG(i).saved='no';
    EEG(i)=pop_saveset(EEG(i), 'savemode', 'resave');
end    
% saving the new dataset


%cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS/derivatives
[STUDY EEG]=pop_savestudy(STUDY, EEG, 'filename', EEG(i).setname);
