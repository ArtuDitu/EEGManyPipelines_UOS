%% 64 PREPROCESSING

%% STEP ONE: Load data + Filtering (high-low), Resampling

% add paths and create a folder
cd D:\Dropbox\Projects\EEGManyPipelines % change directory
addpath(genpath('D:\Dropbox\Projects\EEGManyPipelines')) % add directory with all subfolders to the path
eeglabpath = fileparts(which('eeglab.m')); % create variable storing a path to eeglab toolbox
eeglab;

%EEG = pop_loadset(); % load data: select data from your disk
EEG = pop_loadset('filename', 'EMP01.set', 'filepath','D:\Dropbox\Projects\EEGManyPipelines\eeg_eeglab' ); % load data: specify file names. This can be further improved with sprintf() to loop over different files

% parameters
low_pass = 100;
high_pass = .1;
power_line = 50;

EEG.preprocessing = []; % add field to keep track of changes

%%% DOWNSAMPLING %%%
%EEG = pop_resample(EEG, 512);
%EEG.preprocessing = [EEG.preprocessing 'Resampled,'];

%pop_spectopo(EEG)

%%% HighPass FILTER %%%
EEG = pop_eegfiltnew(EEG, high_pass, []);   % highpass
EEG.preprocessing = [EEG.preprocessing 'Highpass,'];

%%% LowPass FILTER %%%
EEG = pop_eegfiltnew(EEG, [], low_pass);   % lowpass
EEG.preprocessing = [EEG.preprocessing 'Lowpass,'];

%pop_spectopo(EEG)

% remove line noise with zapline 
d_tmp = permute(EEG.data, [2,1]); % change dimensions 
d_tmp = nt_zapline(d_tmp, power_line/EEG.srate); % use zapline
EEG.data = permute(d_tmp,[2,1]); % change dimensions back
EEG.preprocessing = [EEG.preprocessing 'Zaplined,'];

%pop_spectopo(EEG)

% referenca data to average
EEG = pop_reref( EEG, []);

% saving the new dataset
EEG = pop_editset(EEG, 'setname', sprintf('Step1_%s',EEG.setname));
EEG = pop_saveset(EEG, 'filename',EEG.setname,'filepath',fullfile(EEG.filepath, 'preprocessed'));

%% STEP TWO: Selecting channels 

% plot data to scroll and find noisy channels
eegplot(EEG.data,'srate',EEG.srate,'eloc_file',EEG.chanlocs,'events',EEG.event)

% here you enter names of channels to reject
channels_to_reject = {'T8' 'T7'};

%deleting noisy channels
EEG = pop_select(EEG, 'nochannel', channels_to_reject);
EEG.preprocessing = [EEG.preprocessing 'ChannelReject,'];

% save the list of rejected channels
save([EEG.filepath sprintf('\\%s_channels_to_reject.mat',EEG.setname(7:end))], 'channels_to_reject')

% rereference after removing noisy channels
if EEG.nbchan ~= 72
    EEG = pop_reref( EEG, []); 
end

%save file
EEG = pop_editset(EEG, 'setname', sprintf('Step2_%s',EEG.setname(7:end)));
EEG = pop_saveset(EEG, 'filename',EEG.setname,'filepath',fullfile(EEG.filepath,'preprocessed'));
%% STEP THREE: Data Cleaning

% plot data to scroll and reject chunks
eegplot(EEG.data,'command','rej=TMPREJ;','srate',EEG.srate,'eloc_file',EEG.chanlocs,'events',EEG.event);

%save cleaning file
save([EEG.filepath,sprintf('\\%s_cleaningTimes.mat', EEG.setname(7:end))],'tmprej','rej');

%converts rejections into events
tmprej = eegplot2event(rej, -1);
EEG = eeg_eegrej(EEG,tmprej(:,[3 4]));
EEG.preprocessing = [EEG.preprocessing 'Cleaning,'];

% eeg check function after cleaning
EEG = eeg_checkset(EEG,'makeur');

%save file
EEG = pop_editset(EEG, 'setname', sprintf('Step3_%s',EEG.setname(7:end)));
EEG = pop_saveset(EEG, 'filename',EEG.setname,'filepath',fullfile(EEG.filepath));
%% STEP FOUR: ICA (o.o)

%%% HighPass FILTER %%% tmp variable
EEG_tmp = pop_eegfiltnew(EEG, 2, []);   % highpass  2hz to remove slow ocular drifts

% addpath('./amica')
mkdir(fullfile(EEG.filepath,EEG.setname(7:end),'amica'))

% run amica ICA
outDir = what('amica');
dataRank = rank(double(EEG_tmp.data'));

%outDir = fullfile(fullfile(EEG.filepath,EEG.setname(7:end),'amica'));

runamica15(EEG_tmp.data, 'num_chans', EEG_tmp.nbchan,...
            'outdir', outDir.path,...
            'pcakeep', dataRank, 'num_models', 1,...
            'do_reject', 1, 'numrej', 15, 'rejsig', 3, 'rejint', 1);

%% STEP FIVE: Apply ICA wights and Component cleaning


addpath([fullfile(EEG.filepath,'amica')])

%load ICA results
outDir = what('amica');

% load ICA weights
mod = loadmodout15(outDir.path);
         
%apply ICA weights to data
EEG.icasphere = mod.S;
EEG.icaweights = mod.W;
EEG = eeg_checkset(EEG);
EEG.preprocessing = [EEG.preprocessing 'AMICA,'];

%epoch data
window = [-.2 .6];
triggers = {1010 1011 1019 1020 1021 1029 1030 1031 1039 1040 1041 1049 1090 1091 1099 1110 1111 1119 1120 1121 1129 1190 1191 1199 2010 2011 2019 2020 2021 2029 2030 2031 2039 2040 2041 2049 2090 2091 2099 2110 2111 2119 2120 2121 2129 2190 2191 2199}; 
EEG = pop_epoch( EEG,triggers, window, 'epochinfo', 'yes');

% remove baseline
EEG = pop_rmbase(EEG, [-199 0]);

% calculate iclabel classification
EEG = iclabel(EEG);

pop_selectcomps(EEG);
pop_viewprops(EEG, 0)

pop_eegplot(EEG,0,1,1);
                
%reject selected components
comps_to_rej = find(EEG.reject.gcompreject);
EEG = pop_subcomp( EEG, comps_to_rej, 0);
EEG = eeg_checkset(EEG);
EEG.preprocessing = [EEG.preprocessing 'ICACleaned,'];

%save file
EEG = pop_editset(EEG, 'setname', sprintf('Step4_%s',EEG.setname(7:end)));
EEG = pop_saveset(EEG, 'filename',EEG.setname,'filepath',fullfile(EEG.filepath));
save([EEG.filepath,sprintf('\\ICA_%s.mat',EEG.setname)],'comps_to_rej');

%% STEP SIX: Re-reference & Interpolation

load('preprocessed_channels.mat');

%interpolate missing channels
EEG= pop_interp(EEG, full_channels_locs,'spherical');
EEG.preprocessing = [EEG.preprocessing 'channelInterpol'];

%eegplot(EEG.data,'command','rej=TMPREJ;','srate',EEG.srate,'eloc_file',EEG.chanlocs,'events',EEG.event);
%tmprej = eegplot2event(rej, -1);
%EEG = eeg_eegrej(EEG,tmprej(:,[3 4]));
%EEG = eeg_checkset(EEG,'makeur');


% save
EEG = pop_editset(EEG, 'setname', sprintf('Step5_%s',EEG.setname(7:end)));
EEG = pop_saveset(EEG, 'filename',EEG.setname,'filepath',fullfile(EEG.filepath));