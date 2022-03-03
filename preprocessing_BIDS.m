%% BIDS EXTENSION FOR PRE PROCESSING SCRIPT

%% STEP ONE: Load data + Filtering (high-low), Resampling

% add paths and create a folder
cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines % change directory
addpath(genpath('/net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines')) % add directory with all subfolders to the path
eeglabpath = fileparts(which('eeglab.m')); % create variable storing a path to eeglab toolbox
eeglab;

%two plugins have to be installed before loading the dataset is possible:
%bva.io to load Brain Vision Analyzer EEG data files
%bids_tools.io to manage bids files 
%it might be possible that pop_importbids.m has to be adjusted, in my case
%i had to exchange 'File' with 'file' in each exist() function
cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS
[STUDY ALLEEG] = pop_importbids('/net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/Redo'); % load data: select data from your disk
CURRENTSTUDY = 1;
EEG = ALLEEG; 
CURRENTSET = [1:length(EEG)];

eloc=readlocs('chanlocs_ced.txt', 'filetype', 'chanedit'); % load channel locations
% parameters
low_pass = 100;
high_pass = .1;
power_line = 50;


for i=1:length(EEG)
    EEG(i).preprocessing = []; % add field to keep track of changes
    EEG(i)=eeg_checkset(EEG(i), 'loaddata'); %load data for each dataset
    EEG(i).chanlocs=eloc;
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
%% STEP TWO: Selecting channels 
for j=1:length(EEG)
    EEG(j)=eeg_checkset(EEG(j), 'loaddata'); %load data for each dataset
    %EEG(j).chanlocs=eloc;
    % plot data to scroll and find noisy channels
    %eegplot(EEG(j).data,'srate',EEG(j).srate,'eloc_file',EEG(j).chanlocs,'events',EEG(j).event)
   
    %load channels to reject
    path='/net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS/derivatives/%s/eeg/%s_channels_to_reject.mat';
    load(sprintf(path,EEG(j).filename(1:7), EEG(j).filename(1:7)));

    %deleting noisy channels
    EEG(j) = pop_select(EEG(j), 'nochannel', channels_to_reject);
    EEG(j).preprocessing = [EEG(j).preprocessing 'ChannelReject,'];

    % save the list of rejected channels
    %save([EEG(j).filepath sprintf('\\%s_channels_to_reject.mat',EEG(j).filename(1:7))], 'channels_to_reject')

    % rereference after removing noisy channels
    if EEG(j).nbchan ~= 72
        EEG(j) = pop_reref( EEG(j), []); 
    end

    %save file
    EEG(j) = pop_editset(EEG(j), 'setname', sprintf('Step2_%s',EEG(j).setname(7:end)));
    EEG(j).saved='no';
    EEG(j)=pop_saveset(EEG(j), 'savemode', 'resave');
end

[STUDY EEG]=pop_savestudy(STUDY, EEG, 'filename', EEG(j).setname);





%% STEP THREE: Data Cleaning
% BIDS extension
for k=1:length(EEG)
    % plot data to scroll and reject chunks
    EEG(k)=eeg_checkset(EEG(k), 'loaddata'); %load data for each dataset
   %eegplot(EEG(k).data,'command','rej=TMPREJ;','srate',EEG(k).srate,'eloc_file',EEG(k).chanlocs,'events',EEG(k).event);
    %converts rejections into events
    %tmprej = eegplot2event(rej, -1);

    
    
   %load rejected parts
    load(sprintf('%s_cleaningTimes',EEG(k).filename(1:7)));
    %reject 
    EEG(k) = eeg_eegrej(EEG(k),tmprej(:,[3 4]));
    EEG(k).preprocessing = [EEG(k).preprocessing 'Cleaning,'];

    %save cleaning file
    %save([EEG(k).filepath,sprintf('\\%s_cleaningTimes.mat', EEG(k).filename(1:7))],'tmprej','rej');

    % eeg check function after cleaning
    EEG(k) = eeg_checkset(EEG(k),'makeur');

    %save file
    EEG(k) = pop_editset(EEG(k), 'setname', sprintf('Step3_%s',EEG(k).setname(7:end)));
    EEG(k).saved='no';
    EEG(k)=pop_saveset(EEG(k), 'savemode', 'resave');
end
[STUDY EEG]=pop_savestudy(STUDY, EEG, 'filename', EEG(k).setname);

%% STEP FOUR: ICA (o.o)
% The calculated ICA weights can be found in the corresponding subject folder
for l=1:length(EEG)
    %%% HighPass FILTER %%% tmp variable
    EEG_tmp = pop_eegfiltnew(EEG(l), 2, []);   % highpass  2hz to remove slow ocular drifts
    
    %addpath('./amica')
    mkdir(fullfile(EEG.filename(7:11),'amica'));

    % create outDir path
    path='/net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipeline/eeg_BIDS/derivatives/%s/amica'
    participant=EEG(l).filename(1:7);
    outDir = sprintf(path, participant);
    dataRank = rank(double(EEG_tmp.data'));

    %run amica 
    runamica12(EEG_tmp.data, 'num_chans', EEG_tmp.nbchan,...
            'outdir', outDir,...
            'pcakeep', dataRank, 'num_models', 1,...
            'do_reject', 1, 'numrej', 15, 'rejsig', 3, 'rejint', 1);
end


%% STEP FIVE: Apply ICA wights and Component cleaning

addpath(genpath('/net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/Amica'))
for m=1:length(EEG)
    %addpath([fullfile(EEG(m).filepath,'amica')])
    EEG(m)=eeg_checkset(EEG(m), 'loaddata'); 
    %load ICA results
    path='/net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS/derivatives/%s/amica'
    participant=EEG(m).filename(1:7);
    outDir = sprintf(path, participant);


    % load ICA weights
    mod = loadmodout12(outDir);
    
    %apply ICA weights to data
    EEG(m).icasphere = mod.S;
    EEG(m).icaweights = mod.W;
    EEG(m) = eeg_checkset(EEG(m));
    EEG(m).preprocessing = [EEG(m).preprocessing 'AMICA,'];

    %epoch data
    window = [-.2 .8];
    triggers = {1010 1011 1019 1020 1021 1029 1030 1031 1039 1040 1041 1049 1090 1091 1099 1110 1111 1119 1120 1121 1129 1190 1191 1199 2010 2011 2019 2020 2021 2029 2030 2031 2039 2040 2041 2049 2090 2091 2099 2110 2111 2119 2120 2121 2129 2190 2191 2199}; 
    EEG_tmp = pop_epoch( EEG(m),triggers, window, 'epochinfo', 'yes');

    % remove baseline
    EEG_tmp = pop_rmbase(EEG_tmp, [-199 0]);

    % calculate iclabel classification
    EEG_tmp = iclabel(EEG_tmp);

    %pop_selectcomps(EEG_tmp); % open components to reject
    %pop_viewprops(EEG_tmp, 0); % open IC label

    
                
   %load components to reject, which where previously selected 
    tmp_path = sprintf('%s_ICA.mat',EEG(m).filename(1:7));
    comps_to_rej = load(tmp_path);
    
    EEG(m) = eeg_checkset(EEG(m));
    
    %reject selected components
    EEG(m) = pop_subcomp( EEG(m), comps_to_rej.comps_to_rej, 0);
    EEG(m) = eeg_checkset(EEG(m));
    EEG(m).preprocessing = [EEG(m).preprocessing 'ICACleaned,'];

    %save file
    EEG(m).filename=sprintf('ICA_%s',EEG(m).filename); %change filename to prevent overwriting
    STUDY.datasetinfo(m).filename=EEG(m).filename; %change filename also in STUDY 
    EEG(m) = pop_editset(EEG(m), 'setname', 'Step4_');
    %EEG(m) = pop_saveset(EEG(m), 'filename',EEG(m).filename,'filepath',fullfile(EEG(m).filepath));
    %save(sprintf('\\%s_ICA.mat',EEG.filename(1:7)),'comps_to_rej');
end
[STUDY EEG]=pop_savestudy(STUDY, EEG, 'filename', EEG(m).setname);
%% STEP SIX: Re-reference & Interpolation
for n=1:length(EEG)
    %load('preprocessed_channels.mat');
    EEG(n)=eeg_checkset(EEG(n), 'loaddata'); %only necessary if not yet loaded
    %interpolate missing channels
    EEG(n)= pop_interp(EEG(n), EEG(n).chaninfo.removedchans,'spherical');
    EEG(n).preprocessing = [EEG(n).preprocessing 'channelInterpol'];

    %eegplot(EEG(n).data,'command','rej=TMPREJ;','srate',EEG(n).srate,'eloc_file',EEG(n).chanlocs,'events',EEG(n).event);
    %tmprej = eegplot2event(rej, -1);
    %EEG = eeg_eegrej(EEG,tmprej(:,[3 4]));
    %EEG = eeg_checkset(EEG,'makeur');


    % save
    EEG(n).filename=sprintf('final_%s',EEG(n).filename(5:end)); %change filename to not overwrite
    STUDY.datasetinfo(n).filename=EEG(n).filename; %also change filename in STUDY so that the study will find renamed data
    EEG(n) = pop_editset(EEG(n), 'setname', 'Step5_'); %change setname
    %EEG(n) =
    %pop_saveset(EEG(n),'filename',EEG(n).filename,'filepath',fullfile(EEG(n).filepath)); 
    
end
[STUDY EEG]=pop_savestudy(STUDY, EEG, 'filename', EEG(n).setname);
