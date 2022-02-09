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
%it might be possible that pop_importbids.m has to be adjusted, in my case
%i had to exchange 'File' with 'file' in each exist() function
cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS
[STUDY ALLEEG] = pop_importbids('/net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS'); % load data: select data from your disk
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
    % plot data to scroll and find noisy channels
    %eegplot(EEG(j).data,'srate',EEG(j).srate,'eloc_file',EEG(j).chanlocs,'events',EEG(j).event)

   % load channels to reject
    load(sprintf('\\%s_channels_to_reject.mat',EEG(j).filename(1:7)));

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
    
   %eegplot(EEG(k).data,'command','rej=TMPREJ;','srate',EEG(k).srate,'eloc_file',EEG(k).chanlocs,'events',EEG(k).event);
    
   
   %load rejected parts
    load(sprintf('\\%s_cleaningTimes',EEG(k).filename(1:7)));
    %converts rejections into events
    tmprej = eegplot2event(rej, -1);
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
% BIDS extension
for l=1:length(EEG)
    %%% HighPass FILTER %%% tmp variable
    EEG_tmp = pop_eegfiltnew(EEG(l), 2, []);   % highpass  2hz to remove slow ocular drifts

    % addpath('./amica')
    mkdir(fullfile(EEG(l).filepath,EEG(l).filename(1:7),'amica'))

    % run amica ICA
    outDir = what('amica');
    dataRank = rank(double(EEG_tmp.data'));

    %outDir = fullfile(fullfile(EEG.filepath,EEG.setname(7:end),'amica'));

    runamica15(EEG_tmp.data, 'num_chans', EEG_tmp.nbchan,...
            'outdir', outDir.path,...
            'pcakeep', dataRank, 'num_models', 1,...
            'do_reject', 1, 'numrej', 15, 'rejsig', 3, 'rejint', 1);
end


%% STEP FIVE: Apply ICA wights and Component cleaning

% add paths and create a folder
cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines % change directory
addpath(genpath('/net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines')) % add directory with all subfolders to the path
eeglabpath = fileparts(which('eeglab.m')); % create variable storing a path to eeglab toolbox
eeglab;

cd /net/store/nbp/projects/joint_error/EEG_Belt/EEGManyPipelines/eeg_BIDS/derivatives
[STUDY ALLEEG] = pop_loadstudy('Step3_.study')
EEG=ALLEEG;

for m=1:length(EEG)
    addpath([fullfile(EEG(m).filepath,'amica')])

    %load ICA results
    outDir = what('amica');

    % load ICA weights
    mod = loadmodout15(outDir.path);
    
    %apply ICA weights to data
    EEG(m).icasphere = mod.S;
    EEG(m).icaweights = mod.W;
    EEG(m) = eeg_checkset(EEG(m));
    EEG(m).preprocessing = [EEG(m).preprocessing 'AMICA,'];

    %epoch data
    window = [-.2 .6];
    triggers = {1010 1011 1019 1020 1021 1029 1030 1031 1039 1040 1041 1049 1090 1091 1099 1110 1111 1119 1120 1121 1129 1190 1191 1199 2010 2011 2019 2020 2021 2029 2030 2031 2039 2040 2041 2049 2090 2091 2099 2110 2111 2119 2120 2121 2129 2190 2191 2199}; 
    EEG(m) = pop_epoch( EEG(m),triggers, window, 'epochinfo', 'yes');

    % remove baseline
    EEG(m) = pop_rmbase(EEG(m), [-199 0]);

    % calculate iclabel classification
    EEG(m) = iclabel(EEG(m));

    pop_selectcomps(EEG(m)); % open components to reject
    pop_viewprops(EEG(m), 0); % open IC label

    pop_eegplot(EEG(m),0,1,1);
                
    %reject selected components
    comps_to_rej = find(EEG(m).reject.gcompreject);
    EEG(m) = pop_subcomp( EEG(m), comps_to_rej, 0);
    EEG(m) = eeg_checkset(EEG(m));
    EEG(m).preprocessing = [EEG(m).preprocessing 'ICACleaned,'];

    %save file
    EEG(m) = pop_editset(EEG(m), 'setname', sprintf('Step4_%s',EEG(m).setname(1:7)));
    EEG(m) = pop_saveset(EEG(m), 'filename',EEG(m).setname,'filepath',fullfile(EEG(m).filepath));
    save([EEG(m).filepath,sprintf('\\ICA_%s.mat',EEG(m).setname)],'comps_to_rej');
end
[STUDY EEG]=pop_savestudy(STUDY, EEG, 'filename', EEG(m).setname);
%% STEP SIX: Re-reference & Interpolation
for n=1:length(EEG)
    %load('preprocessed_channels.mat');

    %interpolate missing channels
    EEG(n)= pop_interp(EEG(n), full_channels_locs,'spherical');
    EEG(n).preprocessing = [EEG(n).preprocessing 'channelInterpol'];

    eegplot(EEG(n).data,'command','rej=TMPREJ;','srate',EEG(n).srate,'eloc_file',EEG(n).chanlocs,'events',EEG(n).event);
    %tmprej = eegplot2event(rej, -1);
    %EEG = eeg_eegrej(EEG,tmprej(:,[3 4]));
    %EEG = eeg_checkset(EEG,'makeur');


    % save
    EEG(n) = pop_editset(EEG(n), 'setname', sprintf('Step5_%s',EEG(n).setname(1:7)));
    EEG(n) = pop_saveset(EEG(n), 'filename',EEG(n).setname,'filepath',fullfile(EEG(n).filepath));
end
[STUDY EEG]=pop_savestudy(STUDY, EEG, 'filename', EEG(n).setname);
