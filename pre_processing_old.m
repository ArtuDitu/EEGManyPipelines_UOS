%% 64 PREPROCESSING

%% STEP ONE: Load data + Filtering (high-low), Resampling, Deblanking

[filepath filename setname eeglabpath sub] = acz_generate_paths();

%create a target folder
mkdir(fullfile(filepath,'preprocessed'))

if isunix
    acz_EEG = pop_loadeep(fullfile(filepath,filename), 'triggerfile','on');
elseif getenv('computername') == 'DESKTOP-RB7CB67'
    acz_EEG = pop_loadeep_v4(fullfile(filepath,filename), 'triggerfile','on');
end

acz_EEG.preprocessing = [];

%%% DOWNSAMPLING %%%
acz_EEG = pop_resample(acz_EEG, 512);
acz_EEG.preprocessing = [acz_EEG.preprocessing 'Resampled,'];

%%% HighPass FILTER %%%
acz_EEG = pop_eegfiltnew(acz_EEG, 0.1, []);   % highpass  
acz_EEG.preprocessing = [acz_EEG.preprocessing 'Highpass,'];
               
%%% LowPass FILTER %%%
acz_EEG = pop_eegfiltnew(acz_EEG, [], 120);   % lowpass
acz_EEG.preprocessing = [acz_EEG.preprocessing 'Lowpass,'];

% look up location of channels and change location of VEOG from right to
% left
acz_EEG=pop_chanedit(acz_EEG, 'lookup',fullfile(eeglabpath,'plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp'));
acz_EEG.chanlocs(65).theta = -27;
pop_spectopo(acz_EEG) % pwelsh?
%pwelch(acz_EEG.data(1,:))

%%% reduce line noise in the data (note: may not completely eliminate, re-referencing helps at the end as well)
acz_EEG = pop_cleanline(acz_EEG, 'bandwidth',2,'chanlist',1:acz_EEG.nbchan,'computepower',1,'linefreqs',[50 100] ,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',0,'winsize',4,'winstep',1, 'ComputeSpectralPower','False');
acz_EEG = eeg_checkset( acz_EEG );
acz_EEG.preprocessing = [acz_EEG.preprocessing 'Cleanlined,'];

%%% DEBLANKING %%%
x=size(acz_EEG.event);
for t=1:max(x)
    acz_EEG.event(t).type = deblank(acz_EEG.event(t).type);
end

acz_EEG.preprocessing = [acz_EEG.preprocessing 'Deblanked,'];

% saving the new dataset
acz_EEG = pop_editset(acz_EEG, 'setname', sprintf('2_%s_bandpass_resample_deblank',setname));
acz_EEG = pop_saveset(acz_EEG, 'filename',sprintf('2_%s_bandpass_resample_deblank',setname),'filepath',fullfile(filepath,'preprocessed'));

%% STEP TWO: Selecting channels 
%loading new dataset
if ~exist('acz_EEG','var')
    [filepath filename setname eeglabpath sub] = acz_generate_paths();
    acz_EEG = pop_loadset(sprintf('2_%s_bandpass_resample_deblank.set',setname),fullfile(filepath,'preprocessed'));
elseif isempty(~strfind(acz_EEG.preprocessing,'Resampled,Highpass,Lowpass,Cleanlined,Deblanked,'))
    [filepath filename setname eeglabpath sub] = acz_generate_paths();
    acz_EEG = pop_loadset(sprintf('2_%s_bandpass_resample_deblank.set',setname),fullfile(filepath,'preprocessed'));
end

%deleting unused channels for I amplifier 
acz_EEG = pop_select(acz_EEG, 'nochannel', {'BIP2' 'BIP3' 'BIP4' 'AUX1' 'AUX2' 'AUX3' 'AUX4'}); % measurements with one amplifier


%renaming BIP1
acz_EEG.chanlocs(65).labels = 'VEOG'; % for I amplifier

% look up location of channels and change location of VEOG from right to
% left
acz_EEG=pop_chanedit(acz_EEG, 'lookup',fullfile(eeglabpath,'plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp'));
acz_EEG.chanlocs(65).theta = -27;

% save to use later before interpolation
full_channels_locs = acz_EEG.chanlocs;
save([filepath sprintf('65_channels.mat',setname)], 'full_channels_locs')



if exist(fullfile(filepath,sprintf('3_%s_channels_to_reject.mat',setname)))
    tmp_channels_to_reject = load(fullfile(filepath,sprintf('3_%s_channels_to_reject.mat',setname)));
    channels_to_reject = tmp_channels_to_reject.channels_to_reject;
else
    %plotting data scroll to visually detect noisy channels
    eegplot(acz_EEG.data,'srate',acz_EEG.srate,'eloc_file',acz_EEG.chanlocs,'events',acz_EEG.event)
    channels_to_reject = {'CPz'};
end


%deleting noisy channels
acz_EEG = pop_select(acz_EEG, 'nochannel', channels_to_reject);


if acz_EEG.nbchan ~= 65
    acz_EEG = pop_reref( acz_EEG, []); %Participants’ averages were then re-referenced to a common average reference. (Rossion & Caharel, 2011)
end

%save file
acz_EEG.preprocessing = [acz_EEG.preprocessing 'ChannelReject,'];
acz_EEG = pop_editset(acz_EEG, 'setname', sprintf('3_%s_channelrej',setname));
acz_EEG = pop_saveset(acz_EEG, 'filename',sprintf('3_%s_channelrej',setname),'filepath',fullfile(filepath,'preprocessed'));
save([filepath sprintf('3_%s_channels_to_reject.mat',setname)], 'channels_to_reject')
%% STEP THREE: Data Cleaning
%load files
if ~exist('acz_EEG','var')
    [filepath filename setname eeglabpath sub] = acz_generate_paths();
    acz_EEG = pop_loadset(sprintf('3_%s_channelrej.set',setname),fullfile(filepath,'preprocessed'));
elseif isempty(~strfind(acz_EEG.preprocessing,'Resampled,Highpass,Lowpass,Cleanlined,Deblanked,ChannelReject,'))
    [filepath filename setname eeglabpath sub] = acz_generate_paths();
    acz_EEG = pop_loadset(sprintf('3_%s_channelrej.set',setname),fullfile(filepath,'preprocessed'));
end


if exist (fullfile(filepath,sprintf('4_%s_cleaningTimes.mat',setname)))
    load(fullfile(filepath,sprintf('4_%s_cleaningTimes.mat',setname)),'tmprej','rej');
else
    %plotting for visual inspection
    eegplot(acz_EEG.data,'command','rej=TMPREJ;','srate',acz_EEG.srate,'eloc_file',acz_EEG.chanlocs,'events',acz_EEG.event);
end

%converts rejections into events
tmprej = eegplot2event(rej, -1);
acz_EEG = eeg_eegrej(acz_EEG,tmprej(:,[3 4]));

acz_EEG = eeg_checkset(acz_EEG,'makeur');

%save file
save([filepath,sprintf('4_%s_cleaningTimes.mat',setname)],'tmprej','rej');
acz_EEG.preprocessing = [acz_EEG.preprocessing 'Cleaning,'];
acz_EEG = pop_editset(acz_EEG, 'setname', sprintf('4_%s_Clean',setname));
acz_EEG = pop_saveset(acz_EEG, 'filename',sprintf('4_%s_Clean',setname),'filepath',fullfile(filepath,'preprocessed'));
%% STEP FOUR: ICA (o.o)
%load files
if ~exist('acz_EEG','var')
    [filepath filename setname eeglabpath sub] = acz_generate_paths();
    acz_EEG = pop_loadset(sprintf('4_%s_Clean.set',setname),fullfile(filepath,'preprocessed'));
elseif isempty(~strfind(acz_EEG.preprocessing,'Resampled,Highpass,Lowpass,Cleanlined,Deblanked,ChannelReject,Cleaning,'))
    [filepath filename setname eeglabpath sub] = acz_generate_paths();
    acz_EEG = pop_loadset(sprintf('4_%s_Clean.set',setname),fullfile(filepath,'preprocessed'));
end


%%% HighPass FILTER %%%
acz_EEG = pop_eegfiltnew(acz_EEG, 2, []);   % highpass  

% addpath('./amica')
mkdir(fullfile(filepath,'preprocessed','amica'))
                
outDir = fullfile(filepath,'preprocessed','amica');
runamica12(acz_EEG.data,'outdir',outDir,'use_queue','nbp.q','qsubname',['ica_VP' num2str(sub)]);
%% STEP FIVE: Apply ICA wights and Component cleaning
%load data
if ~exist('acz_EEG','var')
    [filepath filename setname eeglabpath sub] = acz_generate_paths();
    acz_EEG = pop_loadset(sprintf('4_%s_Clean.set',setname),fullfile(filepath,'preprocessed'));
elseif isempty(~strfind(acz_EEG.preprocessing,'Resampled,Highpass,Lowpass,Cleanlined,Deblanked,ChannelReject,Cleaning,'))
    [filepath filename setname eeglabpath sub] = acz_generate_paths();
    acz_EEG = pop_loadset(sprintf('4_%s_Clean.set',setname),fullfile(filepath,'preprocessed'));
end


addpath([fullfile(filepath,'preprocessed') '/amica'])

%load ICA results
icapath = fullfile(filepath,'preprocessed','amica');

mod = loadmodout12(icapath);
                
%apply ICA weights to data
acz_EEG.icasphere = mod.S;
acz_EEG.icaweights = mod.W;
acz_EEG = eeg_checkset(acz_EEG);
acz_EEG.preprocessing = [acz_EEG.preprocessing 'AMICA,'];

acz_EEG_temp = acz_EEG;

%EPOCH %%%
window=[-0.2 1];
[acz_EEG_temp indices] = pop_epoch( acz_EEG_temp, {4 5 6 7 11 12 13 14}, window, 'epochinfo', 'yes');
acz_EEG_temp.orig_indices = indices;


% remove baseline
acz_EEG_temp = pop_rmbase( acz_EEG_temp, [-200 0]);

acz_EEG_temp = pop_chanedit(acz_EEG_temp, 'lookup',fullfile(eeglabpath,'plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp'));
veog_channel = size(acz_EEG_temp.icachansind,2);
acz_EEG_temp.chanlocs(veog_channel).theta = -27;

EEG = acz_EEG_temp;
pop_selectcomps(EEG);
pop_eegplot(EEG,0,1,1);
                
%reject selected components
comps_to_rej = find(EEG.reject.gcompreject);
acz_EEG = pop_subcomp( acz_EEG, comps_to_rej, 0);

%save file
acz_EEG.preprocessing = [acz_EEG.preprocessing 'ICACleaned,'];
acz_EEG = pop_editset(acz_EEG, 'setname', sprintf('5_%s_ICAclean',setname));
acz_EEG = pop_saveset(acz_EEG, 'filename',sprintf('5_%s_ICAclean',setname),'filepath',fullfile(filepath,'preprocessed'));
save([filepath,sprintf('5_%s_ICA.mat',setname)],'comps_to_rej');
%% STEP SIX: Re-reference & Interpolation
%load files
if ~exist('acz_EEG','var')
    [filepath filename setname eeglabpath sub] = acz_generate_paths();
    acz_EEG = pop_loadset(sprintf('5_%s_ICAclean.set',setname),fullfile(filepath,'preprocessed'));
elseif isempty(~strfind(acz_EEG.preprocessing,'Resampled,Highpass,Lowpass,Cleanlined,Deblanked,ChannelReject,Cleaning,AMICA,ICACleaned,'))
    [filepath filename setname eeglabpath sub] = acz_generate_paths();
    acz_EEG = pop_loadset(sprintf('5_%s_ICAclean.set',setname),fullfile(filepath,'preprocessed'));
end

load(fullfile(filepath,sprintf('65_channels.mat',setname)));


acz_EEG = pop_reref( acz_EEG, []); %Participants’ averages were then re-referenced to a common average reference. (Rossion & Caharel, 2011)
acz_EEG.preprocessing = [acz_EEG.preprocessing 'Rereference,'];

%interpolate missing channels
acz_EEG=pop_chanedit(acz_EEG, 'lookup',fullfile(eeglabpath,'plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp'));
acz_EEG= pop_interp(acz_EEG, full_channels_locs,'spherical');


acz_EEG.preprocessing = [acz_EEG.preprocessing 'channelInterpol'];
acz_EEG = pop_editset(acz_EEG, 'setname', sprintf('6_%s_RerefInterp',setname));
acz_EEG = pop_saveset(acz_EEG, 'filename',sprintf('6_%s_RerefInterp',setname),'filepath',fullfile(filepath,'preprocessed'));

%% STEP SEVEN: epoching and recleaning
%load data
if ~exist('acz_EEG','var')
    [filepath filename setname eeglabpath sub] = acz_generate_paths();
    acz_EEG = pop_loadset(sprintf('6_%s_RerefInterp.set',setname),fullfile(filepath,'preprocessed'));
elseif isempty(~strfind(acz_EEG.preprocessing,'Resampled,Highpass,Lowpass,Cleanlined,Deblanked,ChannelReject,Cleaning,AMICA,ICACleaned,Rereference,channelInterpol'))
    [filepath filename setname eeglabpath sub] = acz_generate_paths();
    acz_EEG = pop_loadset(sprintf('6_%s_RerefInterp.set',setname),fullfile(filepath,'preprocessed'));
end

%EPOCH %%%
window=[-0.2 1];
[acz_EEG indices] = pop_epoch( acz_EEG, {4 5 6 7 11 12 13 14}, window, 'epochinfo', 'yes');
acz_EEG.orig_indices = indices;

% % remove baseline
acz_EEG = pop_rmbase( acz_EEG, [-200 0]);
% 
% % recleaning
eegplot(acz_EEG.data, 'command', 'rej_epoch=TMPREJ','srate',acz_EEG.srate,'eloc_file',acz_EEG.chanlocs, 'events',acz_EEG.event);
% 
% %converts rejections into epoch
tmprej_epoch = eegplot2trial(rej_epoch,acz_EEG.pnts, acz_EEG.trials);
acz_EEG = pop_rejepoch(acz_EEG,tmprej_epoch);
acz_EEG = eeg_checkset(acz_EEG,'makeur');

acz_EEG.preprocessing = [acz_EEG.preprocessing 'EPOCHED'];
acz_EEG = pop_editset(acz_EEG, 'setname', sprintf('7_%s_EPOCHED',setname));
acz_EEG = pop_saveset(acz_EEG, 'filename', sprintf('7_%s_EPOCHED',setname),'filepath',fullfile(filepath,'preprocessed'));
%% STEP EIGHT epoch and average (applied after all files are pre processed!!!!!
x= [ 4 5 6 7 11 12 13 14];
y = [400 403 404 405 406 408 410 413 416 417 418 419 420 422 424 425 426 427 428 432 433 434 435 436 437 439];


for i=y
    % open eeglab and load clean data set
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadset('filename',sprintf('7_%d_EPOCHED.set',i),'filepath',sprintf('/net/store/nbp/projects/joint_error/all/Data/study/EEG/sub%d/preprocessed/',i));
    EEG=pop_chanedit(EEG, 'lookup','/net/store/nbp/projects/joint_error/all/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp');
    EEG = pop_select( EEG,'channel',{'VEOG' 'Fp1' 'Fpz' 'Fp2' 'F7' 'F3' 'Fz' 'F4' 'F8' 'FC5' 'FC1' 'FC2' 'FC6' 'M1' 'T7' 'C3' 'Cz' 'C4' 'T8' 'M2' 'CP5' 'CP1' 'CP2' 'CP6' 'P7' 'P3' 'Pz' 'P4' 'P8' 'POz' 'O1' 'Oz' 'O2' 'AF7' 'AF3' 'AF4' 'AF8' 'F5' 'F1' 'F2' 'F6' 'FC3' 'FCz' 'FC4' 'C5' 'C1' 'C2' 'C6' 'CP3' 'CPz' 'CP4' 'P5' 'P1' 'P2' 'P6' 'PO5' 'PO3' 'PO4' 'PO6' 'FT7' 'FT8' 'TP7' 'TP8' 'PO7' 'PO8'});
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    for k=x
        EEG_epoch = pop_epoch( EEG, {  k  }, [-0.2 1], 'newname', sprintf('Tr_%d_%d_RerefInterp epochs',k,i), 'epochinfo', 'yes');
        %[ALLEEG EEG_epoch CURRENTSET] = pop_newset(ALLEEG, EEG_epoch, 1,'savenew',sprintf('/net/store/nbp/projects/joint_error/all/Data/study/EEG/sub%d/Tr_%d_%d_new.set',i,k,i),'gui','off'); 
        EEG_epoch = eeg_checkset( EEG_epoch );
        [ALLEEG EEG_epoch CURRENTSET] = pop_newset(ALLEEG, EEG_epoch,11,'savenew',sprintf('/net/store/nbp/projects/joint_error/all/Data/study/EEG/sub%d/Tr_%d_%d.set',i,k,i),'gui','off');
    end
end

for k=x
     [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
     for l=y
     EEG = pop_loadset('filename',sprintf('Tr_%d_%d.set',k,l),'filepath',sprintf('/net/store/nbp/projects/joint_error/all/Data/study/EEG/sub%d',l));
     EEG.avgdata = [];
     EEG.sem = [];
     EEG.sd = [];
        for i = 1:size(EEG.data,1)
            for j = 1:size(EEG.data,2)
                EEG.avgdata(i,j) = mean(EEG.data(i,j,:));
                EEG.sem(i,j) = (std(EEG.data(i,j,:))) ./ (sqrt(length(EEG.data(i,j,:))));
                EEG.sd(i,j) = std(EEG.data(i,j,:));
            end
            
        end
     EEG.data = EEG.avgdata
     EEG.trial_number = EEG.trials;
     EEG = pop_editset(EEG, 'setname', sprintf('8_%d_%d_avg',k,l));
     EEG = pop_saveset(EEG, 'filename',sprintf('8_%d_%d_avg',k,l),'filepath',sprintf('/net/store/nbp/projects/joint_error/all/Data/study/EEG/sub%d',l));
     end
end


