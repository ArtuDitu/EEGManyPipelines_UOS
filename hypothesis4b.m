%Matlab version: R2020a
% eeglab version: 2020_0

%addpath D:\Dropbox\fieldtrip-20220426
%ft_defaults
%addpath(genpath('D:\Dropbox\MATLAB_Tools'))
%change the path to the folder with scripts and add the folder to the path
cd D:\Dropbox\Projects\EEGManyPipelines\Data\final\
addpath(genpath('D:\Dropbox\Projects\EEGManyPipelines'))
eeglab;

list_of_files = dir('**/final*');

for eeg_file = 1:size(list_of_files)
    %load file
    cd D:\Dropbox\Projects\EEGManyPipelines\Data\final\
    EEG = pop_loadset(list_of_files(eeg_file).name);
    EEG = pop_selectevent( EEG, 'type',[1010 1020 1030 1040 1090 1110 1120 1190 2010 2020 2030 2040 2090 2110 2120 2190] ,'renametype','forgotten','deleteevents','off','deleteepochs','off','invertepochs','off');
    EEG = pop_selectevent( EEG, 'type',[1011 1021 1031 1041 1091 1111 1121 1191 2011 2021 2031 2041 2091 2111 2121 2191] ,'renametype','remembered','deleteevents','off','deleteepochs','off','invertepochs','off');
    window_epoch=[-.8, 1.3];
    triggers = {'forgotten', 'remembered'};
    EEG = pop_epoch( EEG, triggers, window_epoch, 'epochinfo', 'yes');
    
    cd D:\Dropbox\Projects\EEGManyPipelines\Data\4b\
    EEG = pop_saveset(EEG, 'filename',EEG.filename);
end

ALLEEG = [];
list_of_files = dir('*.set*');
for i = 1:size(list_of_files)
 
    EEG = pop_loadset(list_of_files(i).name), 'loadmode', 'info'; 
    % Store the current EEG to ALLEEG.
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);
end

STUDY=[];
[STUDY ALLEEG]=std_editset(STUDY, ALLEEG, 'name', 'ERSP_Hypothesis4b',...
    'task', '4b', 'filename', 'ERSP_Hypothesis4b.study', 'filepath', ...
   'D:\Dropbox\Projects\EEGManyPipelines\Data\4b\');

STUDY = std_makedesign(STUDY, ALLEEG, 1, 'name','4b','delfiles','off','defaultdesign','off','variable1','type','values1',{'forgotten','remembered'},'vartype1','categorical','subjselect',{'sub-001','sub-002','sub-003','sub-004','sub-005','sub-006','sub-007','sub-008','sub-009','sub-010','sub-011','sub-012','sub-013','sub-014','sub-015','sub-016','sub-017','sub-018','sub-019','sub-020','sub-021','sub-022','sub-023','sub-024','sub-025','sub-026','sub-027','sub-028','sub-029','sub-030','sub-031','sub-032','sub-033'});
[STUDY, ALLEEG] = std_precomp(STUDY, ALLEEG, 'channels', 'savetrials','on','interp','off','erp', 'off', 'ersp','on','erspparams',{'cycles',[3 0.8],'nfreqs',100, 'ntimesout',60,'freqs', [3 40]});

STUDY = pop_statparams(STUDY, 'condstats', 'on', 'method', 'perm');
for i =1:70
    std_erspplot(STUDY, ALLEEG, 'channels', {ALLEEG(1).chanlocs(1:70).labels}, 'design', 1, 'plotconditions', 'together', 'timerange', [-200 800], 'alpha',0.05, 'mcorrect', 'fdr');
    annotation('textbox', [0.45, 1, 0, 0], 'string', ALLEEG(1).chanlocs(i).labels);
end