addpath(genpath('D:\Dropbox\MATLAB_Tools'))
%change the path to the folder with scripts and add the folder to the path
cd D:\Dropbox\Projects\EEGManyPipelines\Data\final
addpath(genpath('D:\Dropbox\Projects\EEGManyPipelines\Data'))


list_of_files = dir('**/final*');


for eeg_file = 1:size(list_of_files)
    %load file
    cd D:\Dropbox\Projects\EEGManyPipelines\Data\derivatives
    EEG = pop_loadset(list_of_files(eeg_file).name);
    EEG = pop_selectevent( EEG, 'type',[1110 1111 1119 2110 2111 2119] ,'renametype','hit','deleteevents','off','deleteepochs','off','invertepochs','off');
    EEG = pop_selectevent( EEG, 'type',[1120 1121 1129 2120 2121 2129] ,'renametype','miss','deleteevents','off','deleteepochs','off','invertepochs','off');
    window_epoch=[-.8, 1.3];
    triggers = {'hit', 'miss'};
    EEG = pop_epoch( EEG, triggers, window_epoch, 'epochinfo', 'yes');
    EEG = pop_saveset(EEG, 'filename',EEG.filename);
end