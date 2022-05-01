addpath D:\Dropbox\fieldtrip-20220426
ft_defaults
addpath(genpath('D:\Dropbox\MATLAB_Tools'))
%change the path to the folder with scripts and add the folder to the path
cd D:\Dropbox\Projects\EEGManyPipelines\Data\derivatives\
addpath(genpath('D:\Dropbox\Projects\EEGManyPipelines'))





list_of_files = dir('**/final*');


for eeg_file = 28:size(list_of_files)
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



% history
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
[STUDY ALLEEG] = std_editset( STUDY, [], 'name','Hypothesis3b','commands',{{'index',1,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-001_task-xxxx_eeg.set'},{'index',2,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-002_task-xxxx_eeg.set'},{'index',3,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-003_task-xxxx_eeg.set'},{'index',4,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-004_task-xxxx_eeg.set'},{'index',5,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-005_task-xxxx_eeg.set'},{'index',6,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-006_task-xxxx_eeg.set'},{'index',7,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-007_task-xxxx_eeg.set'},{'index',8,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-008_task-xxxx_eeg.set'},{'index',9,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-009_task-xxxx_eeg.set'},{'index',10,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-010_task-xxxx_eeg.set'},{'index',11,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-011_task-xxxx_eeg.set'},{'index',12,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-012_task-xxxx_eeg.set'},{'index',13,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-013_task-xxxx_eeg.set'},{'index',14,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-014_task-xxxx_eeg.set'},{'index',15,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-015_task-xxxx_eeg.set'},{'index',16,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-016_task-xxxx_eeg.set'},{'index',17,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-017_task-xxxx_eeg.set'},{'index',18,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-018_task-xxxx_eeg.set'},{'index',19,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-019_task-xxxx_eeg.set'},{'index',20,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-020_task-xxxx_eeg.set'},{'index',21,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-021_task-xxxx_eeg.set'},{'index',22,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-022_task-xxxx_eeg.set'},{'index',23,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-023_task-xxxx_eeg.set'},{'index',24,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-024_task-xxxx_eeg.set'},{'index',25,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-025_task-xxxx_eeg.set'},{'index',26,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-026_task-xxxx_eeg.set'},{'index',27,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-027_task-xxxx_eeg.set'},{'index',28,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-028_task-xxxx_eeg.set'},{'index',29,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-029_task-xxxx_eeg.set'},{'index',30,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-030_task-xxxx_eeg.set'},{'index',31,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-031_task-xxxx_eeg.set'},{'index',32,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-032_task-xxxx_eeg.set'},{'index',33,'load','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\derivatives\\final_sub-033_task-xxxx_eeg.set'}},'updatedat','on','savedat','on','rmclust','on' );
[STUDY ALLEEG] = std_checkset(STUDY, ALLEEG);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
[STUDY, ALLEEG] = std_precomp(STUDY, ALLEEG, {},'savetrials','on','interp','on','recompute','on','ersp','on','erspparams',{'cycles',[3 0.8] ,'nfreqs',30,'ntimesout',60},'itc','on');
STUDY = std_makedesign(STUDY, ALLEEG, 1, 'name','STUDY.design 1','delfiles','off','defaultdesign','off','variable1','type','values1',{'hit','miss'},'vartype1','categorical','subjselect',{'sub-001','sub-002','sub-003','sub-004','sub-005','sub-006','sub-007','sub-008','sub-009','sub-010','sub-011','sub-012','sub-013','sub-014','sub-015','sub-016','sub-017','sub-018','sub-019','sub-020','sub-021','sub-022','sub-023','sub-024','sub-025','sub-026','sub-027','sub-028','sub-029','sub-030','sub-031','sub-032','sub-033'});
EEG = eeg_checkset( EEG );
[STUDY EEG] = pop_savestudy( STUDY, EEG, 'filename','test_hit_miss.study','filepath','D:\\Dropbox\\Projects\\EEGManyPipelines\\Data\\');
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
pop_limo(STUDY, ALLEEG, 'method','WLS','measure','dattimef','freqlim',[1 25] ,'erase','on','splitreg','off','interaction','off');

