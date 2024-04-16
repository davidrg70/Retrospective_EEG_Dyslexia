% dg_SaveProcStamps_afterNewProc_gluedEEGs

% This script opens conf/subject_cache, copies the subject.stamps.processing_ANTS9Years3T field, opens the
% processed/selected trials, copies this field, saves it, and removes the subject.stamps.processing_nmriSUMAcanon field.
% This allows proper export with the ANTS9Years3T headmodel, and not the old nmriSUMA headmodel, because export needs
% this newer stamp as well as a newer a subject.loaded_from field (Vigilance!)

clear all; close all;
analysis_dir = 'group dir';
subject_list = {'pseudonyms'}; % (Previously "glued" and INCLUDED EEGs)

cd(analysis_dir);
% addpath(genpath(fullfile(analysis_dir, '/scripts')));
analysis_params; %load params

structs = {'subject_cache', 'selected_trials_'}; % list of structs to correct

for i = 1:length(subject_list)
    clear subject data
    folderpath = fullfile(analysis_dir, subject_list{i}, '/eeg/'); % folderpath to explore data in "eeg" folder
    folderpath = fullfile(folderpath, '**');
    filelist   = dir(folderpath);
    name       = {filelist.name};
    name       = name(~strncmp(name, '.', 1));
    
    eeg_num1 = 'EEG.1'; eeg_num2 = 'EEG.2'; eeg_num3 = 'EEG.3';
    
    % if-statement to call the correct subject structure (number 1, 2, or 3)
    if contains(name(1), eeg_num1)
        structure = 'eeg_EEG_1'; substruct = '.EEG.1';
    elseif contains(name(1), eeg_num2)
        structure = 'eeg_EEG_2'; substruct = '.EEG.2';
    elseif contains(name(1), eeg_num3)
        structure = 'eeg_EEG_3'; substruct = '.EEG.3';
    else
        error('No EEG number identified, it could be EEG.4 or superior')
    end
    
    % open conf/subject_cache
    folderpath = fullfile(analysis_dir, subject_list{i}, '/conf/', structure); % folderpath re-named to call "processed" folder
    file = fullfile(folderpath, [structs{1}, '.mat']);
    [~, filename1,~] = fileparts(file);
    load(file);
    
    % save new processing stamp and loaded_from fields
    new_proc_stamp = subject.stamps.processing_ANTS9Years3T_fs7_openmeeg;
    new_loaded_from = subject.loaded_from;
    
    clear subject data max_timestamp source timestamp
    % open processed/selected trials_
    folderpath = fullfile(analysis_dir, subject_list{i}, '/processed/');
    file = fullfile(folderpath, [structs{2}, subject_list{i}, '_', structure, '.mat']);
    [~, filename2,~] = fileparts(file);
    load(file);
    
    % paste fields
    subject.loaded_from = new_loaded_from;
    subject.stamps.processing_ANTS9Years3T_fs7_openmeeg = new_proc_stamp;
    
    % remove previous stamps
    if isfield(subject.stamps, 'processing_nmriSUMA150_fs6_openmeeg')
        subject.stamps = rmfield(subject.stamps, 'processing_nmriSUMA150_fs6_openmeeg');
    end
    
    if isfield(subject.stamps, 'hdmSUMA_nmriSUMA150_fs6_dipoli')
        subject.stamps = rmfield(subject.stamps, 'hdmSUMA_nmriSUMA150_fs6_dipoli');
    end
    % save changes
    save(file, 'subject');
    
    fprintf('Stamps updated (%s ) \n', filename2);
end