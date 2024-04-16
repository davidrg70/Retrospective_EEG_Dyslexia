% dg_repair_glued_sampleinfo
clear all; close all;
analysis_dir = 'group dir';
subject_list = {'pseudonyms'}; % (Previously "glued" EEGs)

cd(analysis_dir);
% addpath(genpath(fullfile(analysis_dir, '/scripts')));
analysis_params; %load params

structs = {'clean_', 'cleanICA_', 'dws_filt_'}; % list of structs to correct-redefine

% the loops redefines one struct of all subjects at a time, before passing to the next struct
for a = 1:length(structs)
    clear subject data new_sampleinfo
    for i = 1:length(subject_list)
        clear subject data new_sampleinfo
        folderpath = fullfile(analysis_dir, subject_list{i}, '/eeg/'); % folderpath to explore data in "eeg" folder
        folderpath = fullfile(folderpath, '**');
        filelist   = dir(folderpath);
        name       = {filelist.name};
        name       = name(~strncmp(name, '.', 1));
        
        eeg_num1 = 'EEG.1'; eeg_num2 = 'EEG.2'; eeg_num3 = 'EEG.3';
        
        % if-statement to call the correct subject structure (number 1, 2, or 3)
        if contains(name(1), eeg_num1)
            structure = '_eeg_EEG_1'; substruct = '.EEG.1';
        elseif contains(name(1), eeg_num2)
            structure = '_eeg_EEG_2'; substruct = '.EEG.2';
        elseif contains(name(1), eeg_num3)
            structure = '_eeg_EEG_3'; substruct = '.EEG.3';
        else
            error('No EEG number identified, it could be EEG.4 or superior')
        end
        
        folderpath = fullfile(analysis_dir, subject_list{i}, '/processed/'); % folderpath re-named to call "processed" folder
        file = fullfile(folderpath, [structs{a}, subject_list{i}, structure, '.mat']);
        [~, filename ,~] = fileparts(file);
        load(file);
        
        % REPAIR
        try % check if data.sampleinfo is actually disorganized
            ft_fetch_data(data)
        catch ME
            if (strcmp(ME.identifier, 'FieldTrip:ft_fetch_data:line166'))
                fprintf('Problem identifier: %s, in %s \n', ME.identifier, filename);
                fprintf('Problem message: %s, in %s \n', ME.message, filename);
                repair = 1;
            end
        end
        
        if repair == 1
            nTrials = length(data.trial); % -----IT TAKES TOTAL NUMBER OF TRIALS TO REPAIR EVERYTHING-----
            trial_length = data.sampleinfo(1,2) / data.fsample; % in seconds
            dlim_trial = 1; % minimal sampleinfo of every trial
            uplim_trial = data.fsample * trial_length; % maximal duration of every trial
            
            % repair data.sampleinfo
            uplim = uplim_trial * nTrials;
            col1 = (dlim_trial:uplim_trial:uplim)';
            col2 = (uplim_trial:uplim_trial:uplim)';
            new_sampleinfo = cat(2, col1, col2);
            
            % now repair data.trial_markings_sampleinfo
            col_j1 = [];
            for k = 1:length(new_sampleinfo)
                col_j1{k} = cat(2,col1(k), col2(k));
            end
            col_j1 = col_j1'; % trial_markings_sampleinfo(A)
            
            dlim_trial = 1/data.fsample; % minimal trial_markings_sampleinfo(b1)
            uplim_trial = trial_length; % maximal trial_markings_sampleinfo(b1)
            uplim = uplim_trial * nTrials;
            col_j2 = (dlim_trial:uplim_trial:uplim)';
            col_j3 = (uplim_trial:uplim_trial:uplim)';
            
            col_j4 = [];
            for k = 1:length(new_sampleinfo)
                col_j4{k} = cat(2,col_j2(k), col_j3(k));
            end
            col_j4 = col_j4'; % trial_markings_sampleinfo(B)
            col_j5 = cat(2,col_j1,col_j4); % concatenate!
            new_trial_markings_sampleinfo = col_j5;
            
            % SAVE new data!
            data.sampleinfo = new_sampleinfo;
            data.trial_markings_sampleinfo = new_trial_markings_sampleinfo;
            subject.data_info.trial_markings_sampleinfo = new_trial_markings_sampleinfo;
            save(file, 'subject', 'data');
        end
    end
    fprintf('Sampleinfo repaired in %s \n', filename);
end