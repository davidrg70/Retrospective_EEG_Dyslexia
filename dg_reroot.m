% Reroots all subjects folders-structures with a new folder given

%% Reroots all structures in "preprocessed"
clear all; close all;
new_analysis_dir = 'group dir';

subject_list = {'pseudonyms'};

struct_to_correct = {'clean_', 'cleanICA_', 'dws_filt_', 'electrodes_aligned_canon_nmriSUMA150_fs6_', ...
    'hdm_only_canon_', 'lead_only_canon_', 'selected_trials_'}; % 'ICA_comp_' was omitted, since it has just comp and selected fields (no subject)

% indicate if to delete selected_trials structure, mri folder and related processing structs
fprintf('Indicate if Selected_trials structure (in processed folder) must be deleted \n');
prompt = {'Remove selected_trials?', 'Remove mri folder and related processing structs?'};
dlgtitle = ['Remove selected_trials, mri folder and canonical processing structs?', ' -> Answers (1 Yes, 2 No)'];
dims = [1 100];
answer = inputdlg(prompt, dlgtitle, dims);
[del_sel_trials, del_mri_folder] = deal(answer{:});
del_sel_trials = str2double(del_sel_trials);
del_mri_folder = str2double(del_mri_folder);

% the loops takes every struct and re-roots all subjects, before passing to the next struct
for a = 1:length(struct_to_correct)
    clear subject data
    for i = 1:length(subject_list)
        clear subject data
        folderpath = fullfile(new_analysis_dir, subject_list{i}, '/eeg/'); % folderpath to explore data in "eeg" folder
        folderpath = fullfile(folderpath, '**');
        filelist   = dir(folderpath);
        name       = {filelist.name};
        name       = name(~strncmp(name, '.', 1));
        
        eeg_num1 = 'EEG.1'; eeg_num2 = 'EEG.2'; eeg_num3 = 'EEG.3';
        
        % if-statement to call the correct subject structure (number 1, 2, or 3)
        if contains(name(1), eeg_num1)
            structure = '_eeg_EEG_1';
            substruct = '.EEG.1';
        elseif contains(name(1), eeg_num2)
            structure = '_eeg_EEG_2';
            substruct = '.EEG.2';
        elseif contains(name(1), eeg_num3)
            structure = '_eeg_EEG_3';
            substruct = '.EEG.3';
        else
            error('No EEG number identified, it could be EEG.4 or superior')
        end
        
        folderpath = fullfile(new_analysis_dir, subject_list{i}, '/processed/'); % folderpath re-named to call "processed" folder
        
        % check structure to call, either 1 to 7, and because the 'hdm_only_canon_nmriSUMA150_' and 'lead_only_canon' structures must be called with additional text info
        if a == 1 || a == 2 || a == 3 % 'clean_', 'cleanICA_', 'dws_filt_'
            % Load those to re-root!
            subject_structure = fullfile(folderpath, [struct_to_correct{a} subject_list{i} structure '.mat']);
            load(subject_structure);
            % Deal with stamps, according to previous selection of mri folders and related processing structs
            fields = {'hdmSUMA_nmriSUMA150_fs6_openmeeg','processing_nmriSUMA150_fs6_openmeeg','processing_sensor'};
            if isfield(subject, 'stamps')
                if isfield(subject.stamps, fields{1}) && del_mri_folder == 1
                    subject.stamps = rmfield(subject.stamps, fields{1});
                elseif isfield(subject.stamps, fields{2}) && ~exist(fullfile(new_analysis_dir, subject_list{i}, 'stats'))
                    subject.stamps = rmfield(subject.stamps, fields{2});
                elseif isfield(subject.stamps, fields{3}) && ~exist(fullfile(new_analysis_dir, subject_list{i}, 'stats'))
                    subject.stamps = rmfield(subject.stamps, fields{3});
                end
            else
            end
            % Re-root!
            [subject] = nmri_reroot_subject(subject, new_analysis_dir);
            % Save changes!
            save(fullfile(folderpath, [struct_to_correct{a} subject_list{i} structure '.mat']), 'subject', 'data');
        elseif a == 4 && del_mri_folder == 1 % delete 'electrodes_aligned_canon_', if indicated
            subject_structure = fullfile(folderpath, [struct_to_correct{a} subject_list{i} structure '.mat']);
            if exist(subject_structure)
                delete(subject_structure)
                fprintf('Deleted %s struct of %s \n', struct_to_correct{a}, subject_list{i});
            else
                warning('Struct %s from %s does not exist', struct_to_correct{a}, subject_list{i});
            end
        elseif a == 5 && del_mri_folder == 1 % delete 'hdm_only_canon_', if indicated
            subject_structure = fullfile(folderpath, [struct_to_correct{a} subject_list{i} structure '_nmriSUMA150_fs6_openmeeg' '.mat']);
            if exist(subject_structure)
                delete(subject_structure)
                fprintf('Deleted %s struct of %s \n', struct_to_correct{a}, subject_list{i});
            else
                warning('Struct %s from %s does not exist', struct_to_correct{a}, subject_list{i});
            end
        elseif a == 6 && del_mri_folder == 1 % delete 'lead_only_canon', if indicated
            subject_structure = fullfile(folderpath, [struct_to_correct{a} subject_list{i} structure, '_nmriSUMA150_fs6_openmeeg', '.mat']);
            if exist(subject_structure)
                delete(subject_structure)
                fprintf('Deleted %s struct of %s \n', struct_to_correct{a}, subject_list{i});
            else
                warning('Struct %s from %s does not exist', struct_to_correct{a}, subject_list{i});
            end
        elseif a == 6 && del_mri_folder == 2 % if mri folders and 'lead_only_canon' shouldn't be deleted, then re-root them!
            subject_structure = fullfile(folderpath, [struct_to_correct{a} subject_list{i} structure, '_nmriSUMA150_fs6_openmeeg', '.mat']);
            if exist(subject_structure)
                load(subject_structure);
                [subject] = nmri_reroot_subject(subject, new_analysis_dir);
                save(subject_structure, 'bnd', 'canonical_subject', 'leadfield', 'subject');
            else
                warning('Struct %s from %s does not exist', struct_to_correct{a}, subject_list{i});
            end
        elseif a == 7 && del_sel_trials == 1 % delete 'selected_trials', if indicated
            subject_structure = fullfile(folderpath, [struct_to_correct{a} subject_list{i} structure '.mat']);
            if exist(subject_structure)
                delete(subject_structure)
                fprintf('Deleted %s struct of %s \n', struct_to_correct{a}, subject_list{i});
            else
                warning('Struct %s from %s does not exist', struct_to_correct{a}, subject_list{i});
            end
        elseif a == 7 && del_sel_trials == 2 % if selected_trials exist shouldn't be deleted, then re-root it!
            % Load selected_trials
            subject_structure = fullfile(folderpath, [struct_to_correct{a} subject_list{i} structure '.mat']);
            % DEAL WITH STAMPS IN SELECTED_TRIALS STRUCT
            if exist(subject_structure)
                load(subject_structure);
                % remove HdmL and Proc stamps! Depending on previous selection and if STATS folder was erased before
                fields = {'hdmSUMA_nmriSUMA150_fs6_openmeeg','processing_nmriSUMA150_fs6_openmeeg','processing_sensor'};
                if isfield(subject.stamps, fields{1}) && del_mri_folder == 1
                    subject.stamps = rmfield(subject.stamps, fields{1});
                elseif isfield(subject.stamps, fields{2}) && ~exist(fullfile(new_analysis_dir, subject_list{i}, 'stats'))
                    subject.stamps = rmfield(subject.stamps, fields{2});
                elseif isfield(subject.stamps, fields{3}) && ~exist(fullfile(new_analysis_dir, subject_list{i}, 'stats'))
                    subject.stamps = rmfield(subject.stamps, fields{3});
                end
                % Re-root!
                [subject] = nmri_reroot_subject(subject, new_analysis_dir);
                % Save changes!
                save(fullfile(folderpath, [struct_to_correct{a} subject_list{i} structure '.mat']), 'subject');
            elseif ~exist(subject_structure)
                warning('Struct %s from %s does not exist', struct_to_correct{a}, subject_list{i});
            end
        end
        fprintf('Rerooting for %s_%s successful \n', struct_to_correct{a}, subject_list{i});
    end
end
fprintf('All "processed" structures rerooted in <<%s>>\n', new_analysis_dir);

%% Reroots all structs in conf (subject_cache and subject_info) and deals with stamps in them
clearvars -except new_analysis_dir subject_list del_sel_trials del_mri_folder

for i = 1:length(subject_list)
    folderpath = fullfile(new_analysis_dir, subject_list{i}, '/eeg/'); % folderpath to explore data in "eeg" folder
    folderpath = fullfile(folderpath, '**');
    filelist   = dir(folderpath);
    name       = {filelist.name};
    name       = name(~strncmp(name, '.', 1));
    
    eeg_num1 = 'EEG.1'; eeg_num2 = 'EEG.2'; eeg_num3 = 'EEG.3';
    
    % if-statement to call the correct subject structure (number 1, 2, or 3)
    if contains(name(1), eeg_num1)
        structure = 'eeg_EEG_1';
        substruct = '.EEG.1';
    elseif contains(name(1), eeg_num2)
        structure = 'eeg_EEG_2';
        substruct = '.EEG.2';
    elseif contains(name(1), eeg_num3)
        structure = 'eeg_EEG_3';
        substruct = '.EEG.3';
    else
        error('No EEG number identified, it could be EEG.4 or superior')
    end
    
    folderpath = fullfile(new_analysis_dir, subject_list{i}, '/conf/', structure); % folderpath re-named to call "conf" folder
    conf_structs = {'subject_cache', 'subject_info'};
    
    for a = 1:length(conf_structs)
        clear subject timestamp source max_timestamp
        load(fullfile(folderpath, conf_structs{a}));
        if a == 1 % remove HdmL and Proc stamps in Subject_cache! Depending on previous selection and if STATS folder was erased before
            fields = {'hdmSUMA_nmriSUMA150_fs6_openmeeg','processing_nmriSUMA150_fs6_openmeeg','processing_sensor'};
            if isfield(subject.stamps, fields{1}) && del_mri_folder == 1
                subject.stamps = rmfield(subject.stamps, fields{1});
            end
            if isfield(subject.stamps, fields{2}) && ~exist(fullfile(new_analysis_dir, subject_list{i}, 'stats'))
                subject.stamps = rmfield(subject.stamps, fields{2});
            end
            if isfield(subject.stamps, fields{3}) && ~exist(fullfile(new_analysis_dir, subject_list{i}, 'stats'))
                subject.stamps = rmfield(subject.stamps, fields{3});
            end
        end
        % Re-root!
        [subject] = nmri_reroot_subject(subject, new_analysis_dir);
        % Save changes!
        if a == 1
            save(fullfile(folderpath, conf_structs{a}), 'subject', 'max_timestamp', 'source', 'timestamp');
        elseif a == 2
            save(fullfile(folderpath, conf_structs{a}), 'subject');
        end
        fprintf('Rerooting for CONF %s_%s successful \n', conf_structs{a}, subject_list{i});
    end
end
fprintf('"HdmL" and "Proc" stamps found in conf/subject_cache modified, in <<%s>>\n', new_analysis_dir);

%% Removes "mri" folder, if requested
clearvars -except new_analysis_dir subject_list del_sel_trials del_mri_folder

if del_mri_folder == 1
    % mri_struct = {'suma_surface_canon_nmriSUMA150_fs6'};
    for i = 1:length(subject_list)
        folderpath = fullfile(new_analysis_dir, subject_list{i}, '/eeg/'); % folderpath to explore data in "eeg" folder
        folderpath = fullfile(folderpath, '**');
        filelist   = dir(folderpath);
        name       = {filelist.name};
        name       = name(~strncmp(name, '.', 1));
        
        eeg_num1 = 'EEG.1'; eeg_num2 = 'EEG.2'; eeg_num3 = 'EEG.3';
        
        % if-statement to call the correct subject structure (number 1, 2, or 3)
        if contains(name(1), eeg_num1)
            structure = '_eeg_EEG_1';
            substruct = '.EEG.1';
        elseif contains(name(1), eeg_num2)
            structure = '_eeg_EEG_2';
            substruct = '.EEG.2';
        elseif contains(name(1), eeg_num3)
            structure = '_eeg_EEG_3';
            substruct = '.EEG.3';
        else
            error('No EEG number identified, it could be EEG.4 or superior')
        end
        
        mri_folder = fullfile(new_analysis_dir, subject_list{i}, '/mri/'); % folderpath re-named to call "mri" folder
        if exist(mri_folder)
            rmdir(mri_folder, 's');
            fprintf('mri-folder from %s removed successfully \n', subject_list{i});
        else
            warning('mri folder from %s does not exist', subject_list{i});
        end
    end
    fprintf('All "mri" structures removed in <<%s>>\n', new_analysis_dir);
end

%% Reroots all structures in "logs"
clearvars -except new_analysis_dir subject_list del_sel_trials del_mri_folder

log_to_correct = {'artifactrejection', 'artifactrejectionICA', 'artifactrejectionICA_estimation', ...
    'hdmSUMA_nmriSUMA150_fs6_openmeeg', 'preproc', 'processing_nmriSUMA150_fs6_openmeeg', 'processing_sensor', ...
    'readingevents', 'vigilance'};

% indicate if to delete Selected_trials structure
fprintf('Indicate if hdmSUMA, processing_nmriSUMA, or vigilance logs should be deleted \n');
prompt = {'Remove hdmSUMA log?', 'Remove processing_nmriSUMA log?', 'Remove processing_sensor log?', 'Remove vigilance log?'};
dlgtitle = 'Answers (1 Yes, 2 No)';
dims = [1 45];
answer = inputdlg(prompt, dlgtitle, dims);
[del_hdm, del_proc, del_procS, del_vig] = deal(answer{:});
del_hdm = str2double(del_hdm);
del_proc = str2double(del_proc);
del_procS = str2double(del_procS);
del_vig = str2double(del_vig);

% the loops take first every subject and then re-root all its structures, before to pass to the next subject
for a = 1:length(log_to_correct)
    clear subject params
    for i = 1:length(subject_list)
        clear subject params
        folderpath = fullfile(new_analysis_dir, subject_list{i}, '/eeg/'); % folderpath to explore data in "eeg" folder
        folderpath = fullfile(folderpath, '**');
        filelist   = dir(folderpath);
        name       = {filelist.name};
        name       = name(~strncmp(name, '.', 1));
        
        eeg_num1 = 'EEG.1'; eeg_num2 = 'EEG.2'; eeg_num3 = 'EEG.3';
        
        % if-statement to call the correct subject structure (number 1, 2, or 3)
        if contains(name(1), eeg_num1)
            structure = '_eeg_EEG_1';
            substruct = '.EEG.1';
        elseif contains(name(1), eeg_num2)
            structure = '_eeg_EEG_2';
            substruct = '.EEG.2';
        elseif contains(name(1), eeg_num3)
            structure = '_eeg_EEG_3';
            substruct = '.EEG.3';
        else
            error('No EEG number identified, it could be EEG.4 or superior')
        end
        
        folderpath = fullfile(new_analysis_dir, subject_list{i}, '/logs/', structure(2:end)); % folderpath re-named
        log_structure = fullfile(folderpath, [log_to_correct{a} '.mat']);
        
        % Check if structure was created before
        if exist(log_structure)
            if a == 4 && del_hdm == 1 % delete 'hdmSUMA' log
                delete(log_structure);
                fprintf('Deleted %s log of %s \n', log_to_correct{a}, subject_list{i});
            elseif a == 6 && del_proc == 1
                delete(log_structure) % delete 'processing_nmriSUMA150' log
                fprintf('Deleted %s log of %s \n', log_to_correct{a}, subject_list{i});
            elseif a == 7 && del_procS == 1
                delete(log_structure) % delete 'processing_sensor' log
                fprintf('Deleted %s log of %s \n', log_to_correct{a}, subject_list{i});
            elseif a == 9 && del_vig == 1
                delete(log_structure) % delete 'vigilance' log
                fprintf('Deleted %s log of %s \n', log_to_correct{a}, subject_list{i});
            else
                % Load structure
                load(log_structure);
                % Deal with stamps, according to previous selection of mri folders and related processing structs
                fields = {'hdmSUMA_nmriSUMA150_fs6_openmeeg','processing_nmriSUMA150_fs6_openmeeg','processing_sensor'};
                if isfield(subject, 'stamps')
                    if isfield(subject.stamps, fields{1}) && del_mri_folder == 1
                        subject.stamps = rmfield(subject.stamps, fields{1});
                    end
                    if isfield(subject.stamps, fields{2}) && ~exist(fullfile(new_analysis_dir, subject_list{i}, 'stats'))
                        subject.stamps = rmfield(subject.stamps, fields{2});
                    end
                    if isfield(subject.stamps, fields{3}) && ~exist(fullfile(new_analysis_dir, subject_list{i}, 'stats'))
                        subject.stamps = rmfield(subject.stamps, fields{3});
                    end
                end
                % Re-root!
                [subject] = nmri_reroot_subject(subject, new_analysis_dir);
                if exist('params')  % Check if a params field is present in the loaded structure
                    % then save subject and params fields
                    save(fullfile(folderpath, [log_to_correct{a} '.mat']), 'subject', 'params');
                    fprintf('Rerooting for %s_%s finished \n', log_to_correct{a}, subject_list{i});
                else
                    % then save only subject field
                    save(fullfile(folderpath, [log_to_correct{a} '.mat']), 'subject');
                    fprintf('Rerooting for %s_%s finished \n', log_to_correct{a}, subject_list{i});
                end
            end
        elseif ~exist(log_structure)
            fprintf('Log structrure not existent for %s_%s \n', log_to_correct{a}, subject_list{i});
        end
    end
end
fprintf('All "logs" rerooted in <<%s>>\n', new_analysis_dir);