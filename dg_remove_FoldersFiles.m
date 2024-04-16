clear all; close all;

analysis_dir = 'group dir';
subject_list = {'pseudonyms'};
group = 'group name';

% PREPARE THINGS TO REMOVE
fprintf('PLEASE INDICATE WHAT TO REMOVE FROM SUBJECTS´ FOLDERS \n');
prompt = {'Remove conf?', 'Remove logs?', 'Remove mri?', 'Remove processed?', 'Remove QC?', 'Remove stats?'};
dlgtitle = 'What to remove? (1 Yes, 2 No)';
dims = [1 45];
answer = inputdlg(prompt, dlgtitle, dims);
[rem_conf, rem_logs, rem_mri, rem_processed, rem_QC, rem_stats] = deal(answer{:});

rem_conf    = str2num(rem_conf);
rem_logs    = str2num(rem_logs);
rem_mri     = str2num(rem_mri);
rem_processed = str2num(rem_processed);
rem_QC      = str2num(rem_QC);
rem_stats   = str2num(rem_stats);

if rem_processed == 2
    fprintf('Please, indicate if you want to still remove something SPECIFIC from processed folders \n');
    warning('Removal of dws_filt, clean, cleanICA, and ICA structs is not possible yet!');
    warning('If removal of processed folder was chosen before, dws_filt, clean, cleanICA, and ICA structs will be removed');
    prompt = {'Remove electrodes_aligned_canon?', 'Remove hdm_lead/only_canon?', 'Remove lead_only_canon?'};
    dlgtitle = 'What to remove from processed folder? (1 Yes, 2 No)';
    dims = [1 60];
    answer = inputdlg(prompt, dlgtitle, dims);
    [rem_elec, rem_hdm, rem_lead] = deal(answer{:});
    
    rem_elec = str2num(rem_elec);
    rem_hdm  = str2num(rem_hdm);
    rem_lead = str2num(rem_lead);
elseif rem_processed == 1
    rem_elec = 2;
    rem_hdm  = 2;
    rem_lead = 2;
    fprintf('ALL Processed folders will be removed! ... \n');
end

%% REMOVAL
for i = 1:length(subject_list)
    clearvars status folderpath filelist name; % clear message; clear messageid;
    
    % commands to get the correct numerated EEG file
    folderpath = fullfile(analysis_dir, subject_list{i}, '/eeg/'); % folderpath to explore data in "eeg" folder
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
    
    % REMOVE CORRESPONDING CONF
    if rem_conf == 1
        conf = fullfile(analysis_dir, subject_list{i}, '/conf/');
        [status, ~, ~] = rmdir(conf, 's');
        if status == 1
            fprintf('conf-folder from %s removed successfully \n', subject_list{i});
        elseif status == 0 || rem_conf == 2
            fprintf('conf-folder from %s not removed \n', subject_list{i});
        end  
    else
    end
    
    % REMOVE CORRESPONDING LOGS
    if rem_logs == 1
        logs = fullfile(analysis_dir, subject_list{i}, '/logs/');
        [status, ~, ~] = rmdir(logs, 's');
        if status == 1
            fprintf('logs-folder from %s removed successfully \n', subject_list{i});
        elseif status == 0
            fprintf('logs-folder from %s not removed \n', subject_list{i});
        end      
    elseif rem_logs == 2 && rem_elec == 1 && rem_hdm == 1 && rem_lead == 1 % if logs should not be removed but headmodels, then remove logs related to headmodel works
        logs = fullfile(analysis_dir, subject_list{i}, '/logs/', structure(2:end));
        hdm_log = fullfile(logs, 'hdmSUMA_nmriSUMA150_fs6_openmeeg.mat');
        proc_hdm_log = fullfile(logs, 'processing_nmriSUMA150_fs6_openmeeg.mat');
        
        if exist(hdm_log)
            delete(hdm_log)
        elseif ~exist(hdm_log)
            warning(['hdm_log was not found for ', subject_list{i}]);
        end
        
        if exist(proc_hdm_log)
            delete(proc_hdm_log)
        elseif ~exist(proc_hdm_log)
            warning(['proc_hdm_log was not found for ', subject_list{i}]);
        end
    end
    
    % REMOVE MRI FOLDER
    if rem_mri == 1
        mri = char(fullfile(analysis_dir, subject_list{i}, 'mri'));
        [status, ~, ~] = rmdir(mri, 's');
        
        if status == 1
            fprintf('mri-folder from %s removed successfully \n', subject_list{i});
        elseif status == 0 || rem_mri == 2
            fprintf('mri-folder from %s not removed \n', subject_list{i});
        end
    else
    end
    
    % REMOVE PROCESSED
    folderpath = fullfile(analysis_dir, subject_list{i}, '/processed/'); % folderpath to explore data in "eeg" folder
    if rem_processed==1 && rem_elec==2 && rem_hdm==2 && rem_lead==2
        proc = fullfile(analysis_dir, subject_list{i}, '/processed/');
        [status, ~, ~] = rmdir(proc, 's');
        if status == 1
            fprintf('proc-folder from %s removed successfully \n', subject_list{i});
        elseif status == 0 || rem_mri == 2
            fprintf('proc-folder from %s not removed \n', subject_list{i});
        end
    elseif rem_processed==2 && rem_elec==1
        % REMOVE ELECTRODES_ALIGNED_CANON
        aligned = fullfile(folderpath, ['electrodes_aligned_canon_' 'nmriSUMA150_fs6_' subject_list{i} structure '.mat']);
        if exist(aligned)
            delete(aligned)
        elseif ~exist(aligned)
            warning(['electrodes_aligned_canon file was not found for ', subject_list{i}]);
        end
    elseif rem_processed==2 && rem_hdm==1
        % REMOVE HDM_ONLY_CANON FILE
        hdm = fullfile(folderpath, ['hdm_only_canon_', subject_list{i}, structure, '_nmriSUMA150_fs6_openmeeg.mat']);
        if exist(hdm)
            delete(hdm)
        elseif ~exist(hdm)
             warning(['hdm_only_canon_ file was not found for ', subject_list{i}]);
        end
    elseif rem_processed==2 && rem_lead==1
        hdm1 = fullfile(folderpath, ['hdm_lead_canon_', subject_list{i}, structure, '_nmriSUMA150_fs6_openmeeg.mat']);
        % REMOVE LEAD_ONLY_CANON FILE
        if exist(hdm1)
            delete(hdm1)
        elseif ~exist(hdm1)
            warning(['lead_only_canon file was not found for ', subject_list{i}]);
        elseif ~exist(hdm) && ~exist(hdm1)
            warning(['hdm_only_canon file was not found for ', subject_list{i}]);
        end
    else
    end
    
    % REMOVE QC FOLDER
    if rem_QC == 1
        qc = char(fullfile(analysis_dir, subject_list{i}, 'QC'));
        [status, ~, ~] = rmdir(qc, 's');
        if status == 1
            fprintf('QC-folder from %s removed successfully \n', subject_list{i});
        elseif status == 0
            fprintf('QC-folder from %s not removed successfully \n', subject_list{i});
        end
    else
    end
    
    % REMOVE STATS FOLDER
    if rem_stats == 1
        folder_1 = char(fullfile(analysis_dir, subject_list{i}, 'stats'));
        [status, ~, ~] = rmdir(folder_1, 's');
        if status == 1
            fprintf('stats-folder from %s removed successfully \n', subject_list{i});
        elseif status == 0
            fprintf('stats-folder from %s not removed successfully \n', subject_list{i});
        end
    else
    end
end

fprintf('Done with removal of selected folders in %s group \n', group);