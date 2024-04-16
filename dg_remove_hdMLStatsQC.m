clearvars;
analysis_dir = 'group dir';

subject_list = {'pseudonyms'};
group = 'group name';

for i = 1:length(subject_list)
    clear status; % clear message; clear messageid;
    
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
    
    folderpath = fullfile(analysis_dir, subject_list{i}, '/processed/'); % folderpath re-named to call "processed" folder
    
%     % REMOVE HDM_ONLY_CANON FILE
%     hdm = fullfile(folderpath, ['hdm_only_canon_', subject_list{i}, structure, '_nmriSUMA150_fs6_openmeeg.mat']);
%     hdm1 = fullfile(folderpath, ['hdm_lead_canon_', subject_list{i}, structure, '_nmriSUMA150_fs6_openmeeg.mat']);
%     
%     if exist(hdm)
%         delete(hdm)
%     elseif exist(hdm1)
%         delete(hdm1)
%     elseif ~exist(hdm) && ~exist(hdm1)
%         warning(['hdm_only_canon file was not found for ', subject_list{i}]);
%     end
%     
%     % REMOVE LEAD_ONLY_CANON FILE
%     lead = fullfile(folderpath, ['lead_only_canon_', subject_list{i}, structure, '_nmriSUMA150_fs6_openmeeg.mat']);
%     
%     if exist(lead)
%         delete(lead)
%     elseif ~exist(lead)
%         warning(['lead_only_canon file was not found for ', subject_list{i}]);
%     end
%     
%     % REMOVE ELECTRODES_ALIGNED_CANON
%     aligned = fullfile(folderpath, ['electrodes_aligned_canon_' 'nmriSUMA150_fs6_' subject_list{i} structure '.mat']);
%     
%     if exist(aligned)
%         delete(aligned)
%     elseif ~exist(aligned)
%         warning(['electrodes_aligned_canon file was not found for ', subject_list{i}]);
%     end
%     
    % REMOVE STATS FOLDER
    folder_1 = char(fullfile(analysis_dir, subject_list{i}, 'stats'));
    
    if exist(folder_1)
        [status, ~, ~] = rmdir(folder_1, 's');
    elseif ~exist(folder_1)
        fprintf('stats folder not found');
        status = 0;
    end
    if status == 1
        fprintf('stats-folder from %s removed successfully \n', subject_list{i});
    elseif status == 0
        fprintf('stats-folder from %s not removed \n', subject_list{i});
    end
%     
%     % REMOVE QC FOLDER
%     folder_2 = char(fullfile(analysis_dir, subject_list{i}, 'QC'));
%     [status, ~, ~] = rmdir(folder_2, 's');
%     
%     if status == 1
%         fprintf('QC-folder from %s removed successfully \n', subject_list{i});
%     elseif status == 0
%         fprintf('QC-folder from %s not removed successfully \n', subject_list{i});
%     end
%     
%     % REMOVE MRI FOLDER
%     mri = char(fullfile(analysis_dir, subject_list{i}, 'mri'));
%     [status, ~, ~] = rmdir(mri, 's');
%     
%     if status == 1
%         fprintf('mri-folder from %s removed successfully \n', subject_list{i});
%     elseif status == 0
%         fprintf('mri-folder from %s not removed successfully \n', subject_list{i});
%     end
%     
%     % REMOVE CORRESPONDING LOGS
%     folderpath = fullfile(analysis_dir, subject_list{i}, '/logs/', structure(2:end)); % folderpath re-named to call "logs" folder
%     hdm_log = fullfile(folderpath, ['hdmSUMA_nmriSUMA150_fs6_openmeeg.mat']);
%     proc_hdm_log = fullfile(folderpath, ['processing_nmriSUMA150_fs6_openmeeg.mat']);
%         
%     if exist(hdm_log)
%         delete(hdm_log)
%     elseif ~exist(hdm_log)
%         warning(['hdm_log was not found for ', subject_list{i}]);
%     end
%     
%     if exist(proc_hdm_log)
%         delete(proc_hdm_log)
%     elseif ~exist(proc_hdm_log)
%         warning(['proc_hdm_log was not found for ', subject_list{i}]);
%     end
end

fprintf('Done with removal of HdML files, stats and QC folders in %s group \n', group);