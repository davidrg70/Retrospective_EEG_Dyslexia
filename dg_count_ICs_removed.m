% dg_count_ICs_removed - counts number of components removed during ICA (preprocessing)
clear all;

analysis_dir1 = 'patients dir';
analysis_dir2 = 'controls dir';

% patients Dys_no_IEDs
patients_list = {'pseudonyms'};

% controls
controls_list = {'pseudonyms'};

%% count for patients
for i = 1:length(patients_list)
    clear selected; clear comp;
    % commands to get the correct numerated EEG file
    folderpath = fullfile(analysis_dir1, patients_list{i}, '/eeg/'); % folderpath to explore data in "eeg" folder
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
    
    folderpath = fullfile(analysis_dir1, patients_list{i}, '/processed/'); % folderpath to upload ICA_comp struct
    load([folderpath 'ICA_comp_' patients_list{i} structure '.mat']); % loads subject struct
    unselected_patients(i) = nnz(~selected); % counts number of removed ICs in/after ICA
end

%% count for controls
for i = 1:length(controls_list)
    clear selected; clear comp;
    % commands to get the correct numerated EEG file
    folderpath = fullfile(analysis_dir2, controls_list{i}, '/eeg/'); % folderpath to explore data in "eeg" folder
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
    
    folderpath = fullfile(analysis_dir2, controls_list{i}, '/processed/'); % folderpath to upload ICA_comp struct
    load([folderpath 'ICA_comp_' controls_list{i} structure '.mat']); % loads subject struct
    unselected_controls(i) = nnz(~selected); % counts number of removed ICs in/after ICA
end

%% averages
ics_removed_patients = mean(unselected_patients);
ics_removed_controls = mean(unselected_controls);
ics_removed_all = (ics_removed_patients + ics_removed_controls)/2;

ics_std_patients = std(unselected_patients);
ics_std_controls = std(unselected_controls);
ics_std_all = (ics_std_patients + ics_std_controls)/2;