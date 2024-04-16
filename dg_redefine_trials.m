% dg_try_redefine_trials
clear all; close all;

fprintf('Select group to redefine trials: Dyslexia or Controls \n');
prompt = {'Answers (1 Dyslexia, 2 Controls)'};
dlgtitle = 'Redefine trials';
dims = [1 45];
answer = inputdlg(prompt, dlgtitle, dims);
[sel_group] = deal(answer{:});
sel_group = str2double(sel_group);

if sel_group == 1
    analysis_dir = 'patients dir';
    %     subject_list = {'pseudonyms'};
    subject_list = {'pseudonyms'}; % (EXCLUDED, conflicting redefinition; likely due to previous "glueing")
    group = 'Dyslexia_ST_PedTemp';
elseif sel_group == 2
    analysis_dir = 'controls dir';
    subject_list = {'pseudonyms'};
    group = 'Controls_ST_Dyslexia_PedTemp';
else
    error('Wrong selection');
end

cd(analysis_dir);
% addpath(genpath(fullfile(analysis_dir, '/scripts')));
analysis_params; %load params

%% define NEW TRIAL LENGTH and structs to correct, and run REDEFINE
structs = {'clean_', 'cleanICA_', 'dws_filt_', 'selected_trials_', 'subject_cache'}; % list of structs to correct-redefine
cfg.length = params.trial_length; % 2 seconds

% the loops redefines one struct of all subjects at a time, before passing to the next struct
for a = 1:length(structs)
    clear subject data new_subject new_data
    for i = 1:length(subject_list)
        clear subject data new_subject new_data
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
        
        % check struct to change and modify trials and data accordingly
        if a == 1  || a == 2 || a == 3  % 'clean_', 'clean_ICA', or % 'dws_filt_'
            folderpath = fullfile(analysis_dir, subject_list{i}, '/processed/'); % folderpath re-named to call "processed" folder
            file = fullfile(folderpath, [structs{a}, subject_list{i}, structure, '.mat']);
            [~, filename ,~] = fileparts(file);
            load(file);
            
            if isfield(data, 'bad_channels') % if bad channels were marked, then keep them!
                bad_chs = data.bad_channels;
                [ new_subject, new_data ] = nmri_redefine_trials( cfg, subject, data ); % redefine trials and related data
                subject = new_subject; data = new_data;
                data.bad_channels = bad_chs; % save bad channels again!
                save(file, 'data', 'subject');
                fprintf('Trials length and corresponding data redefined in: %s \n', filename);
            elseif ~isfield(data, 'bad_channels')
                [ new_subject, new_data ] = nmri_redefine_trials( cfg, subject, data ); % redefine trials and related data
                subject = new_subject; data = new_data;
                save(file, 'data', 'subject');
            end
        elseif a == 4 % 'selected_trials_'
            folderpath = fullfile(analysis_dir, subject_list{i}, '/processed/'); % folderpath re-named to call "processed" folder
            file = fullfile(folderpath, [structs{a}, subject_list{i}, structure, '.mat']);
            [~, filename ,~] = fileparts(file);
            if exist(file)
                load(file);
                % cleaned_s = load(fullfile(folderpath, ['cleanICA_', subject_list{i}, structure, '.mat'])); % open respective cleanICA struct
                % to check for markings and take the CORRECT TRIALS!
                trial_max = 10; % take the previous trial length as parameter to modify the data
                trial_min = cfg.length; % take the new trial length as parameter
                
                % change evt_goodTrials
                change1 = [];
                for j = 1:size(subject.evt_goodTrials, 2)
                    num = subject.evt_goodTrials(j);
                    dlim = (num*(trial_max/trial_min)) - ((trial_max/trial_min)-1);
                    change1{j} = dlim:1:(num*(trial_max/trial_min));
                end
                change1 = cat(2,change1{:}); % reshapes in ascendent order and 1-D array
                
                %             % but look for trials with bad technical markings and/or stimulation (already redefiend with nmri_redefine_trials)
                %             gd1 = []; gd2 = [];
                %             for tm = 1:size(cleaned_s.data.trial_markings, 1)
                %                 if cleaned_s.data.trial_markings{tm,2} == 0
                %                     gd1(tm) = tm;
                %                 end
                %
                %                 if strcmp(cleaned_s.data.trial_markings{tm,4},'eyes_open') || strcmp(cleaned_s.data.trial_markings{tm,4},'PS') || strcmp(cleaned_s.data.trial_markings{tm,4},'HV') || strcmp(cleaned_s.data.trial_markings{tm,4},'postHV')
                %                     gd2(tm) = tm;
                %                 end
                %             end
                %             gd = union(gd1,gd2);
                %
                %             if any(gd) % removes if there's a "zero" trial
                %                 gd(gd == 0) = [];
                %             end
                %             change1 = setdiff(change1, gd); % remove those trials that had more recently "bad" technical markings and/or stimulation (new goodTrials)
                %
                % change evt_badTrials
                change2 = [];
                for j = 1:size(subject.evt_badTrials, 2)
                    num = subject.evt_badTrials(j);
                    dlim = (num*(trial_max/trial_min)) - ((trial_max/trial_min)-1);
                    change2{j} = dlim:1:(num*(trial_max/trial_min));
                end
                change2 = cat(2,change2{:}); % reshapes in ascendent order and 1-D array
                
                % change SelectedTrials
                change3 = [];
                for j = 1:size(subject.SelectedTrials, 2)
                    num = subject.SelectedTrials(j);
                    dlim = (num*(trial_max/trial_min)) - ((trial_max/trial_min)-1);
                    change3{j} = dlim:1:(num*(trial_max/trial_min));
                end
                change3 = cat(2,change3{:}); % reshapes in ascendent order and 1-D array
                
                % save new selected_trials struct, with all changed/corrected fields
                subject.params = params; % replace params in subject for group params parameters
                subject.evt_goodTrials = change1;
                subject.evt_badTrials  = change2;
                subject.SelectedTrials = change3;
                save(file, 'subject');
                fprintf('Trials length and corresponding data redefined in: %s \n', filename);
            end
        elseif a == 5 % conf/subject_cache (since it has also selected-trials fields)
            structure = structure(2:end); % remove a '_' in the folder/struct name
            folderpath = fullfile(analysis_dir, subject_list{i}, '/conf/', structure); % folderpath re-named to call "processed" folder
            file = fullfile(folderpath, [structs{a} '.mat']);
            [~, filename ,~] = fileparts(file);
            if exist(file)
                load(file);
                % to check for markings and take the CORRECT TRIALS!
                trial_max = 10; % take the previous trial length as parameter to modify the data
                trial_min = cfg.length; % take the new trial length as parameter
                
                % change SelectedTrials (if field exists)
                if isfield(subject, 'SelectedTrials')
                    change3 = [];
                    for j = 1:size(subject.SelectedTrials, 2)
                        num = subject.SelectedTrials(j);
                        dlim = (num*(trial_max/trial_min)) - ((trial_max/trial_min)-1);
                        change3{j} = dlim:1:(num*(trial_max/trial_min));
                    end
                    change3 = cat(2,change3{:}); % reshapes in ascendent order and 1-D array
                    
                    % save new selected_trials struct, with all changed/corrected fields
                    subject.params = params; % replace params in subject for group params parameters
                    subject.SelectedTrials = change3;
                    save(file, 'subject');
                end
            end
        end
        fprintf('Trials length and corresponding data redefined in: %s \n', filename);
    end
end

% remove path loaded
rmpath(genpath(fullfile(analysis_dir, '/scripts')));
fprintf('Finished with trials redefinition in %s folder \n', group);