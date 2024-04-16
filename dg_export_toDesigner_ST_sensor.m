% dg_export_ST_sensor
analysis_folder = 'main dir';

%% Controls cell
pwd = 'pwd';
addpath(genpath('script subfolder'));
[all_subjects] = nmri_all_subjects(pwd);
% filters.id = {'pseudonyms'};
filters.id = {'pseudonyms'};
[controls_ST_sensor] = nmri_filter_subjects(all_subjects, filters);
save('controls_ST_sensor.mat', 'controls_ST_sensor');
rmpath(genpath(pwd));

%% Dyslexia without IEDs cell
clearvars;
pwd = 'pwd';
addpath(genpath('script subfolder'));[all_subjects] = nmri_all_subjects(pwd);
filters.id = {'pseudonyms'};

[dys_noIEDs_ST_sensor] = nmri_filter_subjects(all_subjects, filters);
save('dys_noIEDs_ST_sensor.mat', 'dys_noIEDs_ST_sensor');
rmpath(genpath(pwd));

%% EXPORT FOR SENSOR LEVEL (USING SHORT TRIALS!) - WHOLE GROUPS: 70 Dys vs 50 Controls
clearvars;
pwd = 'pwd';
controls_str = fullfile(pwd, 'controls_ST_sensor.mat');
load(controls_str);
dys_noIEDs_str = fullfile(pwd, 'dys_noIEDs_ST_sensor.mat');
load(dys_noIEDs_str);
all_subjects = [controls_ST_sensor; dys_noIEDs_ST_sensor]; % Concatenate!

% NOW EXPORT (all_subjects cell and export)
addpath(genpath('script subfolder'));
analysis_params;

opt = [];
opt.metrics = {'coh_img','power','wpli_debiased'};
opt.scale = {'none'}; % no global scaling
opt.global = {'abs'}; % take absolute values
opt.dim_reduction = {'mean'}; % use the mean to reduce to one value per vertex
opt.save_grp_mean = true; % save a mean map
opt.interpolate = true;
opt.hdm_class = 'interp_sensor';
opt.output = fullfile(pwd,'export1_Dys_Controls_ST_sensor');
nmri_export_metrics_sensor(all_subjects, opt, params);

%% SUB EXPORT FOR SENSOR LEVEL (USING SHORT TRIALS!): DYSLEXIA AND CONTROLS, ONLY THOSE IN 1-2 & 3-4 SCHOOL GRADES
clearvars;
pwd = 'pwd';
addpath(genpath('script subfolder'));[all_subjects] = nmri_all_subjects(pwd);
filters.id = {'pseudonyms'};
[controls_1to4_ST_sensor] = nmri_filter_subjects(all_subjects, filters);
save('controls_1to4_ST_sensor.mat', 'controls_1to4_ST_sensor');
rmpath(genpath(pwd));

clearvars;
pwd = 'pwd';
addpath(genpath('script subfolder'));
[all_subjects] = nmri_all_subjects(pwd);
filters.id = {'pseudonyms'};
[dys_1to4_ST_sensor] = nmri_filter_subjects(all_subjects, filters);
save('dys_1to4_ST_sensor.mat', 'dys_1to4_ST_sensor');
rmpath(genpath(pwd));

clearvars;
pwd = 'pwd';
controls_str = fullfile(pwd, 'controls_1to4_ST_sensor.mat');
load(controls_str);
dys_noIEDs_str = fullfile(pwd, 'dys_1to4_ST_sensor.mat');
load(dys_noIEDs_str);
all_subjects = [controls_1to4_ST_sensor; dys_1to4_ST_sensor]; % Concatenate!

% NOW EXPORT (all_subjects cell and export)
addpath(genpath('script subfolder'));
analysis_params;

opt = [];
opt.metrics = {'coh_img','power','wpli_debiased'};
opt.scale = {'none'}; % no global scaling
opt.global = {'abs'}; % take absolute values
opt.dim_reduction = {'mean'}; % use the mean to reduce to one value per vertex
opt.save_grp_mean = true; % save a mean map
opt.interpolate = true;
opt.hdm_class = 'interp_sensor';
opt.output = fullfile(pwd,'export2_Dys_Controls_1to4_ST_sensor');
nmri_export_metrics_sensor(all_subjects, opt, params);

%% SUB EXPORT FOR SENSOR LEVEL (USING SHORT TRIALS!): DYSLEXIA AND CONTROLS, ONLY THOSE IN 5-6 & 7-8 SCHOOL GRADES SCHOOL GRADES
clearvars;
pwd = '/pwd';
addpath(genpath('script subfolder'));
[all_subjects] = nmri_all_subjects(pwd);
filters.id = {'pseudonyms'};
[controls_5to8_ST_sensor] = nmri_filter_subjects(all_subjects, filters);
save('controls_5to8_ST_sensor.mat', 'controls_5to8_ST_sensor');
rmpath(genpath(pwd));

clearvars;
pwd = 'pwd';
addpath(genpath('script subfolder'));
[all_subjects] = nmri_all_subjects(pwd);
filters.id = {'pseudonums'};
[dys_5to8_ST_sensor] = nmri_filter_subjects(all_subjects, filters);
save('dys_5to8_ST_sensor.mat', 'dys_5to8_ST_sensor');
rmpath(genpath(pwd));

clearvars;
pwd = 'pwd';
controls_str = fullfile(pwd, 'controls_5to8_ST_sensor.mat');
load(controls_str);
dys_noIEDs_str = fullfile(pwd, 'dys_5to8_ST_sensor.mat');
load(dys_noIEDs_str);
all_subjects = [controls_5to8_ST_sensor; dys_5to8_ST_sensor]; % Concatenate
% NOW EXPORT (all_subjects cell and export)
addpath(genpath('script subfolder'));
analysis_params;

opt = [];
opt.metrics = {'coh_img','power','wpli_debiased'};
opt.scale = {'none'}; % no global scaling
opt.global = {'abs'}; % take absolute values
opt.dim_reduction = {'mean'}; % use the mean to reduce to one value per vertex
opt.save_grp_mean = true; % save a mean map
opt.interpolate = true;
opt.hdm_class = 'interp_sensor';
opt.output = fullfile(pwd,'export3_Dys_Controls_5to8_ST_sensor');
nmri_export_metrics_sensor(all_subjects, opt, params);

%% Check that everything is "in order"
for i = 1:120
    fprintf('Subject: %s \n', all_subjects{i,1}.loaded_from);
end

for i = 1:120
    fprintf('Subject: %d \n', length(all_subjects{i,1}.SelectedTrials));
end

for i = 1:96
    fprintf('Subject: %s \n', all_subjects{i,1}.loaded_from);
end

for i = 1:96
    fprintf('Subject: %d \n', length(all_subjects{i,1}.SelectedTrials));
end

for i = 1:24
    fprintf('Subject: %s \n', all_subjects{i,1}.loaded_from);
end

for i = 1:24
    fprintf('Subject: %d \n', length(all_subjects{i,1}.SelectedTrials));
end