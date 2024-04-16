% dg_export_ST_source
analysis_folder = 'main dir';

%% Controls cell
pwd = 'pwd';
addpath(genpath('script subfolder'));
[all_subjects] = nmri_all_subjects(pwd);
% filters.id = {'pseudonyms'};
filters.id = {'pseudonyms'};
[controls_ST] = nmri_filter_subjects(all_subjects, filters);
save('controls_ST.mat', 'controls_ST');
rmpath(genpath(pwd));

%% Dyslexia without IEDs cell
clearvars;
pwd = 'pwd';
addpath(genpath('script subfolder'));
[all_subjects] = nmri_all_subjects(pwd);
filters.id = {'pseudonyms'};

[dys_noIEDs_ST] = nmri_filter_subjects(all_subjects, filters);
save('dys_noIEDs_ST.mat', 'dys_noIEDs_ST');
rmpath(genpath(pwd));

%% EXPORT FOR SOURCE LEVEL (USING NEW PEDIATRIC TEMPLATE AND SHORT TRIALS!) - WHOLE GROUPS: 70 Dys vs 50 Controls
clearvars;
pwd = 'pwd';
controls_str = fullfile(pwd, 'controls_ST.mat');
load(controls_str);
dys_noIEDs_str = fullfile(pwd, 'dys_noIEDs_ST.mat');
load(dys_noIEDs_str);
all_subjects = [controls_ST; dys_noIEDs_ST]; % Concatenate!

% NOW EXPORT (all_subjects cell and export)
addpath(genpath('script subfolder'));
analysis_params;

opt = [];
opt.metrics = {'coh_img','power','wpli_debiased'};
opt.scale = {'none'}; % no global scaling
opt.global = {'abs'}; % take absolute values
opt.dim_reduction = {'mean'}; % use the mean to reduce to one value per vertex
opt.save_grp_mean = true; % save a mean map
opt.hdm_class = 'ANTS9Years3T_fs7_openmeeg';
opt.output = fullfile(pwd,'export2_Dys_Controls_ST_source');
nmri_export_metrics(all_subjects, opt, params);

%% SUB EXPORT: DYSLEXIA AND CONTROLS, ONLY THOSE IN 1-2 & 3-4 SCHOOL GRADES (USING NEW PEDIATRIC TEMPLATE AND SHORT TRIALS!)
% (Controls 6.98-10.80 y.o. &  Dyslexia 7.02-10.95 y.o.)
clearvars;
pwd = 'pwd';
addpath(genpath('script subfolder'));
[all_subjects] = nmri_all_subjects(pwd);
filters.id = {'pseudonyms'};
[controls_1to4_ST] = nmri_filter_subjects(all_subjects, filters);
save('controls_1to4_ST.mat', 'controls_1to4_ST');
rmpath(genpath(pwd));

clearvars;
pwd = 'pwd';
addpath(genpath('script subfolder'));
[all_subjects] = nmri_all_subjects(pwd);
filters.id = {'pseudonums'};
[dys_1to4_ST] = nmri_filter_subjects(all_subjects, filters);
save('dys_1to4_ST.mat', 'dys_1to4_ST');
rmpath(genpath(pwd));

clearvars;
pwd = 'pwd';
controls_str = fullfile(pwd, 'controls_1to4_ST.mat');
load(controls_str);
dys_noIEDs_str = fullfile(pwd, 'dys_1to4_ST.mat');
load(dys_noIEDs_str);
all_subjects = [controls_1to4_ST; dys_1to4_ST]; % Concatenate!

addpath(genpath('script subfolder'));
analysis_params;
opt = [];
opt.metrics = {'coh_img','power','wpli_debiased'};
opt.scale = {'none'}; % no global scaling
opt.global = {'abs'}; % take absolute values
opt.dim_reduction = {'mean'}; % use the mean to reduce to one value per vertex
opt.save_grp_mean = true; % save a mean map
opt.hdm_class = 'ANTS9Years3T_fs7_openmeeg';
opt.output = fullfile(pwd,'export3_Dys_Controls_1to4_ST_source');
nmri_export_metrics(all_subjects, opt, params);

%% SUB EXPORT: DYSLEXIA AND CONTROLS, ONLY THOSE IN 5-6 & 7-8 SCHOOL GRADES (USING NEW PEDIATRIC TEMPLATE AND SHORT TRIALS!)
% (Controls 10.29-13.59 y.o. &  Dyslexia 7.02-10.95 y.o.)
clearvars;
pwd = 'pwd';
addpath(genpath('script subfolder'));
[all_subjects] = nmri_all_subjects(pwd);
filters.id = {'pseudonyms'};
[controls_5to8_ST] = nmri_filter_subjects(all_subjects, filters);
save('controls_5to8_ST.mat', 'controls_5to8_ST');
rmpath(genpath(pwd));

clearvars;
pwd = 'pwd';
addpath(genpath('script subfolder'));
[all_subjects] = nmri_all_subjects(pwd);
filters.id = {'pseudonyms'};
[dys_5to8_ST] = nmri_filter_subjects(all_subjects, filters);
save('dys_5to8_ST.mat', 'dys_5to8_ST');
rmpath(genpath(pwd));

clearvars;
pwd = 'main dir...';
controls_str = fullfile(pwd, 'controls_5to8_ST.mat');
load(controls_str);
dys_noIEDs_str = fullfile(pwd, 'dys_5to8_ST.mat');
load(dys_noIEDs_str);
all_subjects = [controls_5to8_ST; dys_5to8_ST]; % Concatenate!

addpath(genpath('script subfolder'));
analysis_params;
opt = [];
opt.metrics = {'coh_img','power','wpli_debiased'};
opt.scale = {'none'}; % no global scaling
opt.global = {'abs'}; % take absolute values
opt.dim_reduction = {'mean'}; % use the mean to reduce to one value per vertex
opt.save_grp_mean = true; % save a mean map
opt.hdm_class = 'ANTS9Years3T_fs7_openmeeg';
opt.output = fullfile(pwd,'export4_Dys_Controls_5to8_ST_source');
nmri_export_metrics(all_subjects, opt, params);

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