% Created by David Garnica, david.garnica@med.uni-goettingen.de
% 2022, Universitätsmedizin Göttingen, Neurology Department

% dg_Dys_cognitiveData_SignVertices_averaged
clear all; close all; clc;
SCRIPT = 'SignVertices';

cd /main dir...
load([pwd 'data dir.../conf/suma-all-fsaverage-10.mat']); % load SUMA for Desikan-Killiany
suma = suma_all;
params_dir = 'params dir...';
run([params_dir '/analysis_params']);

analysis_dir = 'main dir...';

% Define metric
fprintf('Select ANALYSIS and METRIC to work with! \n');
metric_dir = uigetdir(analysis_dir);

if contains(metric_dir, 'coh_img')
    metric = 'ImCoh'; fprintf('%s metric selected. \n', metric);
elseif contains(metric_dir, 'power')
    metric = 'Power'; fprintf('%s metric selected. \n', metric);
elseif contains(metric_dir, 'wpli_debiased')
    metric = 'wPLI'; fprintf('%s metric selected. \n', metric);
end

% GET EXPORT FILE, TO CHECK LATER EEG METRICS VALUES AND Z-SCORES MATCHING
export_file = [analysis_dir, '/group_filters/', '7_regressors_controls_dys_noIEDs.csv'];
export = 'export_7';

% Open export file
delimiter = ';';
startRow = 2;
formatSpec = '%s%C%f%C%[^\n\r]';
fileID = fopen(export_file,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
export_data = table(dataArray{1:end-1}, 'VariableNames', {'basename','dx_group','age','sex'});
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

data_subjects = string(table2array(export_data(:,1)));
data_diagnoses = string(table2array(export_data(:,2)));
log_patients = (contains(data_diagnoses, 'dyslexia'));
log_controls = (contains(data_diagnoses, 'control'));
n_patients = sum(log_patients(:) == 1);
n_controls = sum(log_controls(:) == 1);

%% Load all raw_data files

for ff = 1:length(params.freqsNames)
    fname_rawValues = [metric_dir, '/data/', params.freqsNames{ff}, '_merged.mgh']; % calls every mgz file (one for every frequency band)
    allRaws{ff} = load_mgh(fname_rawValues);
end

%% Load PALM results
%(loads connectivity/power values after group comparison), and check which ones have significant data

c1_text = 'Dyslexia < Controls';
c2_text = 'Dyslexia > Controls';

for ff = 1:length(params.freqsNames)
    fname_resVals_c1 = [metric_dir, '/palm_surface/', 'palm_out_tfce_tstat_fwep_m', num2str(ff), '_c1.mgz'];
    fname_resVals_c2 = [metric_dir, '/palm_surface/', 'palm_out_tfce_tstat_fwep_m', num2str(ff), '_c2.mgz'];
    
    resVals_c1{ff} = load_mgh(fname_resVals_c1);
    resVals_c2{ff} = load_mgh(fname_resVals_c2);
end

% logical check for constrast 1 and contrast 2
for ff = 1:length(params.freqsNames)
    check1{ff} = (resVals_c1{1,ff}(:,:,1,:) >(-log10(0.05)));
    check2{ff} = (resVals_c2{1,ff}(:,:,1,:) >(-log10(0.05)));
    
    if all(check1{1,ff} == 0)
        freqB_c1{ff} = 0; c1{ff} = 0;
    else
        freqB_c1{ff} = 1; c1{ff} = 1;
    end
    
    if all(check2{1,ff} == 0)
        freqB_c2{ff} = 0; c2{ff} = 0;
    else
        freqB_c2{ff} = 1; c2{ff} = 1;
    end
end

% get indices, of those frequency bands values that are significant
if (find(cell2mat(c1) == 1)) == (find(cell2mat(freqB_c1) == 1))
    idx_c1 = find(cell2mat(c1) == 1);
else
    idx_c1 = 0;
end

if (find(cell2mat(c2) == 1)) == (find(cell2mat(freqB_c2) == 1))
    idx_c2 = find(cell2mat(c2) == 1);
else
    idx_c2 = 0;
end

%% Create masks
% with previous indexes of significant differences

% inform about those frequency bands with significant differences per contrast
if idx_c1 ~= 0 & idx_c2 == 0
    for i = 1:length(idx_c1)
        fprintf('Contrast 1 (%s) showed significant results in %s \n', c1_text, string(params.freqsNames(idx_c1(i))));
    end
    fprintf('Contrast 2 showed no significant results! \n');
elseif idx_c2 ~= 0 & idx_c1 == 0
    for i = 1:length(idx_c2)
        fprintf('Contrast 2 (%s) showed significant results in %s \n', c2_text, string(params.freqsNames(idx_c2(i))));
    end
    fprintf('Contrast 1 showed no significant results! \n');
elseif idx_c1 ~= 0 & idx_c2 ~= 0
    for i = 1:length(idx_c1)
        fprintf('Contrast 1 (%s) showed significant results in %s \n', c1_text, string(params.freqsNames(idx_c1(i))));
    end
    for i = 1:length(idx_c2)
        fprintf('Contrast 2 (%s) showed significant results in %s \n', c2_text, string(params.freqsNames(idx_c2(i))));
    end
elseif (idx_c1 == 0) & (idx_c2 == 0)
    error('No significant results in the 2 contrasts, with 5 frequency bands each. \n')
end

% Obtain masks (logical binaries of those significant results above threshold)
if idx_c1 ~= 0
    for k = 1:length(resVals_c1)
        valuesSign_c1 = cell2mat(resVals_c1(idx_c1));
    end
    mask_c1 = valuesSign_c1 >=(-log10(0.05));
end

if idx_c2 ~= 0
    for k = 1:length(resVals_c2)
        valuesSign_c2 = cell2mat(resVals_c2(idx_c2));
    end
    mask_c2 = valuesSign_c2 >=(-log10(0.05));
end

% RAW VALUES (Indexation) IN CONTRAST 1, using mask, and Frequency bands with signficant results
% NOTE: by doing this, the order of the suma_annotations is lost!
if idx_c1 ~= 0
    for k = 1:length(idx_c1)
        maskedVertices_c1{k} = allRaws{1,idx_c1(k)}(mask_c1(:,k) == 1,1,1,n_controls+1:end);
    end
    for i = 1:size(maskedVertices_c1,2) % report number of significant vertices
        fprintf('%i significant vertices in %s in Contrast 1 \n', size(maskedVertices_c1{1,i},1), string(params.freqsNames(idx_c1(i))));
    end
end

% RAW VALUES (Indexation) IN CONTRAST 2, using mask, and Frequency bands with signficant results
% NOTE: by doing this, the order of the suma_annotations is lost!
if idx_c2 ~= 0
    for k = 1:length(idx_c2)
        maskedVertices_c2{k} = allRaws{1,idx_c2(k)}(mask_c2(:,k) == 1,1,1,n_controls+1:end);
    end
    for i = 1:size(maskedVertices_c2,2) % report number of significant vertices
        fprintf('%i significant vertices in %s in Contrast 2 \n', size(maskedVertices_c2{1,i},1), string(params.freqsNames(idx_c2(i))));
    end
end

% Now average raw & significant values, ONE average per subject - Contrast 1 (any number of significant frequency bands results')
if idx_c1 ~= 0
    for k = 1:length(idx_c1)
        for i = 1:70
            valuesSign_averaged_c1{k,i} = mean((maskedVertices_c1{1,k}(:,:,1,i)));
        end
    end
end

% Now average raw & significant values, ONE average per subject - Contrast 2 (any number of significant frequency bands results')
if idx_c2 ~= 0
    for k = 1:length(idx_c2)
        for i = 1:70
            valuesSign_averaged_c2{k,i} = mean((maskedVertices_c2{1,k}(:,:,1,i)));
        end
    end
end

%% Subjects correction
% WAIT! NOT ALL SUBJECTS WITH POWER/CONNECTIVITY VALUES (OR IN PALM ANALYSIS) HAD THE CORRESPONDING COGNITIVE TESTS!
% --> I programmed a correction here, indexing the correct subjects!
load('data dir.../Dys_noIEDs_pseudonyms_by_test.mat');
load('data dir.../Dys_noIEDs_cognitive_data.mat');
load('data dir.../Dys_noIEDs_age_by_test.mat');
subscores = {'FSIQ', 'VIQ', 'WMI', 'PSI', 'WordsR', 'PseudoWR', 'TextSR', 'SpellingW'};

export_subjects = table2cell(export_data(:,1));
export_controls = string(export_subjects(1:n_controls));
export_patients = string(export_subjects(n_controls+1:end));

% remove NaNs in cognitive subscores data
cognitive = struct;
for j = 1:size(cogdata,2)
    cognitive.(subscores{j}) = cogdata(:,j);     % fills field with every subscore data
    temp = cognitive.(subscores{j});             % logical
    cognitive.(subscores{j})(isnan(temp)) = [];  % removes NaN elements in every subscore field
end

% remove NaNs in age_by_test data
agebytest = struct;
for j = 1:size(cogdata_age_by_test,2)
    agebytest.(subscores{j}) = cogdata_age_by_test(:,j); % fills field with every subscore data
    temp = agebytest.(subscores{j});                     % logical
    agebytest.(subscores{j})(isnan(temp)) = [];          % removes NaN elements in every subscore field
end

% get correct z-scores
z_scores = struct;
for i = 1:length(subscores)
    pseuds = pseudonyms(:,i);
    pseuds(~any(cellfun('isempty',pseuds),2),:);              % remove empty lines from strings arrays of pseudonyms
    [~,ia,~] = intersect(pseuds, export_patients);            % I obtain those patients in original data base that were actually in the export!
    z_scores.(subscores{i}) = cognitive.(subscores{i})(ia,1); % ia is index of z-scores to index from cog_data!
end                                                           % and in this way I get the z-scores of patients that only were exported

% get correct power/FC values, per test
if idx_c1 ~= 0
    valuesSign_averaged_c1 = cell2mat(valuesSign_averaged_c1)';
    Values_c1 = struct;
    for i = 1:size(valuesSign_averaged_c1,2)
        Values = struct;
        for j = 1:length(subscores)
            [~,~,ib] = intersect(pseudonyms(:,j), export_patients);   % ib is index of power/FC values to index from all values averaged!
            Values.(subscores{j}) = valuesSign_averaged_c1(ib,i);     % index those that are only in both arrays!
        end                                                           % and in this way I get the power/FC values of patients that only were exported
        Values_c1.(string(params.freqsNames(idx_c1(i)))) = Values;    % stores correct values per frequency band that showed significant results
    end
end

if idx_c2 ~= 0
    valuesSign_averaged_c2 = cell2mat(valuesSign_averaged_c2)';
    Values_c2 = struct;
    for i = 1:size(valuesSign_averaged_c2,2)
        Values = struct;
        for j = 1:length(subscores)
            [~,~,ib] = intersect(pseudonyms(:,j), export_patients);  % ib is index of power/FC values to index from all values averaged!
            Values.(subscores{j}) = valuesSign_averaged_c2(ib,i);    % index those that are only in both arrays!
        end
        Values_c2.(string(params.freqsNames(idx_c2(i)))) = Values;   % stores correct values per frequency band that showed significant results
    end
end

%% Correlations (for whole group, disregard school grades)
% STOP HERE AND SELECT CORRELATION TO CALCULATE
warning('Select the correlation coefficient to calculate: ');
fprintf('Pearson, for parametric/interval variables. \n');
fprintf('Spearman, for non-parametric/ordinal variables. \n');
fprintf('Partial, to control for effect of age. \n');
prompt = {'1 Pearson, 2 Spearman, 3 Partial RANK'};
dlgtitle = 'Correlation coefficient to calculate';
dims = [1 50];
answer = inputdlg(prompt, dlgtitle, dims);
[correlation] = deal(answer{:});

if idx_c1 ~= 0 % CONTRAST 1
    R_c1 = struct; PVAL_c1 = struct;
    fn = fieldnames(Values_c1);
    for i = 1:length(fn) % calculate correlations for every frequency band that showed significant group differences
        for j = 1:length(subscores)
            if contains(correlation, '1')
                [r, pval] = corrcoef(Values_c1.(fn{i}).(subscores{j}), z_scores.(subscores{j}));
                r = r(1,2); pval = pval(1,2);            % only take relevant values
                PVAL_c1.(fn{i}).(subscores{j}) = pval;   % add values, in case there is more than one significant difference in frequency bands
                R_c1.(fn{i}).(subscores{j}) = r;
            elseif contains(correlation, '2')
                [rho, pval] = corr(Values_c1.(fn{i}).(subscores{j}), z_scores.(subscores{j}), 'Type','Spearman');
                PVAL_c1.(fn{i}).(subscores{j}) = pval;
                R_c1.(fn{i}).(subscores{j}) = rho;
            elseif contains(correlation, '3')
                [rho,pval] = partialcorr(Values_c1.(fn{i}).(subscores{j}), z_scores.(subscores{j}), agebytest.(subscores{j}), 'Type','Spearman');
                PVAL_c1.(fn{i}).(subscores{j}) = pval;
                R_c1.(fn{i}).(subscores{j}) = rho;
            end
        end
    end
end

if idx_c2 ~= 0 % CONTRAST 2
    R_c2 = struct; PVAL_c2 = struct;
    fn = fieldnames(Values_c2);
    for i = 1:length(fn) % calculate correlations for every frequency band that showed significant group differences
        for j = 1:length(subscores)
            if contains(correlation, '1')
                [r, pval] = corrcoef(Values_c2.(fn{i}).(subscores{j}), z_scores.(subscores{j}));
                r = r(1,2); pval = pval(1,2);            % only take relevant values
                PVAL_c2.(fn{i}).(subscores{j}) = pval;   % add values, in case there is more than one significant difference in frequency bands
                R_c2.(fn{i}).(subscores{j}) = r;
            elseif contains(correlation, '2')
                [rho, pval] = corr(Values_c2.(fn{i}).(subscores{j}), z_scores.(subscores{j}), 'Type','Spearman');
                PVAL_c2.(fn{i}).(subscores{j}) = pval;
                R_c2.(fn{i}).(subscores{j}) = rho;
            elseif contains(correlation, '3')
                [rho,pval] = partialcorr(Values_c2.(fn{i}).(subscores{j}), z_scores.(subscores{j}), agebytest.(subscores{j}), 'Type','Spearman');
                PVAL_c2.(fn{i}).(subscores{j}) = pval;
                R_c2.(fn{i}).(subscores{j}) = rho;
            end
        end
    end
end

%% FDR correction
% calculate number of corrections to do (number of frequency bands that showed significant group differences, FROM BOTH CONSTRASTS!)
% https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh
% [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals, q, method, report);

% the pvalues are unified and brought as a matrix to the FDR function
if idx_c1 ~= 0 & idx_c2 == 0                            % In case that c1 showed significant differences but c2 did not
    fn = fieldnames(PVAL_c1);
    PValues = [];
    for i = 1:length(fn)
        for j = 1:length(subscores)
            PValues{j,i} = PVAL_c1.(fn{i}).(subscores{j});
        end
    end
    PValues = cell2mat(PValues);
elseif idx_c1 == 0 & idx_c2 ~= 0                        % In case that c1 showed no significant differences but c2 did
    fn = fieldnames(PVAL_c2);
    PValues = [];
    for i = 1:length(fn)
        for j = 1:length(subscores)
            PValues{j,i} = PVAL_c2.(fn{i}).(subscores{j});
        end
    end
    PValues = cell2mat(PValues);
elseif idx_c1 ~= 0 & idx_c2 ~= 0                        % In case that c1 and c2 have both significant group differences,
    fn = fieldnames(PVAL_c1);
    werte1 = [];
    for i = 1:length(fn)
        for j = 1:length(subscores)
            werte1{j,i} = PVAL_c1.(fn{i}).(subscores{j});
        end
    end
    werte1 = cell2mat(werte1);
    
    fn = fieldnames(PVAL_c2);
    werte2 = [];
    for i = 1:length(fn)
        for j = 1:length(subscores)
            werte2{j,i} = PVAL_c2.(fn{i}).(subscores{j});
        end
    end
    werte2 = cell2mat(werte2);
    
    PValues = horzcat(werte1,werte2);
end

% Calculate FDR correction, accordingly
fprintf('--> FDR correction for the number of frequency bands with sign. group differences (C1, C2, or C1 & C2)! \n');
if idx_c1 ~= 0 & idx_c2 == 0
    for i = 1:length(subscores)
        [~, critP.(subscores{i}), adjCi.(subscores{i}), adjustedP.(subscores{i})] = fdr_bh(PValues(i,:), 0.05, 'pdep', 'yes');
    end
elseif idx_c1 == 0 & idx_c2 ~= 0
    for i = 1:length(subscores)
        [~, critP.(subscores{i}), adjCi.(subscores{i}), adjustedP.(subscores{i})] = fdr_bh(PValues(i,:), 0.05, 'pdep', 'yes');
    end
elseif idx_c1 ~= 0 & idx_c2 ~= 0
    for i = 1:length(subscores)
        [~, critP.(subscores{i}), adjCi.(subscores{i}), adjustedP.(subscores{i})] = fdr_bh(PValues(i,:), 0.05, 'pdep', 'yes');
    end
end

%% Data distribution plotting
% [plots] = dg_masked_Plotting_forCorrelations_wholeGroup(SCRIPT, adjustedP, subscores, Values_c1, idx_c1, idx_c2, cognitive, correlation, metric, params);
% fprintf('%s \n', plots);
