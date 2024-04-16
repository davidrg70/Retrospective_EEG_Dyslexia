% Created by David Garnica, david.garnica@med.uni-goettingen.de
% 2022, Universitätsmedizin Göttingen, Neurology Department

% dg_Dys_cognitiveData_SignRegions_averaged
clear all; close all; clc;
SCRIPT = 'SignRegions';

cd /main dir...
load([pwd '/data dir.../suma-all-fsaverage-10.mat']); % load SUMA for Desikan-Killiany
suma = suma_all;
params_dir = 'params dir';
run([params_dir '/analysis_params']);

analysis_dir = 'main dir';

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
if contains(metric_dir, '7_')
    export_file = [analysis_dir, '/group_filters/', '7_regressors_controls_dys_noIEDs.csv'];
    export = 'export_7';
elseif contains(metric_dir, '11_')
    export_file = [analysis_dir, '/group_filters/', '9_regressors_controls_dys_1to4.csv'];
    export = 'export_9';
end

% Open export file
if contains(export, '7')
    delimiter = ';';
elseif contains(export, '9')
    delimiter = ',';
end
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

%% Report significant vertices per frequency band and contrast

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

%% Create a logical mask of significant vertices in group comparison + report their number

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

%% Take RAW-VALUES from all subjects, only vertices+frequency bands with sign differences in group comparison

vertsThreshold = 10; % defines a threshold of a minimal number of significant vertices per region
% if exist('mask_c1','var')
%     for i = 1:size(mask_c1,2)
%         template_c1{i} = suma_all.annot(mask_c1(:,i)); % suma annotations of significant vertices! in contrast 1
%     end
% end
%
% if exist('mask_c2','var')
%     for i = 1:size(mask_c1,2)
%         template_c2{i} = suma_all.annot(mask_c2(:,i)); % suma annotations of significant vertices! in contrast 2
%     end
% end

% CONTRAST 1
if idx_c1 ~= 0
    RawMaskedSV_allSubjs_ff_c1 = []; reg_id_c1 = [];
    for ff = 1:length(idx_c1) % iterates across the frequency bands that showed significant differences between groups
        reg_lab_idx =  find(mask_c1(:,ff) == 1);   % positions in suma_all.annot of the significantly different vertices
        reg_id_c1{ff}  = suma_all.annot(reg_lab_idx);  % REGIONAL suma annotations of significant vertices, per frequency band!
        RawMaskedSV_allSubjs_c1 = []; % pre-allocate
        for hemiI=1:2
            for iID = 1:length(reg_id_c1{1,ff}) % iIDs refers to the vertices (reg_id refers to suma codes of the vertices that had significantly diff results)
                use_idcs = [];       % pre-allocate
                % finds vertices corresponding to every Region, Hemisphere, and MASK!
                use_idcs = find(suma_all.annot == reg_id_c1{1,ff}(iID) & suma_all.hemi==hemiI & mask_c1(:,ff))';
                
                tmp_ValsAll = []; % pre-allocate
                for iSubj = 1:(n_controls+n_patients) % iterates across subjects analyzed
                    tmp_ValsAll = squeeze(allRaws{1,idx_c1(ff)}(:,:,:,iSubj)); % All vertices' values from every subject (squeezed)
                    % INDEXES RAW CONNECTIVITY/POWER VALUES, PER SUBJECT!
                    temp_vertices = tmp_ValsAll(use_idcs,:);            % selects the values of those vertices belonging to that region!
                    RawMaskedSV_allSubjs_c1(iID, hemiI, iSubj) = nanmean(temp_vertices); % averages masked/significant vertices per region!
                end
                if length(use_idcs) < vertsThreshold % If region has < threshold sign vertices, save its annotation it to remove it later
                    regs_toRemove_c1(iID, hemiI) = double(reg_id_c1{1,ff}(iID)); % stores the region annotation
                end
            end
        end
        regs_toRemove_ff_c1{ff} = regs_toRemove_c1;
        RawMaskedSV_allSubjs_ff_c1{ff} = RawMaskedSV_allSubjs_c1;
    end
end

% CONTRAST 2
if idx_c2 ~= 0
    RawMaskedSV_allSubjs_ff_c2 = []; reg_id_c2 = [];
    for ff = 1:length(idx_c2) % iterates across the frequency bands that showed significant differences between groups
        reg_lab_idx = find(mask_c2(:,ff) == 1);   % positions in suma_all.annot of the significantly different vertices
        reg_id_c2{ff}  = suma_all.annot(reg_lab_idx);  % suma annotations of significant vertices, per frequency band!
        RawMaskedSV_allSubjs_c2 = []; % pre-allocate
        for hemiI=1:2
            for iID = 1:length(reg_id_c2{1,ff}) % iIDs refers to the vertices (reg_id refers to suma codes of the vertices that had significantly diff results)
                use_idcs = [];       % pre-allocate
                % finds vertices corresponding to every Region, Hemisphere, and MASK!
                use_idcs = find(suma_all.annot == reg_id_c2{1,ff}(iID) & suma_all.hemi==hemiI & mask_c2(:,ff))';
                
                tmp_ValsAll = []; % pre-allocate
                for iSubj = 1:120 % iterates across patients analyzed
                    tmp_ValsAll = squeeze(allRaws{1,idx_c2(ff)}(:,:,:,iSubj)); % All vertices' values from every subject (squeezed)
                    % INDEXES RAW CONNECTIVITY/POWER VALUES, PER SUBJECT!
                    temp_vertices = tmp_ValsAll(use_idcs,:);            % selects the values of those vertices belonging to that region!
                    RawMaskedSV_allSubjs_c2(iID, hemiI, iSubj) = nanmean(temp_vertices); % averages masked/significant vertices per region!
                end
                if length(use_idcs) < vertsThreshold % If region has < threshold sign vertices, save its annotation it to remove it later
                    regs_toRemove_c2(iID, hemiI) = double(reg_id_c2{1,ff}(iID)); % stores the region annotation
                end
            end
        end
        regs_toRemove_ff_c2{ff} = regs_toRemove_c2;
        RawMaskedSV_allSubjs_ff_c2{ff} = RawMaskedSV_allSubjs_c2;
    end
end

%% Take RAW-VALUES of patients and regions with >= threshold sign vertices
% SELECT DATA FROM PATIENTS ONLY
if idx_c1 ~= 0
    for i = 1:length(idx_c1)
        RawMaskedSV_Patients_ff_c1{i} = RawMaskedSV_allSubjs_ff_c1{1,i}(:,:,n_controls+1:end);
    end
end

if idx_c2 ~= 0
    for i = 1:length(idx_c2)
        RawMaskedSV_Patients_ff_c2{i} = RawMaskedSV_allSubjs_ff_c2{1,i}(:,:,n_controls+1:end);
    end
end

% Save DATA of regions with >= threshold sign vertices (both Hemispheres and Contrasts!)
% ONE AVERAGE PER REGION!
% CONTRAST 1
if idx_c1 ~= 0
    for i = 1:length(idx_c1)
        keys_c1_lh{i} = regs_toRemove_ff_c1{1,i}(:,1); % get indices of regional codes from LH regions to remove
    end
    regs_toRemove_ff_c1_lh = cellfun(@(x) unique(x),keys_c1_lh,'UniformOutput',0); % unique regional codes from regions to remove
    % ONE COLUMN PER FREQUENCY BAND
    
    for i = 1:length(idx_c1)
        keys_c1_rh{i} = regs_toRemove_ff_c1{1,i}(:,2); % get indices of regional codes from RH regions to remove
    end
    regs_toRemove_ff_c1_rh = cellfun(@(x) unique(x),keys_c1_rh,'UniformOutput',0); % unique regional codes from regions to remove
    % ONE COLUMN PER FREQUENCY BAND
    
    % Left Hemisphere
    PatientsData_lh = []; Patients_c1_lh =  [];
    for i = 1:length(idx_c1) % for loop to obtain data from every region >= threshold sign vertices, LH, and every patient
        PatientsData_lh = RawMaskedSV_Patients_ff_c1{1,i}(:,1,:); % SELECT DATA FROM LEFT HEMISPHERE
        PatientsData_lh = reshape(PatientsData_lh,[size(RawMaskedSV_Patients_ff_c1{1,i},1),n_patients,1]); % reshape to be able to cat regional codes with regional values
        PatientsData_lh = horzcat(double(reg_id_c1{1,i}), PatientsData_lh); % cat regional codes with regional values
        
        logsb = [];
        for k = 1:length(regs_toRemove_ff_c1_lh{1,i})
            logsb{k} = find(reg_id_c1{1,i} == regs_toRemove_ff_c1_lh{1,i}(k));
        end
        
        logsb = sort(vertcat(logsb{:})); % order those indices
        PatientsData_lh(logsb,:) = []; % REMOVES DATA FROM REGIONS WITH < THRESHOLD SIGNIFICANT VERTICES
        PatientsData_lh = unique(PatientsData_lh(:,:),'rows'); % unique regional codes and values (from all patients)
        Patients_c1_lh{i} = PatientsData_lh; % saves the LH data for every frequency band
    end
    
    % Right Hemisphere
    PatientsData_rh = []; Patients_c1_rh =  [];
    for i = 1:length(idx_c1) % for loop to obtain data from every region >= threshold sign vertices, RH, and every patient
        PatientsData_rh = RawMaskedSV_Patients_ff_c1{1,i}(:,2,:); % SELECT DATA FROM RIGHT HEMISPHERE
        PatientsData_rh = reshape(PatientsData_rh,[size(RawMaskedSV_Patients_ff_c1{1,i},1),n_patients,1]); % reshape to be able to cat regional codes with regional values
        PatientsData_rh = horzcat(double(reg_id_c1{1,i}), PatientsData_rh); % cat regional codes with regional values
        
        logsb = [];
        for k = 1:length(regs_toRemove_ff_c1_rh{1,i})
            logsb{k} = find(reg_id_c1{1,i} == regs_toRemove_ff_c1_rh{1,i}(k));
        end
        
        logsb = sort(vertcat(logsb{:})); % order those indices
        PatientsData_rh(logsb,:) = []; % REMOVES DATA FROM REGIONS WITH < THRESHOLD SIGNIFICANT VERTICES
        PatientsData_rh = unique(PatientsData_rh(:,:),'rows'); % unique regional codes and values (from all patients)
        Patients_c1_rh{i} = PatientsData_rh; % saves the RH data for every frequency band
    end
end

% CONTRAST 2
if idx_c2 ~= 0
    for i = 1:length(idx_c2)
        keys_c2_lh{i} = regs_toRemove_ff_c2{1,i}(:,1); % get indices of regional codes from LH regions to remove
    end
    regs_toRemove_ff_c2_lh = cellfun(@(x) unique(x),keys_c2_lh,'UniformOutput',0); % unique regional codes from regions to remove
    % ONE COLUMN PER FREQUENCY BAND
    
    for i = 1:length(idx_c2)
        keys_c2_rh{i} = regs_toRemove_ff_c2{1,i}(:,2); % get indices of regional codes from RH regions to remove
    end
    regs_toRemove_ff_c2_rh = cellfun(@(x) unique(x),keys_c2_rh,'UniformOutput',0); % unique regional codes from regions to remove
    % ONE COLUMN PER FREQUENCY BAND
    
    % Left Hemisphere
    PatientsData_lh = []; Patients_c2_lh =  [];
    for i = 1:length(idx_c2) % for loop to obtain data from every region >= threshold sign vertices, LH, and every patient
        PatientsData_lh = RawMaskedSV_Patients_ff_c2{1,i}(:,1,:); % SELECT DATA FROM LEFT HEMISPHERE
        PatientsData_lh = reshape(PatientsData_lh,[size(RawMaskedSV_Patients_ff_c2{1,i},1),70,1]); % reshape to be able to cat regional codes with regional values
        PatientsData_lh = horzcat(double(reg_id_c2{1,i}), PatientsData_lh); % cat regional codes with regional values
        
        logsb = [];
        for k = 1:length(regs_toRemove_ff_c2_lh{1,i})
            logsb{k} = find(reg_id_c2{1,i} == regs_toRemove_ff_c2_lh{1,i}(k));
        end
        
        logsb = sort(vertcat(logsb{:})); % order those indices
        PatientsData_lh(logsb,:) = []; % REMOVES DATA FROM REGIONS WITH < THRESHOLD SIGNIFICANT VERTICES
        PatientsData_lh = unique(PatientsData_lh(:,:),'rows'); % unique regional codes and values (from all patients)
        Patients_c2_lh{i} = PatientsData_lh; % saves the LH data for every frequency band
    end
    
    % Right Hemisphere
    PatientsData_rh = []; Patients_c2_rh =  [];
    for i = 1:length(idx_c2) % for loop to obtain data from every region >= threshold sign vertices, RH, and every patient
        PatientsData_rh = RawMaskedSV_Patients_ff_c2{1,i}(:,2,:); % SELECT DATA FROM RIGHT HEMISPHERE
        PatientsData_rh = reshape(PatientsData_rh,[size(RawMaskedSV_Patients_ff_c2{1,i},1),70,1]); % reshape to be able to cat regional codes with regional values
        PatientsData_rh = horzcat(double(reg_id_c2{1,i}), PatientsData_rh); % cat regional codes with regional values
        
        logsb = [];
        for k = 1:length(regs_toRemove_ff_c2_rh{1,i})
            logsb{k} = find(reg_id_c2{1,i} == regs_toRemove_ff_c2_rh{1,i}(k));
        end
        
        logsb = sort(vertcat(logsb{:})); % order those indices
        PatientsData_rh(logsb,:) = []; % REMOVES DATA FROM REGIONS WITH < THRESHOLD SIGNIFICANT VERTICES
        PatientsData_rh = unique(PatientsData_rh(:,:),'rows'); % unique regional codes and values (from all patients)
        Patients_c2_rh{i} = PatientsData_rh; % saves the RH data for every frequency band
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

% get correct agebytest
age_correct = struct;
for i = 1:length(subscores)
    pseuds = pseudonyms(:,i);
    pseuds(~any(cellfun('isempty',pseuds),2),:);                 % remove empty lines from strings arrays of pseudonyms
    [~,ia,~] = intersect(pseuds, export_patients);               % I obtain those patients in original data base that were actually in the export!
    age_correct.(subscores{i}) = agebytest.(subscores{i})(ia,1); % ia is index of agebytest to index from cog_data!
end     
agebytest = age_correct;

% get correct z-scores
z_scores = struct;
for i = 1:length(subscores)
    pseuds = pseudonyms(:,i);
    pseuds(~any(cellfun('isempty',pseuds),2),:);              % remove empty lines from strings arrays of pseudonyms
    [~,ia,~] = intersect(pseuds, export_patients);            % I obtain those patients in original data base that were actually in the export!
    z_scores.(subscores{i}) = cognitive.(subscores{i})(ia,1); % ia is index of z-scores to index from cog_data!
end                                                           % and in this way I get the z-scores of patients that only were exported

% save suma codes of regions!
if idx_c1 ~= 0
    for i = 1:size(Patients_c1_lh,2)
        regsCodes_c1_lh{i} = Patients_c1_lh{1,i}(:,1); % save the suma codes in variable!
        Patients_c1_lh{1,i}(:,1) = []; % and removes regional suma codes from database!
    end
    for i = 1:size(Patients_c1_rh,2)
        regsCodes_c1_rh{i} = Patients_c1_rh{1,i}(:,1);
    end
end

if idx_c2 ~= 0
    for i = 1:size(Patients_c2_lh,2)
        regsCodes_c2_lh{i} = Patients_c2_lh{1,i}(:,1); % save the suma codes in variable!
        Patients_c2_lh{1,i}(:,1) = []; % and removes regional suma codes from database!
    end
    for i = 1:size(Patients_c2_rh,2)
        regsCodes_c2_rh{i} = Patients_c2_rh{1,i}(:,1);
    end
end

% get correct power/FC values, per test, and hemisphere!
if idx_c1 ~= 0
    Patients_vals_c1_lh = struct; % Contrast 1 and LH
    for i = 1:size(Patients_c1_lh,2)
        Values = struct;
        for j = 1:length(subscores)
            [~,~,ib] = intersect(pseudonyms(:,j), export_patients);  % ib is index of power/FC values to index from all values averaged!
            Values.(subscores{j}) = Patients_c1_lh{1,i}(:,ib);
        end
        Patients_vals_c1_lh.(string(params.freqsNames(idx_c1(i)))) = Values;
    end
    
    Patients_vals_c1_rh = struct; % Contrast 1 and RH
    for i = 1:size(Patients_c1_rh,2)
        Values = struct;
        for j = 1:length(subscores)
            [~,~,ib] = intersect(pseudonyms(:,j), export_patients);  % ib is index of power/FC values to index from all values averaged!
            Values.(subscores{j}) = Patients_c1_rh{1,i}(:,ib);
        end
        Patients_vals_c1_rh.(string(params.freqsNames(idx_c1(i)))) = Values;
    end
end

if idx_c2 ~= 0
    Patients_vals_c2_lh = struct; % Contrast 2 and LH
    for i = 1:size(Patients_c2_lh,2)
        Values = struct;
        for j = 1:length(subscores)
            [~,~,ib] = intersect(pseudonyms(:,j), export_patients);  % ib is index of power/FC values to index from all values averaged!
            Values.(subscores{j}) = Patients_c2_lh{1,i}(:,ib);
        end
        Patients_vals_c2_lh.(string(params.freqsNames(idx_c2(i)))) = Values;
    end
    
    Patients_vals_c2_rh = struct; % Contrast 2 and RH
    for i = 1:size(Patients_c2_rh,2)
        Values = struct;
        for j = 1:length(subscores)
            [~,~,ib] = intersect(pseudonyms(:,j), export_patients);  % ib is index of power/FC values to index from all values averaged!
            Values.(subscores{j}) = Patients_c2_rh{1,i}(:,ib);
        end
        Patients_vals_c2_rh.(string(params.freqsNames(idx_c2(i)))) = Values;
    end
end

%% Correlations
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

% STOP HERE AND SELECT HEMISPHERE
fprintf('Please, first select HEMISPHERE! \n');
prompt = {'Hemisphere (1 Left, 2 Right)'};
dlgtitle = 'Select hemisphere to analyze its regions';
dims = [1 55];
answer = inputdlg(prompt, dlgtitle, dims);
[hemisphere] = deal(answer{:});

if contains(hemisphere, '1')
    hem_text = 'LEFT hemisphere';
elseif contains(hemisphere, '2')
    hem_text = 'RIGHT hemisphere';
end

if idx_c1 ~= 0 % CONTRAST 1
    if contains(correlation, '1') && contains(hemisphere, '1') % Pearson and LH
        R_c1 = struct; PVAL_c1 = struct;
        fn = fieldnames(Patients_vals_c1_lh);
        for i = 1:length(fn) % calculate correlations for every frequency band that showed significant group differences
            for j = 1:length(subscores)
                r = []; pval = [];
                for k = 1:size(Patients_vals_c1_lh.(string(params.freqsNames(idx_c1(i)))).(subscores{j}),1)
                    [r_t, pval_t] = corrcoef(Patients_vals_c1_lh.(string(params.freqsNames(idx_c1(i)))).(subscores{j})(k,:), (z_scores.(subscores{j})'));
                    r_t = r_t(1,2);       % only take relevant values
                    pval_t = pval_t(1,2); % only take relevant values
                    r{k} = r_t;
                    pval{k} = pval_t;
                end
                R_c1.(string(params.freqsNames(idx_c1(i)))).(subscores{j}) = r';
                PVAL_c1.(string(params.freqsNames(idx_c1(i)))).(subscores{j}) = pval';
            end
        end
    elseif contains(correlation, '1') && contains(hemisphere, '2') % Pearson and RH
        R_c1 = struct; PVAL_c1 = struct;
        fn = fieldnames(Patients_vals_c1_rh);
        for i = 1:length(fn) % calculate correlations for every frequency band that showed significant group differences
            for j = 1:length(subscores)
                r = []; pval = [];
                for k = 1:size(Patients_vals_c1_rh.(string(params.freqsNames(idx_c1(i)))).(subscores{j}),1)
                    [r_t, pval_t] = corrcoef(Patients_vals_c1_rh.(string(params.freqsNames(idx_c1(i)))).(subscores{j})(k,:), (z_scores.(subscores{j})'));
                    r_t = r_t(1,2);       % only take relevant values
                    pval_t = pval_t(1,2); % only take relevant values
                    r{k} = r_t;
                    pval{k} = pval_t;
                end
                R_c1.(string(params.freqsNames(idx_c1(i)))).(subscores{j}) = r';
                PVAL_c1.(string(params.freqsNames(idx_c1(i)))).(subscores{j}) = pval';
            end
        end
    elseif contains(correlation, '2') && contains(hemisphere, '1') % Spearman and LH
        Rho_c1 = struct; PVAL_c1 = struct;
        fn = fieldnames(Patients_vals_c1_lh);
        for i = 1:length(fn) % calculate correlations for every frequency band that showed significant group differences
            for j = 1:length(subscores)
                rho = []; pval = [];
                for k = 1:size(Patients_vals_c1_lh.(string(params.freqsNames(idx_c1(i)))).(subscores{j}),1)
                    [rho{k}, pval{k}] = corr((Patients_vals_c1_lh.(string(params.freqsNames(idx_c1(i)))).(subscores{j})(k,:)'), z_scores.(subscores{j}), 'Type','Spearman');
                end
                Rho_c1.(string(params.freqsNames(idx_c1(i)))).(subscores{j}) = cell2mat(rho)';
                PVAL_c1.(string(params.freqsNames(idx_c1(i)))).(subscores{j}) = cell2mat(pval)';
            end
        end
    elseif contains(correlation, '2') && contains(hemisphere, '2') % Spearman and RH
        Rho_c1 = struct; PVAL_c1 = struct;
        fn = fieldnames(Patients_vals_c1_rh);
        for i = 1:length(fn) % calculate correlations for every frequency band that showed significant group differences
            for j = 1:length(subscores)
                rho = []; pval = [];
                for k = 1:size(Patients_vals_c1_rh.(string(params.freqsNames(idx_c1(i)))).(subscores{j}),1)
                    [rho{k}, pval{k}] = corr((Patients_vals_c1_rh.(string(params.freqsNames(idx_c1(i)))).(subscores{j})(k,:)'), z_scores.(subscores{j}), 'Type','Spearman');
                end
                Rho_c1.(string(params.freqsNames(idx_c1(i)))).(subscores{j}) = cell2mat(rho)';
                PVAL_c1.(string(params.freqsNames(idx_c1(i)))).(subscores{j}) = cell2mat(pval)';
            end
        end
    elseif contains(correlation, '3') && contains(hemisphere, '1') % Partial and LH
        Rho_c1 = struct; PVAL_c1 = struct;
        fn = fieldnames(Patients_vals_c1_lh);
        for i = 1:length(fn) % calculate correlations for every frequency band that showed significant group differences
            for j = 1:length(subscores)
                rho = []; pval = [];
                for k = 1:size(Patients_vals_c1_lh.(string(params.freqsNames(idx_c1(i)))).(subscores{j}),1)
                    [rho{k}, pval{k}] = partialcorr((Patients_vals_c1_lh.(string(params.freqsNames(idx_c1(i)))).(subscores{j})(k,:)'),...
                        z_scores.(subscores{j}), agebytest.(subscores{j}), 'Type','Spearman');
                end
                Rho_c1.(string(params.freqsNames(idx_c1(i)))).(subscores{j}) = cell2mat(rho)';
                PVAL_c1.(string(params.freqsNames(idx_c1(i)))).(subscores{j}) = cell2mat(pval)';
            end
        end
    elseif contains(correlation, '3') && contains(hemisphere, '2') % Partial and RH
        Rho_c1 = struct; PVAL_c1 = struct;
        fn = fieldnames(Patients_vals_c1_rh);
        for i = 1:length(fn) % calculate correlations for every frequency band that showed significant group differences
            for j = 1:length(subscores)
                rho = []; pval = [];
                for k = 1:size(Patients_vals_c1_rh.(string(params.freqsNames(idx_c1(i)))).(subscores{j}),1)
                    [rho{k}, pval{k}] = partialcorr((Patients_vals_c1_rh.(string(params.freqsNames(idx_c1(i)))).(subscores{j})(k,:)'),...
                        z_scores.(subscores{j}), agebytest.(subscores{j}), 'Type','Spearman');
                end
                Rho_c1.(string(params.freqsNames(idx_c1(i)))).(subscores{j}) = cell2mat(rho)';
                PVAL_c1.(string(params.freqsNames(idx_c1(i)))).(subscores{j}) = cell2mat(pval)';
            end
        end        
    end
end

if idx_c2 ~= 0 % CONTRAST 2
    if contains(correlation, '1') && contains(hemisphere, '1') % Pearson and LH
        R_c2 = struct; PVAL_c2 = struct;
        fn = fieldnames(Patients_vals_c2_lh);
        for i = 1:length(fn) % calculate correlations for every frequency band that showed significant group differences
            for j = 1:length(subscores)
                r = []; pval = [];
                for k = 1:size(Patients_vals_c2_lh.(string(params.freqsNames(idx_c2(i)))).(subscores{j}),1)
                    [r_t, pval_t] = corrcoef(Patients_vals_c2_lh.(string(params.freqsNames(idx_c2(i)))).(subscores{j})(k,:), (z_scores.(subscores{j})'));
                    r_t = r_t(1,2);       % only take relevant values
                    pval_t = pval_t(1,2); % only take relevant values
                    r{k} = r_t;
                    pval{k} = pval_t;
                end
                R_c2.(string(params.freqsNames(idx_c2(i)))).(subscores{j}) = r';
                PVAL_c2.(string(params.freqsNames(idx_c2(i)))).(subscores{j}) = pval';
            end
        end
    elseif contains(correlation, '1') && contains(hemisphere, '2') % Pearson and RH
        R_c2 = struct; PVAL_c2 = struct;
        fn = fieldnames(Patients_vals_c2_rh);
        for i = 1:length(fn) % calculate correlations for every frequency band that showed significant group differences
            for j = 1:length(subscores)
                r = []; pval = [];
                for k = 1:size(Patients_vals_c2_rh.(string(params.freqsNames(idx_c2(i)))).(subscores{j}),1)
                    [r_t, pval_t] = corrcoef(Patients_vals_c2_rh.(string(params.freqsNames(idx_c2(i)))).(subscores{j})(k,:), (z_scores.(subscores{j})'));
                    r_t = r_t(1,2);       % only take relevant values
                    pval_t = pval_t(1,2); % only take relevant values
                    r{k} = r_t;
                    pval{k} = pval_t;
                end
                R_c2.(string(params.freqsNames(idx_c2(i)))).(subscores{j}) = r';
                PVAL_c2.(string(params.freqsNames(idx_c2(i)))).(subscores{j}) = pval';
            end
        end
    elseif contains(correlation, '2') && contains(hemisphere, '1') % Spearman and LH
        Rho_c2 = struct; PVAL_c2 = struct;
        fn = fieldnames(Patients_vals_c2_lh);
        for i = 1:length(fn) % calculate correlations for every frequency band that showed significant group differences
            for j = 1:length(subscores)
                rho = []; pval = [];
                for k = 1:size(Patients_vals_c2_lh.(string(params.freqsNames(idx_c2(i)))).(subscores{j}),1)
                    [rho{k}, pval{k}] = corr((Patients_vals_c2_lh.(string(params.freqsNames(idx_c2(i)))).(subscores{j})(k,:)'), z_scores.(subscores{j}), 'Type','Spearman');
                end
                Rho_c2.(string(params.freqsNames(idx_c2(i)))).(subscores{j}) = cell2mat(rho)';
                PVAL_c2.(string(params.freqsNames(idx_c2(i)))).(subscores{j}) = cell2mat(pval)';
            end
        end
    elseif contains(correlation, '2') && contains(hemisphere, '2') % Spearman and RH
        Rho_c2 = struct; PVAL_c2 = struct;
        fn = fieldnames(Patients_vals_c2_rh);
        for i = 1:length(fn) % calculate correlations for every frequency band that showed significant group differences
            for j = 1:length(subscores)
                rho = []; pval = [];
                for k = 1:size(Patients_vals_c2_rh.(string(params.freqsNames(idx_c2(i)))).(subscores{j}),1)
                    [rho{k}, pval{k}] = corr((Patients_vals_c2_rh.(string(params.freqsNames(idx_c2(i)))).(subscores{j})(k,:)'), z_scores.(subscores{j}), 'Type','Spearman');
                end
                Rho_c2.(string(params.freqsNames(idx_c2(i)))).(subscores{j}) = cell2mat(rho)';
                PVAL_c2.(string(params.freqsNames(idx_c2(i)))).(subscores{j}) = cell2mat(pval)';
            end
        end
        elseif contains(correlation, '3') && contains(hemisphere, '1') % Partial and LH
        Rho_c2 = struct; PVAL_c2 = struct;
        fn = fieldnames(Patients_vals_c2_lh);
        for i = 1:length(fn) % calculate correlations for every frequency band that showed significant group differences
            for j = 1:length(subscores)
                rho = []; pval = [];
                for k = 1:size(Patients_vals_c2_lh.(string(params.freqsNames(idx_c2(i)))).(subscores{j}),1)
                    [rho{k}, pval{k}] = partialcorr((Patients_vals_c2_lh.(string(params.freqsNames(idx_c2(i)))).(subscores{j})(k,:)'),...
                        z_scores.(subscores{j}), agebytest.(subscores{j}), 'Type','Spearman');
                end
                Rho_c2.(string(params.freqsNames(idx_c2(i)))).(subscores{j}) = cell2mat(rho)';
                PVAL_c2.(string(params.freqsNames(idx_c2(i)))).(subscores{j}) = cell2mat(pval)';
            end
        end    
    elseif contains(correlation, '3') && contains(hemisphere, '2') % Partial and RH
        Rho_c2 = struct; PVAL_c2 = struct;
        fn = fieldnames(Patients_vals_c2_rh);
        for i = 1:length(fn) % calculate correlations for every frequency band that showed significant group differences
            for j = 1:length(subscores)
                rho = []; pval = [];
                for k = 1:size(Patients_vals_c2_rh.(string(params.freqsNames(idx_c2(i)))).(subscores{j}),1)
                    [rho{k}, pval{k}] = partialcorr((Patients_vals_c2_rh.(string(params.freqsNames(idx_c2(i)))).(subscores{j})(k,:)'),...
                        z_scores.(subscores{j}), agebytest.(subscores{j}), 'Type','Spearman');
                end
                Rho_c2.(string(params.freqsNames(idx_c2(i)))).(subscores{j}) = cell2mat(rho)';
                PVAL_c2.(string(params.freqsNames(idx_c2(i)))).(subscores{j}) = cell2mat(pval)';
            end
        end        
    end    
end

%% FDR correction
% Calculate number of corrections to do (Apply FD to the Number of REGIONS with >= threshold significant vertices, FROM BOTH CONSTRASTS!)
% https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh
% [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals, q, method, report);

% Calculate FDR correction, accordingly
if idx_c1 ~= 0 & idx_c2 == 0     % In case that c1 showed significant differences but c2 did not
    fn = fieldnames(PVAL_c1);
    for i = 1:length(fn)
        for j = 1:length(subscores)
            [~, critP_c1.(string(params.freqsNames(idx_c1(i)))).(subscores{j}), ~, ...
                adjustedP_c1.(string(params.freqsNames(idx_c1(i)))).(subscores{j})] = ...
                fdr_bh(PVAL_c1.(string(params.freqsNames(idx_c1(i)))).(subscores{j}), 0.05, 'pdep', 'yes');
        end
    end
elseif idx_c1 == 0 & idx_c2 ~= 0    % In case that c1 showed no significant differences but c2 did
    fn = fieldnames(PVAL_c2);
    for i = 1:length(fn)
        for j = 1:length(subscores)
            [~, critP_c2.(string(params.freqsNames(idx_c2(i)))).(subscores{j}), ~, ...
                adjustedP_c2.(string(params.freqsNames(idx_c2(i)))).(subscores{j})] = ...
                fdr_bh(PVAL_c2.(string(params.freqsNames(idx_c2(i)))).(subscores{j}), 0.05, 'pdep', 'yes');
        end
    end
elseif idx_c1 ~= 0 & idx_c2 ~= 0    % In case that c1 and c2 have both significant group differences
    fn = fieldnames(PVAL_c1);
    for i = 1:length(fn)
        for j = 1:length(subscores)
            [~, critP_c1.(string(params.freqsNames(idx_c1(i)))).(subscores{j}), ~, ...
                adjustedP_c1.(string(params.freqsNames(idx_c1(i)))).(subscores{j})] = ...
                fdr_bh(PVAL_c1.(string(params.freqsNames(idx_c1(i)))).(subscores{j}), 0.05, 'pdep', 'yes');
        end
    end
    
    fn = fieldnames(PVAL_c2);
    for i = 1:length(fn)
        for j = 1:length(subscores)
            [~, critP_c2.(string(params.freqsNames(idx_c2(i)))).(subscores{j}), ~, ...
                adjustedP_c2.(string(params.freqsNames(idx_c2(i)))).(subscores{j})] = ...
                fdr_bh(PVAL_c2.(string(params.freqsNames(idx_c2(i)))).(subscores{j}), 0.05, 'pdep', 'yes');
        end
    end
end

%% Regions names, after SUMA annotations
if idx_c1 ~= 0 & ~isempty(regsCodes_c1_lh) & contains(hemisphere, '1') % Contrast 1 and LH
    regs_names = suma_all.annot_key{1,2}(1:76); % 1:76 refers to LH cortical regions
    for i = 1:size(regsCodes_c1_lh,2)
        regs_idxs_lh = []; regs_names_lh = [];
        for j = 1:length(regsCodes_c1_lh{1,i})
            regs_idxs_lh = find(suma_all.annot_key{1,1}(1:76) == regsCodes_c1_lh{1,i}(j));
            regs_names_lh{j} = regs_names(regs_idxs_lh);
        end
        regsNames_c1_lh{i} = regs_names_lh';
    end
    clear regs_names_lh;
elseif idx_c1 ~= 0 & ~isempty(regsCodes_c1_rh) & contains(hemisphere, '2') % Contrast 1 and RH
    regs_names = suma_all.annot_key{1,2}(101:176); % 101:176 refers to RH cortical regions
    for i = 1:size(regsCodes_c1_rh,2)
        regs_idxs_rh = []; regs_names_rh = [];
        for j = 1:length(regsCodes_c1_rh{1,i})
            regs_idxs_rh = find(suma_all.annot_key{1,1}(101:176) == regsCodes_c1_rh{1,i}(j));
            regs_names_rh{j} = regs_names(regs_idxs_rh);
        end
        regsNames_c1_rh{i} = regs_names_rh';
    end
    clear regs_names_rh;
elseif idx_c2 ~= 0 & ~isempty(regsCodes_c2_lh) & contains(hemisphere, '1') % Contrast 2 and LH
    regs_names = suma_all.annot_key{1,2}(1:76); % 1:76 refers to LH cortical regions
    for i = 1:size(regsCodes_c2_lh,2)
        regs_idxs_lh = []; regs_names_lh = [];
        for j = 1:length(regsCodes_c2_lh{1,i})
            regs_idxs_lh = find(suma_all.annot_key{1,1}(1:76) == regsCodes_c2_lh{1,i}(j));
            regs_names_lh{j} = regs_names(regs_idxs_lh);
        end
        regsNames_c2_lh{i} = regs_names_lh';
    end
    clear regs_names_lh;
elseif idx_c2 ~= 0 & ~isempty(regsCodes_c2_rh) & contains(hemisphere, '2') % Contrast 2 and RH
    regs_names = suma_all.annot_key{1,2}(101:176); % 101:176 refers to RH cortical regions
    for i = 1:size(regsCodes_c2_rh,2)
        regs_idxs_rh = []; regs_names_rh = [];
        for j = 1:length(regsCodes_c2_rh{1,i})
            regs_idxs_rh = find(suma_all.annot_key{1,1}(101:176) == regsCodes_c2_rh{1,i}(j));
            regs_names_rh{j} = regs_names(regs_idxs_rh);
        end
        regsNames_c2_rh{i} = regs_names_rh';
    end
    clear regs_names_rh;
end

fprintf('Correlations+FDR finished for %s, masked significant regions/vertices in %s. \n',  metric, hem_text);
fprintf('Run the script again to calculate Correlations+FDR for contralateral regions. \n');

%% Plots
% adjustedP = struct;
% adjustedP.adjustedP_c1 = adjustedP_c1;
% adjustedP.adjustedP_c2 = adjustedP_c2;
% [plots] = dg_masked_Plotting_forCorrelations_wholeGroup(SCRIPT, adjustedP, subscores, Values_c1_lh, idx_c1, idx_c2, cognitive, correlation, metric, params);

%% Report SignVertices per Metric, Contrast, FreqBand

% if metric == 'ImCoh' & correlation == '1' % theta c1 and alpha c2
%     left_theta = length(find(mask_c1 == 1 & suma_all.hemi == 1));
%     right_theta = length(find(mask_c1 == 1 & suma_all.hemi == 2));
%     left_alpha = length(find(mask_c2 == 1 & suma_all.hemi == 1));
%     right_alpha = length(find(mask_c2 == 1 & suma_all.hemi == 2));
% elseif metric == 'Power' & correlation == '1' % delta and theta, c1
%     left_delta = length(find(mask_c1(:,1) == 1 & suma_all.hemi == 1));
%     right_delta = length(find(mask_c1(:,1) == 1 & suma_all.hemi == 2));
%     left_theta = length(find(mask_c1(:,2) == 1 & suma_all.hemi == 1));
%     right_theta = length(find(mask_c1(:,2) == 1 & suma_all.hemi == 2));
% elseif metric == 'wPLI' & correlation == '1' % only theta, c1
%     left = length(find(mask_c1 == 1 & suma_all.hemi == 1));
%     right = length(find(mask_c1 == 1 & suma_all.hemi == 2));
% end

%% show number of vertices per sign. region
wtregion = [];
for a = 1:length(regsCodes_c1_lh{1,2})
    wtregion(a,:) = length(find(suma_all.annot == regsCodes_c1_lh{1,2}(a) & suma_all.hemi == 1 & mask_c1(:,2)));
end
