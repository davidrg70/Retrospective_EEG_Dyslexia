% DG violin plots function for power/FC global data, 10.2022 & 08.2023
% using https://github.com/bastibe/Violinplot-Matlab
clc; clear all; close all;

% Define the METRIC
main_dir = 'main dir';
addpath([main_dir, '/', 'violinplot_master']); % LOADS KEY FUNCTIONS TO PLOT VIOLINS

analysis_dir = '2nd level analysis dir';
run([analysis_dir, '/', 'analysis_params.m']); % LOADS PARAMS CONFIGURATION
fprintf('Please, first select ANALYSIS and METRIC to work with! \n');
fprintf('Within: 2nd level analysis dir \n');
metric_dir = uigetdir(analysis_dir);

% load data
if contains(metric_dir, 'coh_img')
    metric = 'ImCoh'; fprintf('%s metric selected. \n', metric);
elseif contains(metric_dir, 'power')
    metric = 'Power'; fprintf('%s metric selected. \n', metric);
elseif contains(metric_dir, 'wpli_debiased')
    metric = 'wPLI'; fprintf('%s metric selected. \n', metric);
end

freq_bands = params.freqsNames;
load([metric_dir, '/', 'nmri_palm_run_global_cf_removed_cfg.mat']); % opens configuration given to PALM, and gets relevant data
freq_str = {'Delta (1-4 Hz)', 'Theta (4-8 Hz)', 'Alpha (8-12 Hz)', 'Beta1 (12-20 Hz)', 'Beta2 (20-30 Hz)'};
% freq_str = data_labels;

global_data = [];
for i = 1:length(freq_bands)
    filename = fullfile(metric_dir, 'data', [freq_bands{i}, '_global.csv']);
    delimiter = {''};
    formatSpec = '%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    global_data{i} = dataArray{:,1};
    clearvars delimiter formatSpec fileID dataArray ans;
end
global_data = cell2mat(global_data);

%% Get p-values and t-stats given by PALM
% REMEMBER: frequencies are modalities in PALM's nomenclature
tstats = []; tstat = [];
for m = 1:length(freq_bands)
    for c = 1:length(contrast_labels)
        pfile = fullfile(metric_dir, 'palm_global', ['palm_out_dat_tstat_fwep', '_m', num2str(m), '_c', num2str(c), '.csv']); % defines files names to open them
        delimiter = {''};
        formatSpec = '%f%[^\n\r]';
        fileID = fopen(pfile,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
        fclose(fileID);
        tstat(c,:) = [dataArray{1:end-1}];
        clearvars filename delimiter formatSpec fileID dataArray ans;
    end
    tstats(m,:) = tstat;
end
clear tstat;

p_vals = [];
for m = 1:length(freq_bands)
    for c = 1:length(contrast_labels)
        p(c) = power(10, -(tstats(m,c)));
    end
    
    % take the most "significant one"
    if p(1) < p(2)
        p_vals(m) = p(1);
    elseif p(1) > p(2)
        p_vals(m) = p(2);
    end
end
clear p;

%%  Get Cohen's d given by PALM
cohens = [];
for m = 1:length(freq_bands)
    pfile = fullfile(metric_dir, 'palm_global', ['palm_out_dat_cohen', '_m', num2str(m), '_c', num2str(c), '.csv']); % defines files names to open them
    delimiter = {''};
    formatSpec = '%f%[^\n\r]';
    fileID = fopen(pfile,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    cohens(m) = abs([dataArray{1:end-1}]);
    clearvars filename delimiter formatSpec fileID dataArray ans;
end

%% Do FDR correction on p-values! (correct over frequencies-modalities)
[~, ~, ~, p_fdr] = fdr_bh(p_vals, 0.05, 'pdep', 'yes');

%% Violinplots
close all;

% determine min of min and max of max, to define x/y limits for plots
lowlim = min(global_data);
lowlim = min(lowlim);
if lowlim <= 0.02
    lowlim = 0.01;
end

% determine controls and patients dimensions!
reg_fields = fieldnames(reg);
logs_controls = reg.(reg_fields{2});
logs_patients = reg.(reg_fields{3});

num_controls = sum(logs_controls(logs_controls == 1));
num_patients = sum(logs_patients(logs_patients == 1));

figure('Renderer', 'painters', 'Position', [0 450 2200 900]);
hold on;
for k = 1:length(freq_bands)
    % define and prepare data matrix
    controls = global_data(logs_controls,k);
    patients = global_data(logs_patients,k);
    c_p = num_controls - num_patients;
    
    % check for different number of subjects in groups
    if c_p < 0  % if patients > controls
        diff = zeros(abs(c_p),1);
        controls = cat(1,controls, diff); % add difference as zeros to controls array
        controls(controls == 0) = NaN;    % but replace with NaN, to not plot those zeros
    elseif c_p >= 1  % if controls > patients
        diff = zeros(abs(c_p),1);
        patients = cat(1,patients, diff);   % add difference as zeros to patients array
        patients(patients == 0) = NaN;      % but replace with NaN, to not plot those zeros
    end
    
    matrix = cat(2,controls,patients); % defines the matrix to plot
    
    % categories
    expression = '(^|\.)\s*.';
    replace = '${upper($0)}';
    for j = 1:length(design_col_labels)
        str = design_col_labels{j};
        categories{j} = regexprep(str,expression,replace); % makes the first letter upper case
    end
    color = [0.3 0.4 0.8; 0.8 0.2 0.3];
    
    % run func to plot
    subplot(1,length(freq_bands),k);
    vplot = violinplot(matrix, categories, 'GroupOrder', categories, 'ViolinColor', color, ...
        'ShowWhiskers', true, 'ShowNotches', false,'ShowMean', false, 'ShowBox', true, 'ShowMedian', true);
    clear title; % THERE WAS A CONFLICT (variable/function)
    title([freq_str{k}], 'FontSize',14);
   
    % set proper labels and x/y ranges
    if contains('ImCoh', metric)
        upplim = max(global_data);
        for il = 1:length(upplim)                   % this for-loop determines a dynamic upper y-limit for ImCoh
            upplim(il) = round(upplim(il),2);
            upplim(il) = upplim(il)+(upplim(il)/5);
        end
        ylim([lowlim, upplim(k)]);
        xlim([0.5, 2.5]);
        ax1 = gca;
        ax1.FontSize = 16;
        set(gca,'XTick',[], 'FontSize',14);         % removes X-labels but leaves font number 18 for y-ticks
    elseif contains('Power', metric)
        upplim = max(global_data);
        for il = 1:length(upplim)                   % this for-loop determines a dynamic upper y-limit for Power
            upplim(il) = round(upplim(il),2);
            upplim(il) = upplim(il)+(upplim(il)/5);
        end
        ylim([lowlim, upplim(k)]);
        xlim([0.5, 2.5]);
        ax1 = gca;
        ax1.FontSize = 50;
        ax1.YAxis.Exponent = 0;
        set(gca,'XTick',[], 'FontSize',14);         % removes X-labels but leaves font number 18 for y-ticks    
    elseif contains('wPLI', metric)
        upplim = max(global_data);
        for il = 1:length(upplim)                   % this for-loop determines a dynamic upper y-limit for wPLI
            upplim(il) = round(upplim(il),2);
            upplim(il) = upplim(il)+(upplim(il)/5);
        end
        ylim([lowlim, upplim(k)]);
        xlim([0.5, 2.5]);
        ax1 = gca;
        ax1.FontSize = 16;
        set(gca,'XTick',[], 'FontSize',14);         % removes X-labels but leaves font number 18 for y-ticks
    end
    
    % Introduce p/d values according to PALM results ('Automatically')
    % by determining first which are the significant (uncorrected) p-values!
    % it will plot uncorrected and corrected ones, but just only if former are significant!
    if p_vals(k) < 0.001
        str1 = ['p<0.001', ',' ' d=', num2str(round(cohens(k),2))];
        text(0.82, upplim(k)-(upplim(k)*0.05), sprintf(str1), 'FontSize',16);
        str2 = ['**p(FDR)=', num2str(round(p_fdr(k),3))];
        text(0.84, upplim(k)-(upplim(k)*0.08), sprintf(str2), 'FontSize',16);
        errorbar(1.5, upplim(k)-(upplim(k)*0.12), 0.5, 'horizontal', 'color', [0.2, 0, 0], 'LineWidth', 2);
    elseif p_vals(k) > 0.001 && p_vals(k) < 0.05
        str1 = ['p=', num2str(round(p_vals(k),3)), ',' ' d=', num2str(round(cohens(k),2))];
        text(0.78, upplim(k)-(upplim(k)*0.03), sprintf(str1), 'FontSize',16);
        if p_fdr(k) < 0.05
            str2 = ['**p(FDR)=', num2str(round(p_fdr(k),3))];
            text(0.80, upplim(k)-(upplim(k)*0.06), sprintf(str2), 'FontSize',16);
            errorbar(1.5, upplim(k)-(upplim(k)*0.1), 0.5, 'horizontal', 'color', [0.2, 0, 0], 'LineWidth', 2);
        elseif p_fdr(k) >= 0.05
            str2 = ['p(FDR)=', num2str(round(p_fdr(k),3))];
            text(0.92, upplim(k)-(upplim(k)*0.06), sprintf(str2), 'FontSize',16);
        end
    elseif p_vals(k) > 0.05 && p_vals(k) < 0.1
        str1 = ['p=', num2str(round(p_vals(k),3)), ',' ' d=', num2str(round(cohens(k),2))];
        text(0.78, upplim(k)-(upplim(k)*0.05), sprintf(str1), 'FontSize',16);
        str2 = ['p(FDR)=', num2str(round(p_fdr(k),3))];
        text(0.92, upplim(k)-(upplim(k)*0.08), sprintf(str2), 'FontSize',16);    
    elseif p_vals(k) >= 0.1
        fprintf('p-value is not significant for %s, and then it was not plotted \n', freq_bands{k});
    end
   
    % add y label to the first
    if k == 1 && contains(metric, 'Power')
        metric_adj = 'Source-projected power (µV^{2})';
         ylabel(metric_adj, 'FontSize',20);
    elseif k == 1 && contains(metric, 'ImCoh')
        metric_adj = metric;
        ylabel(metric_adj, 'FontSize',20);
    elseif k == 1 && contains(metric, 'wPLI')
        metric_adj = metric;
        ylabel(metric_adj, 'FontSize',20);
    end
    
    % add legend to the last
    if k==length(freq_bands)
        if length(categories) == 2                                              % plot 2 legends if there are 2 groups compared
            handlevec = [vplot.ViolinPlot vplot.ViolinPlot2];
            legend(handlevec,categories, 'FontSize',18);
        end
        
        if length(categories) == 3                                              % plot 3 legends if there are 3 groups compared
            handlevec = [vplot.ViolinPlot vplot.ViolinPlot2 vplot.ViolinPlot3];
            legend(handlevec,categories, 'FontSize',18);
        end
    end
end
fprintf('Done with violin-plotting of %s \n', metric);