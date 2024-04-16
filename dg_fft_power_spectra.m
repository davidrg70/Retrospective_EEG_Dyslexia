clear all; close all; clc;

% DEFINE FILTERS USED (during extraction and filtering -within pipeline!-)
pass_high = 0.5; % Hz
pass_low  = 30; % Hz

% check if Butterworth filters exist in params
analysis_params;
if params.preproc_cfg.hpfiltord == 1 && params.preproc_cfg.hpfiltord == 1
    if params.preproc_cfg.hpfreq > pass_high % checks which is higher (extraction high-pass or pipeline filter)
        pass_high_checked = params.preproc_cfg.hpfreq;
    elseif params.preproc_cfg.hpfreq <= pass_high
        pass_high_checked = pass_high;
    end
    
    if params.preproc_cfg.lpfreq > pass_low % checks which is higher (extraction low-pass or pipeline filter)
        pass_low_checked = pass_low;
    elseif params.preproc_cfg.lpfreq <= pass_low
        pass_low_checked = params.preproc_cfg.lpfreq;
    end
elseif params.preproc_cfg.hpfiltord ~= 1 || params.preproc_cfg.hpfiltord ~= 1
    warning('High-pass Butterworth filter not defined in params. Setting extraction high-pass to plot.');
    warning('Low-pass Butterworth filter not defined in params. Setting extraction low-pass to plot.');
else
    pass_high_checked = pass_high;
    pass_low_checked = pass_low;
end

analysis_dir1 = 'patients dir';
analysis_dir2 = 'controls dir';

% patients Dys_no_IEDs
patients_list = {'pseudonyms'};

% controls
controls_list = {'pseudonyms'};

%% PATIENTS data, fft, and averaging

for i = 1:length(patients_list)
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
    
    folderpath = fullfile(analysis_dir1, patients_list{i}, '/processed/'); % folderpath to upload SelectedTrials
    load([folderpath 'selected_trials_' patients_list{i} structure '.mat']); % loads subject struct
    selected_trials = subject.SelectedTrials; % array of actual trials brought to connectivity analyses
    clear subject;
    
    % load cleanICA data
    load([folderpath 'cleanICA_' patients_list{i} structure '.mat']);
    
    % concatenates trials selected
    eeg_data = cat(2, data.trial{1,selected_trials(:)});
    N = length(eeg_data);
    
    % check which are the bad channels!
    for a = 1:length(data.bad_channels)
        log_check(:,a) = contains(data.label, char(data.bad_channels(a)));
    end
    
    % checks positions of bad channels (and concatenates)
    for a = 1:length(data.bad_channels)
            where_bad(:,a) = find(log_check(:,a) == 1);
    end

    eeg_data(where_bad(:),:) = []; % EXCLUDES BAD CHANNELS FROM EEG DATA!
    
    % FFT!
    eeg_data = eeg_data'; % transpose, as MikeXCohen has con ANTS_Tuesday_prac_1_SOL.m (LINE 112)
    amplitude = 2*abs(fft(eeg_data,7500,1)/N); % FFT of amplitude (normalized with number of timepoints)
    amplitude = mean(amplitude,2);  % Averages all channels
    
    amplitudes_patients(i,:) = amplitude;
end

%% CONTROLS data, fft, and averaging
clearvars -except analysis_dir2 controls_list amplitudes_patients pass_high_checked pass_low_checked

for i = 1:length(controls_list)
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
    
    folderpath = fullfile(analysis_dir2, controls_list{i}, '/processed/'); % folderpath to upload SelectedTrials
    load([folderpath 'selected_trials_' controls_list{i} structure '.mat']); % loads subject struct
    selected_trials = subject.SelectedTrials; % array of actual trials brought to connectivity analyses
    clear subject;
    
    % load cleanICA data
    load([folderpath 'cleanICA_' controls_list{i} structure '.mat']);
    
    % concatenates trials selected
    eeg_data = cat(2, data.trial{1,selected_trials(:)});
    N = length(eeg_data);
    
    % check which are the bad channels!
    for a = 1:length(data.bad_channels)
        log_check(:,a) = contains(data.label, char(data.bad_channels(a)));
    end
    
    % checks positions of bad channels (and concatenates)
    for a = 1:length(data.bad_channels)
            where_bad(:,a) = find(log_check(:,a) == 1);
    end

    eeg_data(where_bad(:),:) = []; % EXCLUDES BAD CHANNELS FROM EEG DATA!
    
    % FFT!
    eeg_data = eeg_data'; % transpose, as MikeXCohen has con ANTS_Tuesday_prac_1_SOL.m (LINE 112)
    amplitude = 2*abs(fft(eeg_data,7500,1)/N); % FFT of amplitude (normalized with number of timepoints)
    amplitude = mean(amplitude,2);  % Averages all channels
    
    amplitudes_controls(i,:) = amplitude;
end

save('amplitudes_patients.mat', 'amplitudes_patients');
save('amplitudes_controls.mat', 'amplitudes_controls');

%% PLOTTING
load('amplitudes_patients.mat')
load('amplitudes_controls.mat');

fsample = 150;
N = 7500;

amplitudes_controls = amplitudes_controls';
amplitudes_patients = amplitudes_patients';

av_ampControls   = mean(amplitudes_controls,2);         % mean of power spectra, controls
var_ampControls1 = prctile(amplitudes_controls, 5);     % variance with 95% confidence interval, controls
var_ampControls2 = prctile(amplitudes_controls, 95);    % 5-95%

av_ampPatients   = mean(amplitudes_patients,2);         % mean of power spectra, patients
var_ampPatients1 = prctile(amplitudes_patients, 5);     % variance with 95% confidence interval, patients
var_ampPatients2 = prctile(amplitudes_patients, 95);    % 5-95%

hz = linspace(0,fsample,N); % create a frequencies vector
                            % freq sample / amplitude's length = frequency resolution
                            % 150 / 7500 = 0.02 Hz

% figure();
% plot(hz,amplitudes_patients)
% set(gca,'xlim',[0 40])
% xlabel('Frequency (Hz)')
% ylabel('Averaged amplitude (µV)')
% title(['Amplitude spectrum of patients, Dyslexia without IEDs'])
% 
% figure();
% plot(hz,amplitudes_controls)
% set(gca,'xlim',[0 40])
% xlabel('Frequency (Hz)')
% ylabel('Averaged amplitude (µV)')
% title(['Amplitude spectrum of Controls'])
% 
% figure();
% plot(hz,av_ampPatients)
% set(gca,'xlim',[0 40])
% xlabel('Frequency (Hz)')
% ylabel('Averaged amplitude (µV)')
% title(['Averaged amplitude spectrum of Dyslexia without IEDs'])
% 
% figure();
% plot(hz,av_ampControls)
% set(gca,'xlim',[0 40])
% xlabel('Frequency (Hz)')
% ylabel('Averaged amplitude (µV)')
% title(['Averaged amplitude spectrum of Controls'])

%% GROUP PLOT
close all
analysis_params;

smoothing_window = params.trial_length;
figure('Position', [ 10 10 2400 1600]);
hold on;

for ff = 1:length(params.freqsNames)
    x_val = params.freqs(ff)-params.tapsmofrq(ff):0.5:params.freqs(ff)+params.tapsmofrq(ff);
    y_val = ones([1 length(x_val)]) * 2;
    ar1 = area(x_val, y_val,'basevalue', 0, 'FaceColor',[0.1 0.1 0.1], 'LineStyle', '-');
    ar1.FaceAlpha = 0.04 + (ff*0.015); % transparency of background ("grey area")
    % ar1.EdgeColor = [0.1 0.1 0.1]*5;    
    ar1.EdgeColor = [1 1 1]; % thickness of lines of frequency bands (edges)
end

p_av_ampPatients = smoothdata(av_ampPatients, 'gaussian', smoothing_window)';
p_av_ampControls = smoothdata(av_ampControls, 'gaussian', smoothing_window)';

% make and plot patch field one (HC)
var_ampPatients1   = prctile(amplitudes_patients', 5);
p_var_ampPatients1 = smoothdata(var_ampPatients1, 'gaussian', smoothing_window);
var_ampPatients2   = prctile(amplitudes_patients', 95);
p_var_ampPatients2 = smoothdata(var_ampPatients2, 'gaussian', smoothing_window);
p1 = patch([hz fliplr(hz)], [p_var_ampPatients1 fliplr(p_var_ampPatients2)], 'r');  % plots the variance area
p1.EdgeAlpha = 0;       % transparency of edges
p1.FaceAlpha = 0.2;

% make and plot patch field two (Pat)
var_ampControls1 = prctile(amplitudes_controls', 5);
p_var_ampControls1 = smoothdata(var_ampControls1, 'gaussian', smoothing_window);
var_ampControls2 = prctile(amplitudes_controls', 95);
p_var_ampControls2 = smoothdata(var_ampControls2, 'gaussian', smoothing_window);
p2 = patch([hz fliplr(hz)], [p_var_ampControls1 fliplr(p_var_ampControls2)], 'b');  % plots the variance area
p2.EdgeAlpha = 0;       % transparency of edges
p2.FaceAlpha = 0.2;

p(1) = plot(hz,p_av_ampPatients, 'r', 'LineWidth', 1);
p(2) = plot(hz,p_av_ampControls, 'b', 'LineWidth', 1);
% set(gca,'xlim',[pass_high_checked pass_low_checked]); % SETTING PLOT LIMITS ACCORDING TO FILTERS!
set(gca,'xlim',[pass_high_checked pass_low_checked]); % SETTING PLOT LIMITS ACCORDING TO FILTERS!
set(gca,'ylim',[0 1.8]);
ax1 = gca;
ax.FontSize = 22;
% set(gca,'XTick',[pass_high_checked 4 8 12 20 29 40],'FontSize',22); % SETTING PLOT LIMITS ACCORDING TO FILTERS!
set(gca,'XTick',[pass_high_checked 2 4 8 12 20 29 40],'FontSize',22); % SETTING PLOT LIMITS ACCORDING TO FILTERS!
set(gca,'YTick',0.2:0.2:1.8,'FontSize',22);
xlabel('Frequency (Hz)','FontSize',26);
ylabel('Power (µV^{2}/Hz)','FontSize',26);
xline(0.5,'k--');
xline(1,'k--');
xline(2,'k--');
% title('Group averaged power at the sensor level');
position_leg = [0.726, 0.75, .1, .05];
% legend(p(:,:), 'Dyslexia', 'Controls', 'FontSize', 26, 'LineWidth', 1, 'Location', 'northeast');
legend(p(:,:), 'Dyslexia', 'Controls', 'FontSize', 26, 'LineWidth', 1, 'Position', position_leg);
t1 = text(1.80,1.74,'Delta', 'FontSize',26);
t2 = text(5,1.74,'Theta', 'FontSize',26);
t3 = text(9,1.74,'Alpha', 'FontSize',26);
t4 = text(15,1.74,'Beta 1', 'FontSize',26);
t5 = text(24,1.74,'Beta 2', 'FontSize',26);
% t6 = text(33,1.60,'Gamma', 'FontSize',18);