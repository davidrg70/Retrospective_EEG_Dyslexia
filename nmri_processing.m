function [ subject ] = nmri_processing(subject, data, params)
%[ subject ] = nmri_processing(subject, data, params)
%
% modified by DG, April 2024

% This function will do the source space processing including
% power @sensor level (for QC)
% source recon
% power @source
% connectivity @source
%
% subject   =   subject structure, see nmri_read_subject
%               could be a struct or .m file or .mat file
% data      =   will return the data (optional)

% written by NF 11/2016 - 06/2021
%
% 06/2021 - 07 /2021 added option to deal with canonical headmodel


% check the call
if (~exist('subject','var') )
    error('Need a valid subject struct or .m/.mat file to work with')
end

% call the subject and params include
nmri_include_read_ps


if ~exist('data','var') || isempty(data)
    % check if we have cleaned dataset (ICA or not)
    if (isfield(params,'useICA_clean') && params.useICA_clean==1)
        if (isfield(subject,'cleanICA_dataset') &&  exist(subject.cleanICA_dataset,'file'))
            input=subject.cleanICA_dataset;
            useICA=true;
        else
            error('ICA cleaned dataset not found - Run nmri_artifactrejection first')
        end
    else
        if (isfield(subject,'clean_dataset') &&  exist(subject.clean_dataset,'file'))
            input=subject.clean_dataset;
            useICA=false;
        else
            error('Cleaned dataset (w/o ICA) not found - Run nmri_artifactrejection first')
        end
    end
    % retain subject info, if different
    load(input,'data');
    if (~exist('data','var') )
        error('Could not load data')
    end
end


%% Get the modality-specific analysis params
[ params ] = nmri_get_modality_params( params, subject.dtype );

%% Check hdm class and electrodes

% now select the right hdm_class
if isfield(subject,'active_hdm_class')
    active_hdm_class=subject.active_hdm_class;
    % take hdm class from active type
    items=strsplit(active_hdm_class,'_');
    if length(items)>2
        params.headmodel=items{end};
    end
elseif isfield(params,'active_hdm_class')
    active_hdm_class=params.active_hdm_class;
    % take hdm class from active type
    items=strsplit(active_hdm_class,'_');
    if length(items)>2
        params.headmodel=items{end};
    end
else
    active_hdm_class=['individual_',params.headmodel];
end

% make sure we do not work with sensor space data
if strcmp(active_hdm_class,'sensor')
    error('Sensor space analyis was selected, but the source space script was called. Use nmri_processing_sensor instead')
end

% check if the have a specific hdm_class for this
if isfield(subject,'hdm_classes') && isfield(subject.hdm_classes,active_hdm_class)
    % then use these files from now on
    fprintf('Using the HDM specific info found for %s.\n',active_hdm_class)
    subject.electrodes_aligned=subject.hdm_classes.(active_hdm_class).electrodes_aligned;
    if isfield(subject.hdm_classes.(active_hdm_class),'hdm') && isfield(subject.hdm_classes.(active_hdm_class),'lead')
        subject.lead=subject.hdm_classes.(active_hdm_class).lead;
        subject.hdm=subject.hdm_classes.(active_hdm_class).hdm;
    else
        subject.hdm_lead=subject.hdm_classes.(active_hdm_class).hdm_lead;
    end
    subject.suma_surface=subject.hdm_classes.(active_hdm_class).suma_surface;
else
    % old style
    fprintf('No HDM specific info found for %s. Take what info is in the main struct.\n',active_hdm_class)
    % make empty struct
    subject.hdm_classes.(active_hdm_class)=[];
end



if (strcmp(subject.dtype,'EEG'))
    if (~isfield(subject,'electrodes_aligned') || ~exist(subject.electrodes_aligned,'file'))
        error('Have not found EEG electrodes - should be generated by the headmodel step')
    end
    load(subject.electrodes_aligned,'elec_aligned','elec_present','elec_missing')
end

%% make dir / files
if (~isfield(subject,'stats_dir'))
    subject.stats_dir=fullfile(subject.analysis_dir,subject.id,'stats');
end
if (~exist(subject.stats_dir,'dir'))
    mkdir(subject.stats_dir)
end

if isfield(subject,'hdm_classes') && isfield(subject.hdm_classes,active_hdm_class) && isfield(subject.hdm_classes.(active_hdm_class),'stats')
    % we can assume that this is the legal stats file
    subject.stats=subject.hdm_classes.(active_hdm_class).stats;
else
    % old style, or set new
    if (~isfield(subject,'stats') || isempty(regexp(subject.stats,[ active_hdm_class '.mat$'])))
        subject.stats=fullfile(subject.stats_dir,['source_stats_' subject.id '_' subject.exam_id '_' active_hdm_class '.mat']);
    end
    % remember hdm specific
    subject.hdm_classes.(active_hdm_class).stats=subject.stats;
end



% This is subject and exam_id specific (leadfield may differ), check if the
% currently requested headmodel type is used

% set and check for new style (hdm, lead seperate)
if (~isfield(subject,'hdm') || isempty(regexp(subject.hdm,[ params.headmodel '.mat$'])))
    % set to current standard
    warning(['Need to reset the hdm  file / struct for HDM=' params.headmodel '. Should usually not happen. Check and only proceed if you know what you are doing.'])
    if exist(fullfile(subject.analysis_dir,subject.id,'processed',['hdm_only_' subject.id '_' subject.exam_id '_' active_hdm_class '.mat']),'file')
        subject.hdm=fullfile(subject.analysis_dir,subject.id,'processed',['hdm_only_' subject.id '_' subject.exam_id '_' active_hdm_class '.mat']);
    else
        subject.hdm=fullfile(subject.analysis_dir,subject.id,'processed',['hdm_only_' subject.id '_' subject.exam_id '_' params.headmodel '.mat']);
    end
end

if (~isfield(subject,'lead') || isempty(regexp(subject.lead,[ params.headmodel '.mat$'])))
    % set to current standard
    %warning(['Need to reset the lead  file / struct for HDM=' params.headmodel '. Should usually not happen. Check and only proceed if you know what you are doing.'])
    if exist(fullfile(subject.analysis_dir,subject.id,'processed',['lead_only_' subject.id '_' subject.exam_id '_' active_hdm_class '.mat']),'file')
        subject.lead=fullfile(subject.analysis_dir,subject.id,'processed',['lead_only_' subject.id '_' subject.exam_id '_' active_hdm_class '.mat']);
    else
        subject.lead=fullfile(subject.analysis_dir,subject.id,'processed',['lead_only_' subject.id '_' subject.exam_id '_' params.headmodel '.mat']);
    end
end


% check for old style hdm-lead combined
if ~exist(subject.hdm,'file') || ~exist(subject.lead,'file')
    % set the hdm_lead then
    if (~isfield(subject,'hdm_lead') || isempty(regexp(subject.hdm_lead,[ params.headmodel '.mat$'])))
        % set to current standard
        warning(['Need to reset the hdm_lead  file / struct for HDM=' params.headmodel '. Should usually not happen. Check and only proceed if you know what you are doing.'])
        if exist(fullfile(subject.analysis_dir,subject.id,'processed',['hdm_lead_' subject.id '_' subject.exam_id '_' active_hdm_class '.mat']),'file')
            subject.hdm_lead=fullfile(subject.analysis_dir,subject.id,'processed',['hdm_lead_' subject.id '_' subject.exam_id '_' active_hdm_class '.mat']);
        else
            subject.hdm_lead=fullfile(subject.analysis_dir,subject.id,'processed',['hdm_lead_' subject.id '_' subject.exam_id '_' params.headmodel '.mat']);
        end
    end
    % and set hdm_lead for all
    subject.hdm=subject.hdm_lead;
    subject.lead=subject.hdm_lead;
end

fprintf('Using HDM file = %s\n',subject.hdm)
fprintf('Using leadfield file = %s\n',subject.lead)



%% Deal with trials selection
if (isfield(params,'nTrials'))
    nTrials=params.nTrials;
else
    nTrials=[]; % take all that are good
end

% check if we have a trial selection already
if (~isfield(subject,'SelectedTrials_file'))
    subject.SelectedTrials_file=fullfile(subject.analysis_dir,subject.id,'processed',['selected_trials_' subject.id '_' subject.exam_id '.mat']);
end
if (exist(subject.SelectedTrials_file,'file'))
    % if we have a file, load this, and be done
    disp('TrialSelection file found - loading')
    tsubj=load(subject.SelectedTrials_file,'subject');
    % check if correct nTrials (may have changed in the mean time)
    if isempty(nTrials) || nTrials==length(tsubj.subject.SelectedTrials)
        % seems still legal, take it
        fprintf('Reusing trial selection from %s\n',subject.SelectedTrials_file)
        subject.SelectedTrials=tsubj.subject.SelectedTrials;
    else
        warning(['There is a mismatch between trials in selection file and currently requested nTrials. Assuming that this is intentional and drawing new trials.'])
        if isfield(subject,'SelectedTrials')
            subject=rmfield(subject,'SelectedTrials');
        end
    end
end

if ~isfield(subject,'SelectedTrials') || isempty(subject.SelectedTrials)
    % determine trials now
    badTrials = [];
    goodTrials = [1:length(data.trial)]; % start with all good
    
    
    if (isfield(params,'rejectEvents') && params.rejectEvents == 1)
        % check in dws_filt dataset not to miss atypically run processing
        if ~isfield(subject,'evt_timings_seconds') && ~isfield(subject,'evt_markerFile_notFound')
            subject_clean=load(subject.dws_filt_dataset,'subject');
            items=fieldnames(subject_clean.subject);
            for i=1:length(items)
                if length(items{i})>4 && strcmp(items{i}(1:4),'evt_')
                    subject.(items{i})=subject_clean.subject.(items{i});
                end
            end
            if isfield(subject_clean.subject,'stamps') && isfield(subject_clean.subject.stamps,'readingevents')
                subject.stamps.readingevents=subject_clean.subject.stamps.readingevents;
            end
            clear subject_clean
        end
    end
    
    % call the central selection function now
    [ goodTrials, badTrials ] = nmri_trial_selector(subject,data,params);
    
    fprintf('\nTotal: GoodTrials, N=%d / BadTrials, N=%d\n',length(goodTrials),length(badTrials))
    subject.evt_goodTrials=goodTrials;
    subject.evt_badTrials=badTrials;
    
    if ~isempty(nTrials)
        % check if sufficient
        if (length(subject.evt_goodTrials)<nTrials)
            error(sprintf('Fewer good trials (%d) [after excluding events] in dataset than given in nTrials(%d)',length(subject.evt_goodTrials),nTrials))
        end
        % draw only from good trials then - make sure we are really random, use a
        % new rand stream
        s = RandStream('mt19937ar','Seed', seconds(round(milliseconds(second(datetime))*1000000)));
        subject.SelectedTrials = sort(datasample(s,subject.evt_goodTrials,nTrials,'Replace',false));
        fprintf('\nHave drawn N=%d random trials from all good trials now\n',length(subject.SelectedTrials))
    else
        % just take all good
        subject.SelectedTrials = goodTrials;
        fprintf('\nSelecting all good trials (N=%d) now\n',length(subject.SelectedTrials))
    end
    
    % and safe selection
    save(subject.SelectedTrials_file,'subject')
end

% we do not need trial markings any more now, remove to avoid Fieldtrip
% warnings
if isfield(data,'trial_markings')
    data=rmfield(data,'trial_markings');
end
if isfield(data,'trial_markings_sampleinfo')
    data=rmfield(data,'trial_markings_sampleinfo');
end

%% select just the wanted trials
cfg           = [];
cfg.trials    = subject.SelectedTrials;

% and also the non-bad channels - as of 01102019
if isfield(data,'bad_channels')
    good_channels={};
    for i=1:length(data.label)
        if ~any(strcmp(data.label{i},data.bad_channels))
            good_channels(end+1)=data.label(i);
        end
    end
    
    % also filter with elec_present and cross-check with headmodel, mismatch
    % can occur with more recent versions of EGI software having a E257 chanel
    % but this may also occur otherwise, if a channel was later
    % excluded/included or for canonical headmodels
    tmp=intersect(good_channels,elec_present.label','stable');
    if length(tmp)~=length(good_channels)
        warning(['There is a mismatch between electrodes present in the headmodel and data. Taking only the common channels.'])
        diffCh=[setdiff(good_channels,tmp) setdiff(tmp,good_channels)];
        fprintf('Difference: %s\n',diffCh{:})
        good_channels=tmp;
    end
    
    cfg.channel=good_channels;
    
    % we do not need bad_channels any more now, remove to avoid Fieldtrip
    % warnings
    if isfield(data,'bad_channels')
        data=rmfield(data,'bad_channels');
    end
    
end
data=ft_selectdata(cfg,data);


% Now demean finally and re-reference if EEG
if (strcmp(subject.dtype,'EEG'))
    cfg          = [];
    cfg.demean   = params.preproc_cfg.demean;
    cfg.reref      = 'yes';
    cfg.refchannel = 'all';
    data        = ft_preprocessing(cfg,data);
end

%% make the QC dir
if (~isfield(subject,'QCdir'))
    subject.QCdir=fullfile(subject.analysis_dir,subject.id,'QC');
end
if (~exist(subject.QCdir,'dir'))
    mkdir(subject.QCdir)
end


%% load the head model and leadfield


if ~isfield(subject,'hdm') || ~exist(subject.hdm,'file')
    error('Could not find headmodel - Make sure to run nmri_make_hdm_suma first')
end

if ~isfield(subject,'lead') || ~exist(subject.lead,'file')
    error('Could not find leadfield - Make sure to run nmri_make_hdm_suma first')
end

if ~isfield(subject,'suma_surface') || ~exist(subject.suma_surface,'file')
    error('Could not find SUMA surface - Make sure to run nmri_make_hdm_suma first')
end

% check if canonical HDM
hdmVars=who('-file',subject.lead);
if ismember('canonical_subject',hdmVars)
    load(subject.lead,'leadfield','canonical_subject');
else
    load(subject.lead,'leadfield');
end

load(subject.hdm,'hdm');

sumaVars=who('-file',subject.suma_surface);
if ismember('canonical_subject',sumaVars)
    load(subject.suma_surface,'suma_all');
    % check is match
    if exist('canonical_subject','var')
        s=load(subject.suma_surface,'canonical_subject');
        if ~strcmp(canonical_subject,s.canonical_subject)
            error(['Headmodel was reporting canonical subject=' canonical_subject ' but suma_all=' s.canonical_subject '. Likely a mismatch of HDM and SUMA.'])
        end
    else
        error(['Headmodel was reporting canonical subject=' canonical_subject ' but could not find info in suma_all. Likely a mismatch of HDM and SUMA.'])
    end
else
    load(subject.suma_surface,'suma_all');
end

% if we have canonical, make sure to have a canonical stats dir
if exist('canonical_subject','var')
    %  if isempty(regexp(subject.stats_dir,[ '_canon$'], 'once'))
    %   [pa,fi]=fileparts(subject.stats_dir);
    %   subject.stats_dir=fullfile(pa,[fi '_canon']);
    %  end
    %
    %
    % make a canonical file label
    ctxt=['canon_'];
    ctit=[ canonical_subject ' (canonical) - '];
    
    %if isempty(regexp(subject.stats,[ '_canon$'], 'once'))
    % subject.stats=fullfile(subject.stats_dir,['source_stats_canon_' ctxt subject.id '_' subject.exam_id '_' params.headmodel '.mat']);
    %end
else
    ctxt='';
    ctit='';
end


% check if leadfield, data and electrodes are consistant
if (length(leadfield.label)~=length(data.label)) || length(setdiff(leadfield.label,data.label))>0 || length(setdiff(data.label,leadfield.label))>0
    error('There is a mismatch between channels in the data and the leadfield. Maybe the data was changed (e.g. channels removed or added). Re-run the headmodel making (delete hdm_lead and electrodes_aligend [for EEG])')
end



%% Now do full Fourier analysis
% preallocate the freq structure
iFreq = 1;
fprintf('\nFourier transform - frequency: %s\n',params.freqsNames{iFreq})
cfg           = [];
cfg.method    = 'mtmfft';
cfg.output    = 'fourier';
cfg.pad       = 'nextpow2'; % recommended for speed
cfg.foi       = params.freqs(iFreq);
%cfg.trials    = subject.SelectedTrials;
cfg.tapsmofrq = params.tapsmofrq(iFreq);
tmp   = ft_freqanalysis(cfg, data);
sensor_fourier = repmat(tmp,1,length(params.freqs));

for iFreq = 2:length(params.freqs)
    fprintf('\nFourier transform - frequency: %s\n',params.freqsNames{iFreq})
    cfg.foi       = params.freqs(iFreq);
    cfg.tapsmofrq = params.tapsmofrq(iFreq);
    sensor_fourier(iFreq)   = ft_freqanalysis(cfg, data);
end


%% perform beamforming using DICS

% make the source frequency structure
source_fourier = sensor_fourier;
source_filters = {};

% calculate power and CSD in sensor space
for ff = 1:length(params.freqs);
    cfg           = [];
    cfg.method    = 'mtmfft';
    cfg.output    = 'powandcsd';
    cfg.pad       = 'nextpow2'; % recommended for speed
    cfg.foi       = sensor_fourier(ff).freq;
    % cfg.trials    = subject.SelectedTrials;
    cfg.tapsmofrq = params.tapsmofrq(ff);
    pow_and_csd_freq(ff)       = ft_freqanalysis(cfg, data);
end

if (isfield(params,'QC_plots') && params.QC_plots==1)
    % make plot
    % save QC images if not present
    if (~exist(fullfile(subject.QCdir,['pow_sens_plot_' subject.id '_' subject.exam_id '.png']),'file'))
        disp('Now doing sensor space power plots for QC')
        allf=length(params.freqs);
        cols=3;
        rows=ceil(allf/cols);
        hFig=figure('Position',[0,0,cols*300,rows*300],'Visible','off');
        if (~isfield(subject,'layout'))
            if (isfield(params,'layout'))
                subject.layout=params.layout;
            else
                error('no layout scheme specified in either subject or param')
            end
        end
        % get our layout
        cfg = [];
        cfg.layout=subject.layout;
        cfg.skipscale='yes';
        cfg.skipcomnt='yes';
        layout=ft_prepare_layout(cfg);
        layout.label=upper(layout.label); % make upper case
        
        cfg=[];
        cfg.gridscale = 200; %more beautiful
        cfg.layout    = layout;
        cfg.comment   = 'xlim';
        for i=1:allf
            subH=subplot(rows,cols,i);
            cfg.figure    = subH; % keep in same figure
            title(params.freqsNames{i},'FontSize',12,'FontWeight','bold')
            ft_topoplotTFR(cfg,pow_and_csd_freq(i))
        end
        saveas(hFig,fullfile(subject.QCdir,['pow_sens_plot_' subject.id '_' subject.exam_id '.png']),'png');
        
        % only show after save
        if isdeployed || ~usejava('desktop')
            % headless, do not show
            close(hFig)
        else
            set(hFig,'Visible','on')
        end
    end
end

%% Now do the source projection
for ff = 1:length(params.freqs)
    
    % make labels for SUMA sources
    source_fourier(ff).label = cellstr(num2str((1:sum(leadfield.inside))'));
    
    % construct an empty matrix for the frequency projection into source
    source_fourier(ff).fourierspctrm = zeros(size(sensor_fourier(ff).fourierspctrm,1),sum(leadfield.inside),size(sensor_fourier(ff).fourierspctrm,3));
    
    % check SUMA and leadfield N's
    if size(suma_all.pos,1)~=sum(leadfield.inside)
        warning('Not all SUMA points are inside the brain, there may be issues with vertex-dipole matching - investigate')
    end
    
    cfg                 = [];
    cfg.channel         = data.label;
    cfg.method          = 'dics';
    cfg.pad             = 'nextpow2'; % recommended for speed
    cfg.frequency       = pow_and_csd_freq(ff).freq;
    cfg.headmodel       = hdm;
    %cfg.grid            = leadfield;
    cfg.sourcemodel    = leadfield; % changed in recent Fieldtrip versions
    cfg.dics.lambda     = '5%';
    cfg.dics.keepfilter = 'yes';
    cfg.dics.fixedori   = 'yes';
    cfg.dics.realfilter = 'yes';
    cfg.dics.projectnoise  = 'yes';
    
    % for EEG need the electrodes
    if (strcmp(subject.dtype,'EEG'))
        cfg.elec            = elec_present;
    end
    
    source_filters{ff}     = ft_sourceanalysis(cfg, pow_and_csd_freq(ff));
    
    % now project the fourier spectrum into source space - NOTE: these are
    % only INSIDE dipoles. i.e there may be a mismatch to SUMA in rare cases
    ss=0;
    for i = 1:length(source_filters{ff}.avg.filter) % loop for each source/filter
        if ~isempty(source_filters{ff}.avg.filter{i})
            % for non-inside the filter will be empty, but still present
            ss=ss+1;
            source_fourier(ff).fourierspctrm(:,ss)=source_filters{ff}.avg.filter{i}*sensor_fourier(ff).fourierspctrm(:,:)';
        end
    end
    if (ss~=sum(leadfield.inside))
        error('Have not found the needed N of inside sources. This should not happen')
    end
    
end

% now safe the stats file (full or compact)
if (~isfield(params,'stats_compact') || params.stats_compact==0)
    if exist('canonical_subject','var')
        save(subject.stats,'sensor_fourier','source_fourier','source_filters','canonical_subject');
    else
        save(subject.stats,'sensor_fourier','source_fourier','source_filters');
    end
else
    source_filters_red=cell(1,length(source_filters));
    for i=1:length(source_filters)
        source_filters_red{i}.avg=source_filters{i}.avg;
        source_filters_red{i}.avg=removefields(source_filters_red{i}.avg,{'ori','eta','filter','label'});
    end
    source_filters=source_filters_red;
    if exist('canonical_subject','var')
        save(subject.stats,'source_filters','canonical_subject');
    else
        save(subject.stats,'source_filters');
    end
end


clear hdm data sensor_fourier

%% Plot power @sources for QC
if (isfield(params,'QC_plots') && params.QC_plots==1)
    % clear figure if still open
    if (exist('hFig','var') && ishandle(hFig))
        close(hFig);
    end
    % make plot
    % save QC images if not present
    powAsource=fullfile(subject.QCdir,['pow_source_plot_' ctxt subject.id '_' subject.exam_id '_' active_hdm_class '_' params.freqsNames{end} '.png']);
    if (~exist(powAsource,'file'))
        disp('Now doing source space power plots for QC')
        
        allf=length(params.freqs);
        
        opt=[];
        opt.format='png';
        opt.colormap='hot';
        opt.opathresh=3;
        opt.per_hemi=1;
        if (isfield(params,'ASEG_use') && params.ASEG_use==1)
            opt.per_cortex=1;
        else
            opt.per_cortex=0;
        end
        % load(fullfile(subject.analysis_dir,'scripts','roi_colors.mat'),'roi_cols');
        % opt.roi_color=roi_cols;
        for i=1:allf
            opt.output=fullfile(subject.QCdir,['pow_source_plot_' ctxt subject.id '_' subject.exam_id '_' active_hdm_class '_' params.freqsNames{i} '.png']);
            opt.title=['Normalized Power - ' ctit subject.id ' - ' subject.exam_id ' - ' params.headmodel ' - ' params.freqsNames{i} ' (' num2str(params.freqs(i)) '+/-' num2str(params.tapsmofrq(i)) ' Hz)'];
            % calculate power normalized by noise
            func=source_filters{i}.avg.pow./source_filters{i}.avg.noise;
            func(~suma_all.cortex)=NaN;
            hFig=nmri_plot_surface_suma(suma_all,func,opt);
            close(hFig);
        end
    end
end

clear source_filters

%% calculating connectivity metric by requested params
for mm=1:length(params.con_method)
    method = params.con_method{mm};
    opt=[];
    opt.format='png';
    opt.per_hemi=1;
    if (isfield(params,'ASEG_use') && params.ASEG_use==1)
        opt.per_cortex=1;
    else
        opt.per_cortex=0;
    end
    
    for ff = 1:length(params.freqs)
        cfg = [];
        cfg.method = method;
        if (strcmp(method, 'coh'))
            cfg.complex = 'complex';
        end
        fprintf('\nDoing connectivity analysis for: %s, Frequency: %s\n', method, params.freqsNames{ff})
        tmp = ft_connectivityanalysis(cfg,source_fourier(ff));
        if (strcmp(method, 'wpli_debiased'))
            metric = tmp.wpli_debiasedspctrm;
        elseif (strcmp(method, 'coh'))
            metric = tmp.cohspctrm;
        end
        
        % map to SUMA if needed
        if size(suma_all.pos,1)~=size(metric,1)
            warning('Mismatch of source-points and SUMA, will try to re-map')
            mapped_metric=zeros(size(suma_all.pos,1),size(suma_all.pos,1));
            mapped_metric(leadfield.inside,leadfield.inside)=metric;
            mapped_metric(~leadfield.inside,:)=NaN;
            mapped_metric(:,~leadfield.inside)=NaN;
            metric=mapped_metric;
            % check again
            if size(suma_all.pos,1)~=size(metric,1)
                error('Could not match size of SUMA and leadfield.inside. There seems something terribly wrong. Stopping here.')
            end
        end
        
        % save metric
        if (isfield(cfg,'complex') && strcmp(cfg.complex,'complex') && ~isreal(metric))
            % safe real
            cmetric=metric; % safe complex metric
            metric=real(cmetric);
            fprintf('Saving metric (real)...\n')
            save(fullfile(subject.stats_dir,[method '_real_' ctxt subject.id '_' subject.exam_id '_' active_hdm_class '_' params.freqsNames{ff} '.mat' ]),'metric')
            % Plot metric @sources for QC
            opt.output=fullfile(subject.QCdir,['con_' method '_real_' ctxt subject.id '_' subject.exam_id '_' active_hdm_class '_' params.freqsNames{ff} '.png']);
            opt.title=[ method ' (real/abs) - ' ctit subject.id ' - ' subject.exam_id ' - ' params.headmodel ' - ' params.freqsNames{ff} ' (' num2str(params.freqs(ff)) '+/-' num2str(params.tapsmofrq(ff)) ' Hz)'];
            if (isfield(params,'QC_plots') && params.QC_plots==1)
                % make plot
                % calculate only within mask
                msk_metric=metric;
                msk_metric(suma_all.msk==0,suma_all.msk==0)=NaN;
                func=nanmean(abs(msk_metric));
                func(suma_all.msk==0)=NaN;
                fprintf('Making surface plot (real)...\n')
                hFig=nmri_plot_surface_suma(suma_all,func,opt);
                close(hFig);
            end
            % safe img
            metric=imag(cmetric);
            fprintf('Saving metric (imaginary)...\n')
            save(fullfile(subject.stats_dir,[method '_img_' ctxt subject.id '_' subject.exam_id '_' active_hdm_class '_' params.freqsNames{ff} '.mat' ]),'metric')
            % Plot metric @sources for QC
            opt.output=fullfile(subject.QCdir,['con_' method '_img_' ctxt subject.id '_' subject.exam_id '_' active_hdm_class '_' params.freqsNames{ff} '.png']);
            opt.title=[ method ' (imaginary/abs) - ' ctit subject.id ' - ' ctxt subject.exam_id ' - ' params.headmodel ' - ' params.freqsNames{ff} ' (' num2str(params.freqs(ff)) '+/-' num2str(params.tapsmofrq(ff)) ' Hz)'];
            if (isfield(params,'QC_plots') && params.QC_plots==1)
                % make plot
                % calculate only within mask
                msk_metric=metric;
                msk_metric(suma_all.msk==0,suma_all.msk==0)=NaN;
                func=nanmean(abs(msk_metric));
                func(suma_all.msk==0)=NaN;
                fprintf('Making surface plot (imaginary)...\n')
                hFig=nmri_plot_surface_suma(suma_all,func,opt);
                close(hFig);
            end
        else
            % safe as is
            fprintf('Saving metric...\n')
            save(fullfile(subject.stats_dir,[method '_' ctxt subject.id '_' subject.exam_id '_' active_hdm_class '_' params.freqsNames{ff} '.mat' ]),'metric')
            % Plot metric @sources for QC
            opt.output=fullfile(subject.QCdir,['con_' method '_' ctxt subject.id '_' subject.exam_id '_' active_hdm_class '_' params.freqsNames{ff} '.png']);
            opt.title=[ method ' - ' subject.id ' - ' ctit subject.exam_id ' - ' params.headmodel ' - ' params.freqsNames{ff} ' (' num2str(params.freqs(ff)) '+/-' num2str(params.tapsmofrq(ff)) ' Hz)'];
            if (isfield(params,'QC_plots') && params.QC_plots==1)
                % make plot
                % calculate only within mask
                msk_metric=metric;
                msk_metric(suma_all.msk==0,suma_all.msk==0)=NaN;
                func=nanmean(msk_metric);
                func(suma_all.msk==0)=NaN;
                fprintf('Making surface plot...\n')
                hFig=nmri_plot_surface_suma(suma_all,func,opt);
                close(hFig);
            end
        end
    end
end

%% if we got this far, place a stamp for completed processing
subject=nmri_stamp_subject(subject,['processing_' active_hdm_class],params);

end % final end