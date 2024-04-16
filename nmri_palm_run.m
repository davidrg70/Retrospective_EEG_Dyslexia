function cfg = nmri_palm_run(cfg)
%cfg = nmri_palm_run(cfg)

% modified by DG, April 2024

%   Function to run surface-based PALM in a compilable way
%   can deal with different analyis variants:
%     surface 3D (e.g. EEG/MEG source space analysis)
%       --> .mgh / .mgz data + surface .gii
%     surface 2D (e.g. EEG/MEG sensor space analysis), .mgh / .mgz data
%       --> .mgh / .mgz data + Fieldtrip layout
%     global 1D (e.g. global/averaged EEG/MEG analysis)
%       --> .csv data, no surface/layout
%
% not (yet) implemented
%      volume 3D
%
% input needs to be passed as a single cfg struct (good for flexible
% compilation)
%
%
% cfg              = config struct
%
% common to all analyis modes:
%  .data           = path to PALM-readable data files, can be a cell array
%                    (for multiple input modalities)
%                    usually .mgh from nmri_export_metrics
%  .data_labels    = optional, lables for each dataset, e.g. frequencies
%  .design         = design matrix for main analyis (without regessors)
%  .contrast       = contrast matrix for main analyis (without regessors)
%  .contrast_labels= cell array of contrast labels
%  .corrcon        = use correction over contrasts (PALM and plots)
%  .corrmod        = use correction over modalities (PALM and plots)
%  .subject_ids    = subject IDs (to match with regressors 'fname')
%  .reg            = regressors (path or struct)
%  .reg_use        = cell arry with info which regessor fields to use
%  .output         = base path for output
%  .concat_prefix  = prefix for concated output, default: 'concat_'
%  .noranktest     = whether to pass the '-noranktest' to PALM, default: 0
%  .exchb_groups   = exchangebility groups to pass to PALM, optional
%
% for surface 3D:
%  .surface        = path to a PALM readable surface file, e.g. suma_all.gii
%  .mask           = path to a PALM reable analyis mask, e.g. all_suma_msk.mgh
%  .plot_surface   = path or struct to a 3D nmri_suma_plot surface file,
%                    optional, will generate an auto-plot if given
%                    note, can be a higher ld version in which case, a
%                    auto-remapping will be done to the ref_surface
%  .ref_surface    = optional in case a different ld plot_surface is used.
%                    This needs to be identical to the analyis space
%                    (.surface)
%  .vis_fwe        = display range for plotting (-logp, FWE), default [1.3 3.5]
%  .vis_unc        = display range for plotting (-logp, uncp), default [1.3 4]
%
%
% for surface 2D / layout:
%  .layout         = path to a fieldtrip layout or loaded Fieldtrip layout
%                    struct. Surface and mask will be auto-generated
%  .use_channels   = Fieldtrip style selector of channels  to use from
%                    layout (see ft_channelselection), default: 'all'
%  .do_plot        = bool, make a plot or not
%
% for global 1D:
%  .do_plot        = 0/1, make a plot or notm default: 1
%  .global_viz_remove = cell array of confounds to remove from plot, default: =reg_use
%  .global_vis_log = plot logarithmic, default: 0
%  .global_vis_lock_scale = 0/1, fix the y-axis for multipel modalities/datasets
%  .global_vis_colors = special plot colors per group, cell array of Matlab
%                       color definition(s)
%  .viz_data_back  = Matlab plot colors per modalitiy / dataset, optional
%
%


if ~exist('cfg','var')
    error('Need a cfg variable')
end

if ~isfield(cfg,'data')
    error('Need dataset(s)')
else
    % deal with datsets
    if ischar(cfg.data)
        cfg.data={cfg.data};
    end
    
    if ~iscell(cfg.data)
        error('Datsets are needed as cell array');
    end
    
    % guess filetype of data
    [pa fi ext]=fileparts(cfg.data{1});
end

if ~isfield(cfg,'output')
    error('Need an output path')
end


% determine mode


if isfield(cfg,'surface') && ( strcmpi(ext,'.mgh') || strcmpi(ext,'.mgz'))
    pmode='surface';
    
    if ~exist(cfg.surface,'file')
        error('Cannot read surface file')
    end
    
    if ~isfield(cfg,'mask')
        error('Need a mask to define the analysis space')
    end
    
    if ~exist(cfg.mask,'file')
        error('Cannot read mask file')
    end
    
elseif isfield(cfg,'layout') && ( strcmpi(ext,'.mgh') || strcmpi(ext,'.mgz'))
    pmode='layout';
    
    % check the layout
    if ischar(cfg.layout) && exist(cfg.layout,'file')
        tmp=load(cfg.layout);
        cfg.layout=tmp.layout;
        clear tmp
    end
    
    if ~isstruct(cfg.layout) || ~isfield(cfg.layout,'label') || ~isfield(cfg.layout,'pos')
        error('Provided layout does not seem to be a valid Fieldtrip layout, or readable path to it')
    end
    
    if ~isfield(cfg,'use_channels')
        cfg.use_channels='all';
    end
    
    if ~isfield(cfg,'do_plot')
        cfg.do_plot=1;
    end
    
elseif strcmpi(ext,'.csv')
    pmode='global';
    
    if ~isfield(cfg,'do_plot')
        cfg.do_plot=1;
    end
    
else
    error('Could not detect the working mode of dataset.')
end

fprintf('Working mode: %s\n',pmode)


if ~isfield(cfg,'subject_ids')
    cfg.subject_ids=[];
end

if ~isfield(cfg,'concat_prefix')
    concatTxt='concat_';
else
    concatTxt=cfg.concat_prefix;
end


if ~isfield(cfg,'design') || ~isfield(cfg,'contrast')
    warning('No design and/or contrasts specified, assume a one-sample test')
    cfg.design=ones(length(cfg.all_subjects),1);
    cfg.contrast=[1;-1];
    cfg.contrast_labels={'Positive Main Effect','Negative Main Effect'};
end


if ~isfield(cfg,'reg')
    cfg.reg=[];
end

if ~isfield(cfg,'reg_use')
    cfg.reg_use={};
end

if ~isfield(cfg,'corrcon')
    cfg.corrcon=1;
end

if ~isfield(cfg,'corrmod')
    cfg.corrmod=1;
end

if ~isfield(cfg,'noranktest')
    cfg.noranktest=0;
end

if ~isfield(cfg,'exchb_groups')
    cfg.exchb_groups=[];
end

if ~isfield(cfg,'overwrite')
    cfg.overwrite=1; % overwrite as default, may unset for testing
end

if ischar(cfg.reg) && exist(cfg.reg,'file')
    % load regressor file
    cfg.reg=load(cfg.reg);
end

if ~isempty(cfg.reg) && ~isstruct(cfg.reg)
    error('Regressors need to be struct or empty')
end


if ~isfield(cfg,'permut')
    cfg.permut=5000;
end


if strcmp(pmode,'surface')
    if ~isfield(cfg,'vis_fwe')
        cfg.vis_fwe=[1.3 3.5];
    end
    
    if ~isfield(cfg,'vis_uncp')
        cfg.vis_uncp=[1.3 4];
    end
    
    if ~isfield(cfg,'plot_surface')
        % this is the SUMA-style 3D / source space plot
        cfg.plot_surface=[];
    end
    
    if ~isfield(cfg,'ref_surface')
        % this is the SUMA-style 3D / source space reference surface (same ld as
        % .surface), needed if resampling is requested
        cfg.ref_surface=[];
    end
    
    
    if ischar(cfg.plot_surface) && exist(cfg.plot_surface,'file')
        % load plot_surface file
        cfg.plot_surface=load(cfg.plot_surface);
    end
    
    if ~isempty(cfg.plot_surface) && ~isstruct(cfg.plot_surface)
        error('Plotting surface needs to be struct or empty')
    end
    
    if ischar(cfg.ref_surface) && exist(cfg.ref_surface,'file')
        % load ref_surface file (for mapping to plot_surface, optional)
        cfg.ref_surface=load(cfg.ref_surface);
    end
    
    
    
end

if strcmp(pmode,'layout')
    if ~isfield(cfg,'vis_fwe')
        cfg.vis_fwe=[1.3 3.5];
    end
    
    if ~isfield(cfg,'vis_uncp')
        cfg.vis_uncp=[1.3 4];
    end
    
    if ~isfield(cfg,'viz_data_back')
        cfg.viz_data_back={};
    end
    
    % fill up colors of datasets
    for i=1:length(cfg.data)
        if length(cfg.viz_data_back)<i || isempty(cfg.viz_data_back{i})
            cfg.viz_data_back{i}=['#FFFFFF']; % white as default
        end
    end
    
end


if ~isfield(cfg,'data_labels')
    % make empty labels
    cfg.data_labels=cell(length(cfg.data),1);
end

% check datasets
for i=1:length(cfg.data)
    if ~exist(cfg.data{i},'file')
        error(['Could not find dataset=' cfg.data{i}])
    else
        if isempty(cfg.data_labels{i})
            [pa fi ext]=fileparts(cfg.data{i});
            cfg.data_labels{i}=fi;
        end
    end
end

% check Ns
N=length(cfg.subject_ids);

if ~isempty(cfg.design) && size(cfg.design,1)~=N
    error('Design matrix does not match subject''s N')
end

% gobal viz presets
if strcmp(pmode,'global')
    if ~isfield(cfg,'global_vis_log')
        cfg.global_vis_log=0; % no log by default
    end
    
    if ~isfield(cfg,'global_vis_colors')
        cfg.global_vis_colors={[0.4 0.4 0.8],...
            [0.8 0.4 0.2],...
            [0.4 0.8 0.2],...
            [0.2 0.8 0.8],...
            [0.8 0.2 0.8],...
            [0.2 0.8 0.8]};
    end
    
    if ~isfield(cfg,'global_vis_log')
        cfg.design_col_labels={};
    end
    
    % fill up design col labes
    for i=1:size(cfg.design,2)
        if length(cfg.design_col_labels)<i || isempty(cfg.design_col_labels{i})
            cfg.design_col_labels{i}=['Group ' num2str(i)];
        end
    end
    
    if ~isfield(cfg,'viz_data_back')
        cfg.viz_data_back={};
    end
    
    % fill up colors of datasets
    for i=1:length(cfg.data)
        if length(cfg.viz_data_back)<i || isempty(cfg.viz_data_back{i})
            cfg.viz_data_back{i}=['#FFFFFF']; % white as default
        end
    end
    
    if ~isfield(cfg,'global_viz_remove')
        % default is to reomve effect of all plots
        cfg.global_viz_remove=cfg.reg_use;
    end
    
    if ~isfield(cfg,'global_viz_plot_labels')
        % shall we plot subject IDs?
        cfg.global_viz_plot_labels=0;
    end
    
    if ~isfield(cfg,'global_vis_lock_scale')
        % identical scaling for all datasets / frequencies
        cfg.global_vis_lock_scale=1; % default is yes
    end
    
    if ~isfield(cfg,'do_plot')
        cfg.do_plot=1;
    end
    
end

%% now prepare the analyis
% regressors
design_mat=[cfg.design zeros(N,length(cfg.reg_use))];
if ~isempty(cfg.reg)
    for n=1:N
        idx=find(strcmp(cfg.subject_ids{n},cfg.reg.fname));
        if ~isempty(idx)
            for r=1:length(cfg.reg_use)
                design_mat(n,size(cfg.design,2)+r)=cfg.reg.(cfg.reg_use{r})(idx);
            end
        else
            warning(['No regressor data found for subject=' cfg.subject_ids{n}])
            design_mat(n,size(cfg.design,2)+1:size(cfg.design,1)+length(cfg.reg_use))=nan(length(cfg.reg_use),1); % not found, so NaN
        end
    end
    % check for de-mean
    for r=1:length(cfg.reg_use)
        if ~islogical(cfg.reg.(cfg.reg_use{r})) && ~strncmp(cfg.reg_use{r},'grp_',4)
            % de-mean all non categorical / logical
            design_mat(:,size(cfg.design,2)+r)=design_mat(:,size(cfg.design,2)+r)-nanmean(design_mat(:,size(cfg.design,2)+r));
        end
    end
end

% make con
con_mat=zeros(size(cfg.contrast,1)+(length(cfg.reg_use)*2),size(design_mat,2));
% add external contrasts
con_mat(1:size(cfg.contrast,1),1:size(cfg.contrast,2))=cfg.contrast;
% now add a main effect per regressor
for r=1:length(cfg.reg_use)
    con_mat((r*2)-1+size(cfg.contrast,1),r+size(cfg.contrast,2))=1;
    con_mat((r*2)+size(cfg.contrast,1),r+size(cfg.contrast,2))=-1;
    % add a label
    txt=cfg.reg_use{r};
    txt=strrep(txt,'grp_','');
    cfg.contrast_labels=[cfg.contrast_labels;{[upper(txt(1)) txt(2:end) ' - Positive Effect']}];
    cfg.contrast_colormap=[cfg.contrast_colormap;{'spring'}];
    cfg.contrast_labels=[cfg.contrast_labels;{[upper(txt(1)) txt(2:end) ' - Negative Effect']}];
    cfg.contrast_colormap=[cfg.contrast_colormap;{'winter'}];
end

if ~exist(cfg.output,'dir')
    mkdir(cfg.output)
end



%% Now run PALM
palmdir=fullfile(cfg.output,['palm_' pmode]);
if ~exist(palmdir,'dir')
    mkdir(palmdir)
end

% setup for layout mode

if strcmp(pmode,'layout')
    % now make a PALM surface and mask of the layout
    selCh=ft_channelselection(cfg.use_channels,cfg.layout.label);
    
    % get index
    selIdx=cellfun(@(x) find(strcmp(x,cfg.layout.label)),selCh);
    
    % and make a selected channels layout
    layout=[];
    layout.pos=cfg.layout.pos(selIdx,:);
    layout.label=cfg.layout.label(selIdx,:);
    if isfield(cfg.layout,'width')
        layout.width=cfg.layout.width(selIdx,:);
    end
    if isfield(cfg.layout,'height')
        layout.height=cfg.layout.height(selIdx,:);
    end
    if isfield(cfg.layout,'mask')
        layout.mask=cfg.layout.mask;
    end
    if isfield(cfg.layout,'outline')
        layout.outline=cfg.layout.outline;
    end
    
    % triangulate
    trl=delaunayTriangulation(layout.pos);
    
    % and make a "conventional" 3D surface
    surf=[];
    surf.pos=[trl.Points,ones([length(selIdx),1])];
    surf.tri=trl.ConnectivityList;
    
    % cross check layout with data (Ns must agree), and make mask
    
    msk=ones([size(surf,1) 1]);
    for i=1:length(cfg.data)
        tmp=load_mgh(cfg.data{i});
        msk=msk&~any(isnan(squeeze(tmp)),2);
        if size(surf.pos,1)~=size(tmp,1)
            error('Mismatch of N selected channels and data. Check use_channels and layout.')
        end
    end
    clear tmp
    
    % write mask
    cfg.mask=fullfile(palmdir,'sensor_msk.mgh');
    nmri_write_mgh(cfg.mask,eye(4),{msk})
    
    
    hFig=figure('Position',[0 0 800 800],'Visible','off','Color','w');
    patch('Vertices', surf.pos, 'Faces', surf.tri,'FaceColor', 'none', 'EdgeColor', 'b');
    export_fig(hFig,fullfile(palmdir,'layout.png'),'-nocrop','-r200')
    close(hFig)
    cfg.surface=fullfile(palmdir,'layout_surface.gii');
    % now save as gifti
    ft_write_headshape(cfg.surface,surf,'format','gifti');
end


% check if we need to run and
% add each dataset / metric
pargs={};
do_run=0;
for i=1:length(cfg.data)
    pargs=[pargs {'-i',cfg.data{i}}];
    if length(cfg.data)>1
        mtxt=['m' num2str(i) '_'];
    else
        mtxt='';
    end
    for ci=1:size(con_mat,1)
        if strcmp(pmode,'surface')
            if ~exist(fullfile(palmdir,['palm_out_tfce_tstat_' mtxt 'c' num2str(ci) '.mgz']),'file')
                do_run=1;
            end
        elseif strcmp(pmode,'layout')
            if ~exist(fullfile(palmdir,['palm_out_tfce_tstat_' mtxt 'c' num2str(ci) '.mgz']),'file')
                do_run=1;
            end
        elseif strcmp(pmode,'global')
            if ~exist(fullfile(palmdir,['palm_out_dat_tstat_' mtxt 'c' num2str(ci) '.csv']),'file')
                do_run=1;
            end
        end
    end
end

if do_run
    fprintf('Runing PALM\n')
    pargs=[ pargs { '-d',fullfile(cfg.output,'design.csv'),'-t',fullfile(cfg.output,'contrast.csv'),...
        '-logp','-saveglm','-savedof','-accel','tail','-quiet',...
        '-n',num2str(cfg.permut),...
        '-o',fullfile(palmdir,'palm_out')}];
    
    if strcmp(pmode,'surface') || strcmp(pmode,'layout')
        % surface options (TFCE 2D, surface + mask)
        pargs=[ pargs { '-s',cfg.surface,'-m',cfg.mask, '-T','-tfce2d'}];
    elseif strcmp(pmode,'volume')
        % surface options (TFCE colume + mask)
        pargs=[ pargs { '-s',cfg.surface,'-m',cfg.mask, '-T'}];
    elseif strcmp(pmode,'global')
        % actually nothing special right now
    end
    
    if length(cfg.data)>1 && cfg.corrmod==1
        % add correct over modalities/datasets
        pargs=[ pargs { '-corrmod' }];
    end
    
    if size(cfg.contrast,1)>1 && cfg.corrcon==1
        % add correct over modalities/datasets
        pargs=[ pargs { '-corrcon' }];
    end
    
    if cfg.noranktest
        pargs=[ pargs { '-noranktest' }];
    end
    
    
    
    %unset freesurfer and FSL dir (PALM will look for .m commands there)
    setenv('FREESURFER_HOME','');
    setenv('FSLDIR','');
    
    
    % write out design + contrasts
    csvwrite(fullfile(cfg.output,'design.csv'),design_mat)
    csvwrite(fullfile(cfg.output,'contrast.csv'),con_mat)
    
    % write subject infos also
    nf_csvwrite(fullfile(cfg.output,'subjects.csv'),[cfg.subject_ids])
    
    % and write exchangeability, if so requested
    if ~isempty(cfg.exchb_groups)
        csvwrite(fullfile(cfg.output,'eb.csv'),cfg.exchb_groups)
        pargs=[ pargs { '-eb' , fullfile(cfg.output,'eb.csv')}];
        if size(cfg.exchb_groups,2) == 1
            % if just one column should be safe to be flexible
            pargs=[ pargs { '-within','-whole' }];
        end
    end
    
    
    
    
    palm(pargs{:})
else
    fprintf('PALM seems to be done already\n')
end


%% make surface viz
if isfield(cfg,'plot_surface') && ~isempty(cfg.plot_surface) && strcmp(pmode,'surface')
    fprintf('Now making plots and reporting results\n')
    
    remap_matrix=[];
    % check if the surface is identical and attempt remapping otherwise
    if isfield(cfg,'ref_surface') && ~isempty(cfg.ref_surface)
        [~, ~, remap_matrix.vertices, remap_matrix.weights ]=nmri_suma_surf2surf_transform(cfg.ref_surface.suma_all,cfg.plot_surface.suma_all, zeros([size(cfg.ref_surface.suma_all.pos,1),1]));
    else
        % no remapping, ref=plot
        cfg.ref_surface=cfg.plot_surface;
    end
    
    
    tbl={};
    
    % loop over thresholds
    pcons={{'fwep','FWE',0.05},{'uncp','uncorr.',0.05}};
    % add cross modality corrected
    
    if length(cfg.data)>1 && cfg.corrmod==1
        pcons=[pcons {{'mfwep','cross-modality FWE',0.05}}];
    end
    
    if size(cfg.contrast,1)>1 && cfg.corrcon==1
        pcons=[pcons {{'cfwep','cross-contrast FWE',0.05}}];
    end
    
    if length(cfg.data)>1 && size(cfg.contrast,1)>1 && cfg.corrmod==1 && cfg.corrcon==1
        pcons=[pcons {{'mcfwep','cross-modality & cross-contrast FWE',0.05}}];
    end
    
    
    for pci=1:length(pcons)
        
        viz_con={};
        viz_con_label={};
        viz_con_mod=[];
        viz_con_mc={};
        
        
        for mi=1:length(cfg.data)
            if length(cfg.data)>1
                mtxt=['m' num2str(mi) '_'];
            else
                mtxt='';
            end
            
            for ci=1:size(con_mat,1)
                viz_con=[viz_con {[pcons{pci}{1} '_' mtxt 'c' num2str(ci)]}];
                viz_con_label=[viz_con_label {[cfg.contrast_labels{ci} ' (' pcons{pci}{2}  ', p<' sprintf('%0.2f',pcons{pci}{3}) ')']} ];
                viz_con_mod=[viz_con_mod mi];
                viz_con_mc=[viz_con_mc {[mtxt 'c' num2str(ci)]}];
            end
        end
        
        
        for vi=1:length(viz_con)
            fprintf('Now making a plot for %s, %s\n',cfg.data_labels{viz_con_mod(vi)},viz_con_label{vi})
            pfile=fullfile(palmdir,['palm_out_tfce_tstat_' viz_con{vi} '.mgz']);
            if exist(pfile,'file')
                p_map=load_mgh(pfile);
                % get corresponding d-map and report
                dfile=fullfile(palmdir,['palm_out_dpv_cohen_' viz_con_mc{vi} '.mgz']);
                if exist(dfile,'file')
                    d_map=load_mgh(dfile);
                else
                    d_map=nan(size(p_map));
                end
                tbl=nmri_report_results(cfg.ref_surface.suma_all,p_map,pcons{pci}{3},d_map);
                save(fullfile(palmdir,['SUMA_results_' viz_con{vi} '.mat']),'tbl');
                nf_csvwrite(fullfile(palmdir,['SUMA_results_' viz_con{vi} '.csv']),tbl);
                
                % only show significant
                p_map(p_map<-log10(pcons{pci}{3}))=NaN; % zero out non-sig. vertices
                % set the range based on FWE or uncp
                if ~isempty(regexp(viz_con{vi},'uncp'))
                    lthr=cfg.vis_uncp(1);
                    uthr=cfg.vis_uncp(2);
                else
                    lthr=cfg.vis_fwe(1);
                    uthr=cfg.vis_fwe(2);
                end
                
                vcfg=[];
                vcfg.thresh=lthr;
                vcfg.fixalpha=0.7;  % fix to 70% alpha
                % extract contrast count from string
                cparts=strsplit(viz_con_mc{vi},'c');
                vcfg.colormap=cfg.contrast_colormap{str2num(cparts{end})};
                vcfg.clim=[lthr uthr];
                vcfg.per_hemi=1;
                vcfg.per_cortex=1;
                vcfg.output=fullfile(palmdir,['SUMA_fig_' viz_con{vi} '.png']);
                if ~exist(vcfg.output,'file') || cfg.overwrite==1
                    % check for remapping
                    if ~isempty(remap_matrix)
                        fprintf('Remapping results to %d vertices surface\n',size(cfg.plot_surface.suma_all.pos,1))
                        % fast reamppaing with pre-determined matrix
                        p_map=nmri_suma_surf2surf_transform(cfg.ref_surface.suma_all,cfg.plot_surface.suma_all, p_map, remap_matrix.vertices, remap_matrix.weights );
                    end
                    hFig=nmri_plot_surface_suma(cfg.plot_surface.suma_all,p_map,vcfg);
                    close(hFig)
                end
            else
                error(['Missing PALM file: ' pfile])
            end
        end
        
        % now merge per contrast
        for ci=1:size(con_mat,1)
            merge_all={};
            mtbl={};
            
            % make title
            txt=[cfg.contrast_labels{ci}];
            if isfield(cfg,'title') && ~isempty(cfg.title)
                txt=[cfg.title ': ' txt];
            end
            
            for mi=1:length(cfg.data)
                if length(cfg.data)>1
                    mtxt=['m' num2str(mi) '_'];
                else
                    mtxt='';
                end
                
                this_img=fullfile(palmdir,['SUMA_fig_' pcons{pci}{1} '_' mtxt 'c' num2str(ci) '.png']);
                this_img_lab=fullfile(palmdir,['SUMA_fig_' pcons{pci}{1} '_' mtxt 'c' num2str(ci) '_lab.png']);
                this_tbl=fullfile(palmdir,['SUMA_results_' pcons{pci}{1} '_' mtxt 'c' num2str(ci) '.mat']);
                
                % select images
                if exist(this_img,'file')
                    if ~exist(this_img_lab,'file') || cfg.overwrite==1
                        cmd=['convert ' this_img ' -font FreeSans -pointsize ' num2str(64) ' -fill black -background white -gravity SouthWest -annotate +16+60 "' cfg.data_labels{mi} '" ' this_img_lab];
                        system(cmd);
                    end
                    merge_all=[merge_all {this_img_lab}];
                end
                
                % now load and concat results table
                if exist(this_tbl,'file')
                    ctbl=load(this_tbl);
                    % add the dataset info and p-sig
                    ctbl.tbl=[[{'Analysis','Dataset','Threshold'};repmat({txt,cfg.data_labels{mi},[ pcons{pci}{2}  ', p<' sprintf('%0.2f',pcons{pci}{3})]},[size(ctbl.tbl,1)-1,1])],ctbl.tbl];
                    if isempty(mtbl)
                        mtbl=ctbl.tbl;
                        % add empty line
                        mtbl=[mtbl;repmat({''},[1,size(mtbl,2)])];
                    else
                        mtbl=[mtbl;ctbl.tbl(2:end,:)]; % merge without headline
                        % add empty line
                        mtbl=[mtbl;repmat({''},[1,size(mtbl,2)])];
                    end
                end
                
            end
            % now concat
            cout=fullfile(cfg.output,[concatTxt 'fig_' pcons{pci}{1} '_' legalize_label(cfg.contrast_labels{ci}) '.png']);
            tout=fullfile(cfg.output,[concatTxt 'results_' pcons{pci}{1} '_' legalize_label(cfg.contrast_labels{ci}) '.csv']);
            toutm=fullfile(cfg.output,[concatTxt 'results_' pcons{pci}{1} '_' legalize_label(cfg.contrast_labels{ci}) '.mat']);
            
            cmd=['montage -gravity South -mode concatenate -tile x' num2str(length(merge_all)) ' -background white ' sprintf('%s ',merge_all{:}) cout];
            system(cmd);
            % and annotate
            cmd=['convert ' cout ' -font FreeSans -fill black -background white -pointsize ' num2str(64) ' -gravity SouthEast -annotate +20+20 "' [ pcons{pci}{2}  ', p<' sprintf('%0.2f',pcons{pci}{3})]  '" -append  -pointsize ' num2str(100) ' label:''' txt ''' -gravity Center +swap -append ' cout];
            system(cmd);
            
            % and save table
            %save(toutm,'mtbl')
            nf_csvwrite(tout,mtbl);
            
        end
    end
    
    % now also plot Cohen's d maps
    viz_con={};
    viz_con_label={};
    viz_con_mod=[];
    
    % read mask
    msk=load_mgh(cfg.mask);
    
    for mi=1:length(cfg.data)
        if length(cfg.data)>1
            mtxt=['m' num2str(mi) '_'];
        else
            mtxt='';
        end
        
        for ci=1:size(con_mat,1)
            viz_con=[viz_con {['dpv_cohen_' mtxt 'c' num2str(ci)]}];
            viz_con_label=[viz_con_label {['Cohen''s d for ' cfg.contrast_labels{ci} ]} ];
            viz_con_mod=[viz_con_mod mi];
        end
    end
    
    
    for vi=1:length(viz_con)
        fprintf('Now making a plot for %s, %s\n',cfg.data_labels{viz_con_mod(vi)},viz_con_label{vi})
        d_map=load_mgh(fullfile(palmdir,['palm_out_' viz_con{vi} '.mgz']));
        % show only pos Cohen's d
        d_map(d_map<0)=0;
        % set masked vertices to NaN (make them transparent)
        d_map(msk==0)=NaN;
        vcfg=[];
        vcfg.fixalpha=0.7;  % fix to 70% alpha
        vcfg.clim=[0 1.2];
        vcfg.per_hemi=1;
        vcfg.per_cortex=1;
        vcfg.output=fullfile(palmdir,['SUMA_fig_' viz_con{vi} '.png']);
        % check for remapping
        if ~isempty(remap_matrix)
            fprintf('Remapping results to %d vertices surface\n',size(cfg.plot_surface.suma_all.pos,1))
            % fast reamppaing with pre-determined matrix
            d_map=nmri_suma_surf2surf_transform(cfg.ref_surface.suma_all,cfg.plot_surface.suma_all, d_map, remap_matrix.vertices, remap_matrix.weights );
        end
        hFig=nmri_plot_surface_suma(cfg.plot_surface.suma_all,d_map,vcfg);
        close(hFig)
    end
    
    % now merge per contrast
    for ci=1:size(con_mat,1)
        merge_all={};
        for mi=1:length(cfg.data)
            if length(cfg.data)>1
                mtxt=['m' num2str(mi) '_'];
            else
                mtxt='';
            end
            
            this_img=fullfile(palmdir,['SUMA_fig_dpv_cohen_' mtxt 'c' num2str(ci) '.png']);
            this_img_lab=fullfile(palmdir,['SUMA_fig_dpv_cohen_' mtxt 'c' num2str(ci) '_lab.png']);
            if exist(this_img,'file')
                cmd=['convert ' this_img ' -font FreeSans -pointsize ' num2str(64) ' -fill black -background white -gravity SouthWest -annotate +16+60 "' cfg.data_labels{mi} '" ' this_img_lab];
                system(cmd);
                merge_all=[merge_all {this_img_lab}];
            end
        end
        % now concat
        cout=fullfile(cfg.output,[concatTxt 'fig_cohen_' legalize_label(cfg.contrast_labels{ci}) '.png']);
        cmd=['montage -gravity South -mode concatenate -tile x' num2str(length(merge_all)) ' -background white ' sprintf('%s ',merge_all{:}) cout];
        system(cmd);
        % and annotate
        txt=[cfg.contrast_labels{ci}];
        if isfield(cfg,'title') && ~isempty(cfg.title)
            txt=[cfg.title ': ' txt];
        end
        cmd=['convert ' cout ' -font FreeSans -fill black -background white -pointsize ' num2str(64) ' -gravity SouthEast -annotate +20+20 "Cohen''s d" -append  -pointsize ' num2str(100) ' label:''' txt ''' -gravity Center +swap -append ' cout];
        system(cmd);
    end
    
end

%% make layout viz
if strcmp(pmode,'layout') && cfg.do_plot==1
    fprintf('Now making plots and reporting results for layout mode\n')
    
    tbl={};
    
    % loop over thresholds
    pcons={{'fwep','FWE',0.05},{'uncp','uncorr.',0.05}};
    % add cross modality corrected
    if length(cfg.data)>1 && cfg.corrmod==1
        pcons=[pcons {{'mfwep','cross-modality FWE',0.05}}];
    end
    
    if size(cfg.contrast,1)>1 && cfg.corrcon==1
        pcons=[pcons {{'cfwep','cross-contrast FWE',0.05}}];
    end
    
    if length(cfg.data)>1 && size(cfg.contrast,1)>1 && cfg.corrmod==1 && cfg.corrcon==1
        pcons=[pcons {{'mcfwep','cross-modality & cross-contrast FWE',0.05}}];
    end
    
    
    for pci=1:length(pcons)
        
        viz_con={};
        viz_con_label={};
        viz_con_mod=[];
        viz_con_mc={};
        
        
        for mi=1:length(cfg.data)
            if length(cfg.data)>1
                mtxt=['m' num2str(mi) '_'];
            else
                mtxt='';
            end
            
            for ci=1:size(con_mat,1)
                viz_con=[viz_con {[pcons{pci}{1} '_' mtxt 'c' num2str(ci)]}];
                viz_con_label=[viz_con_label {[cfg.contrast_labels{ci} ' (' pcons{pci}{2}  ', p<' sprintf('%0.2f',pcons{pci}{3}) ')']} ];
                viz_con_mod=[viz_con_mod mi];
                viz_con_mc=[viz_con_mc {[mtxt 'c' num2str(ci)]}];
            end
        end
        
        
        for vi=1:length(viz_con)
            fprintf('Now making a plot for %s, %s\n',cfg.data_labels{viz_con_mod(vi)},viz_con_label{vi})
            pfile=fullfile(palmdir,['palm_out_tfce_tstat_' viz_con{vi} '.mgz']);
            if exist(pfile,'file')
                p_map=load_mgh(pfile);
                % get corresponding d-map and report
                dfile=fullfile(palmdir,['palm_out_dpv_cohen_' viz_con_mc{vi} '.mgz']);
                if exist(dfile,'file')
                    d_map=load_mgh(dfile);
                else
                    d_map=nan(size(p_map));
                end
                
                tbl=nmri_report_results_sensor(layout,p_map,pcons{pci}{3},d_map);
                save(fullfile(palmdir,['Sensor_results_' viz_con{vi} '.mat']),'tbl');
                nf_csvwrite(fullfile(palmdir,['Sensor_results_' viz_con{vi} '.csv']),tbl);
                
                
                % only show significant
                p_map(p_map<-log10(pcons{pci}{3}))=0; % zero out non-sig. vertices
                % set the range based on FWE or uncp
                if ~isempty(regexp(viz_con{vi},'uncp'))
                    lthr=cfg.vis_uncp(1);
                    uthr=cfg.vis_uncp(2);
                else
                    lthr=cfg.vis_fwe(1);
                    uthr=cfg.vis_fwe(2);
                end
                
                
                % now plot map with Fieldtrip
                this_plot=fullfile(palmdir,['Sensor_fig_' viz_con{vi} '.png']);
                if ~exist(this_plot,'file') || cfg.overwrite==1
                    hFig=figure('Position',[0 0 800 800],'Visible','off','Color','w');
                    % colorbar only for last dataset
                    if viz_con_mod(vi)==length(cfg.data)
                        ft_topoplotTFR(struct('gridscale',200,'layout',cfg.layout,'comment','no','parameter','metric','marker','on','colorbar','yes','zlim',[lthr uthr]),...
                            struct('label',{layout.label},'dimord','chan','metric',p_map))
                    else
                        ft_topoplotTFR(struct('gridscale',200,'layout',cfg.layout,'comment','no','parameter','metric','marker','on','colorbar','no','zlim',[lthr uthr]),...
                            struct('label',{layout.label},'dimord','chan','metric',p_map))
                    end
                    export_fig(hFig,this_plot,'-nocrop','-r200')
                    close(hFig)
                end
            else
                error(['Missing PALM file: ' pfile])
            end
        end
        
        % now merge per contrast
        for ci=1:size(con_mat,1)
            merge_all={};
            mtbl={};
            
            % make title
            txt=[cfg.contrast_labels{ci}];
            if isfield(cfg,'title') && ~isempty(cfg.title)
                txt=[cfg.title ': ' txt];
            end
            
            for mi=1:length(cfg.data)
                if length(cfg.data)>1
                    mtxt=['m' num2str(mi) '_'];
                else
                    mtxt='';
                end
                
                this_img=fullfile(palmdir,['Sensor_fig_' pcons{pci}{1} '_' mtxt 'c' num2str(ci) '.png']);
                this_img_lab=fullfile(palmdir,['Sensor_fig_' pcons{pci}{1} '_' mtxt 'c' num2str(ci) '_lab.png']);
                this_tbl=fullfile(palmdir,['Sensor_results_' pcons{pci}{1} '_' mtxt 'c' num2str(ci) '.mat']);
                
                % select images
                if exist(this_img,'file')
                    if ~exist(this_img_lab,'file') || cfg.overwrite==1
                        % add freq title with same font
                        cmd=['convert ' this_img ' -shave 0x80 -font FreeSans-Bold -pointsize 64 -fill black -background "' cfg.viz_data_back{mi} '" -gravity North -splice 0x80 -annotate +0+2 "' cfg.data_labels{mi} '" ' this_img_lab];
                        %cmd=['convert ' this_img ' -font FreeSans -pointsize ' num2str(64) ' -fill black -background white -gravity SouthWest -annotate +16+60 "' cfg.data_labels{mi} '" ' this_img_lab];
                        system(cmd);
                    end
                    merge_all=[merge_all {this_img_lab}];
                end
                
                % now load and concat results table
                if exist(this_tbl,'file')
                    ctbl=load(this_tbl);
                    % add the dataset info and p-sig
                    ctbl.tbl=[[{'Analysis','Dataset','Threshold'};repmat({txt,cfg.data_labels{mi},[ pcons{pci}{2}  ', p<' sprintf('%0.2f',pcons{pci}{3})]},[size(ctbl.tbl,1)-1,1])],ctbl.tbl];
                    if isempty(mtbl)
                        mtbl=ctbl.tbl;
                        % add empty line
                        mtbl=[mtbl;repmat({''},[1,size(mtbl,2)])];
                    else
                        mtbl=[mtbl;ctbl.tbl(2:end,:)]; % merge without headline
                        % add empty line
                        mtbl=[mtbl;repmat({''},[1,size(mtbl,2)])];
                    end
                end
                
            end
            % now concat
            cout=fullfile(cfg.output,[concatTxt 'fig_' pcons{pci}{1} '_' legalize_label(cfg.contrast_labels{ci}) '.png']);
            tout=fullfile(cfg.output,[concatTxt 'results_' pcons{pci}{1} '_' legalize_label(cfg.contrast_labels{ci}) '.csv']);
            toutm=fullfile(cfg.output,[concatTxt 'results_' pcons{pci}{1} '_' legalize_label(cfg.contrast_labels{ci}) '.mat']);
            
            cmd=['montage -gravity South -tile ' num2str(length(merge_all)) 'x -background white -geometry +10+0 ' sprintf('%s ',merge_all{:}) cout];
            system(cmd);
            % and annotate
            cmd=['convert ' cout ' -font FreeSans -fill black -background white -pointsize ' num2str(64) ' -gravity SouthEast -annotate +20+20 "' [ pcons{pci}{2}  ', p<' sprintf('%0.2f',pcons{pci}{3})]  '" -append  -pointsize ' num2str(100) ' label:''' txt ''' -gravity Center +swap -append ' cout];
            system(cmd);
            
            % and save table
            %save(toutm,'mtbl')
            nf_csvwrite(tout,mtbl);
            
        end
    end
    %%
    % now also plot Cohen's d maps
    viz_con={};
    viz_con_label={};
    viz_con_mod=[];
    
    for mi=1:length(cfg.data)
        if length(cfg.data)>1
            mtxt=['m' num2str(mi) '_'];
        else
            mtxt='';
        end
        
        for ci=1:size(con_mat,1)
            viz_con=[viz_con {['dpv_cohen_' mtxt 'c' num2str(ci)]}];
            viz_con_label=[viz_con_label {['Cohen''s d for ' cfg.contrast_labels{ci} ]} ];
            viz_con_mod=[viz_con_mod mi];
        end
    end
    
    
    for vi=1:length(viz_con)
        fprintf('Now making a plot for %s, %s\n',cfg.data_labels{viz_con_mod(vi)},viz_con_label{vi})
        d_map=load_mgh(fullfile(palmdir,['palm_out_' viz_con{vi} '.mgz']));
        % show only pos Cohen's d
        d_map(d_map<0)=0;
        
        hFig=figure('Position',[0 0 800 800],'Visible','off','Color','w');
        % colorbar only for last dataset
        if viz_con_mod(vi)==length(cfg.data)
            ft_topoplotTFR(struct('gridscale',200,'layout',cfg.layout,'comment','no','parameter','metric','marker','on','colorbar','yes','zlim',[0 1.2]),...
                struct('label',{layout.label},'dimord','chan','metric',d_map))
        else
            ft_topoplotTFR(struct('gridscale',200,'layout',cfg.layout,'comment','no','parameter','metric','marker','on','colorbar','no','zlim',[0 1.2]),...
                struct('label',{layout.label},'dimord','chan','metric',d_map))
        end
        export_fig(hFig,fullfile(palmdir,['Sensor_fig_' viz_con{vi} '.png']),'-nocrop','-r200')
        close(hFig)
    end
    
    % now merge per contrast
    for ci=1:size(con_mat,1)
        merge_all={};
        for mi=1:length(cfg.data)
            if length(cfg.data)>1
                mtxt=['m' num2str(mi) '_'];
            else
                mtxt='';
            end
            
            this_img=fullfile(palmdir,['Sensor_fig_dpv_cohen_' mtxt 'c' num2str(ci) '.png']);
            this_img_lab=fullfile(palmdir,['Sensor_fig_dpv_cohen_' mtxt 'c' num2str(ci) '_lab.png']);
            if exist(this_img,'file')
                cmd=['convert ' this_img ' -shave 0x80 -font FreeSans-Bold -pointsize 64 -fill black -background "' cfg.viz_data_back{mi} '" -gravity North -splice 0x80 -annotate +0+2 "' cfg.data_labels{mi} '" ' this_img_lab];
                system(cmd);
                merge_all=[merge_all {this_img_lab}];
            end
        end
        % now concat
        cout=fullfile(cfg.output,[concatTxt 'fig_cohen_' legalize_label(cfg.contrast_labels{ci}) '.png']);
        cmd=['montage -gravity South -tile ' num2str(length(merge_all)) 'x -background white -geometry +10+0 ' sprintf('%s ',merge_all{:}) cout];
        system(cmd);
        % and annotate
        txt=[cfg.contrast_labels{ci}];
        if isfield(cfg,'title') && ~isempty(cfg.title)
            txt=[cfg.title ': ' txt];
        end
        cmd=['convert ' cout ' -font FreeSans -fill black -background white -pointsize ' num2str(64) ' -gravity SouthEast -annotate +20+20 "Cohen''s d" -append  -pointsize ' num2str(100) ' label:''' txt ''' -gravity Center +swap -append ' cout];
        system(cmd);
    end
    
    
    
end


%% make global viz
if strcmp(pmode,'global') && cfg.do_plot==1
    
    fprintf('Now making plots and reporting results\n')
    % make indices from design
    
    
    % now find the correct contrasts, i.e. which groups shall be compared
    Dgrps=size(cfg.design,2);
    pot=nchoosek(1:Dgrps,2);
    plot_contrasts={};
    plot_groups=[];
    % parse over all options and determine which one are contrasted
    for i=1:size(pot,1)
        % contrast to look for, 1 -1 mode
        cpot=zeros(1,size(con_mat,2));
        cpot(pot(i,1))=1;
        cpot(pot(i,2))=-1;
        % look if we have this contrast
        pos=find(all(con_mat==cpot,2));
        % and the inverse
        neg=find(all(con_mat==-cpot,2));
        if ~isempty(pos) && ~isempty(neg)
            plot_contrasts=[plot_contrasts {[pos,neg]}];
            if ~ismember(pot(i,1),plot_groups)
                plot_groups=[plot_groups pot(i,1)];
            end
            if ~ismember(pot(i,2),plot_groups)
                plot_groups=[plot_groups pot(i,2)];
            end
        else
            % not found in 1 -1 mode, try with single column
            cpot=zeros(1,size(con_mat,2));
            cpot(pot(i,1))=1;
            % look if we have this contrast
            pos=find(all(con_mat==cpot,2));
            % and the inverse
            neg=find(all(con_mat==-cpot,2));
            if ~isempty(pos) && ~isempty(neg)
                if ~any(cellfun(@(x) all(x==[pos,neg]),plot_contrasts))
                    plot_contrasts=[plot_contrasts {[pos,neg]}];
                    if ~ismember(pot(i,1),plot_groups)
                        plot_groups=[plot_groups pot(i,1)];
                    end
                    if ~ismember(pot(i,2),plot_groups)
                        plot_groups=[plot_groups pot(i,2)];
                    end
                end
            end
        end
    end
    
    
    % make filter to assign subjects to groups
    SubjFilt=zeros(size(cfg.design,1),1);
    plot_groups=sort(plot_groups);
    Ngrps=size(plot_groups,2);
    for i=1:Ngrps
        f=cfg.design(:,i)==1;
        if sum(f)>0
            % assign straight
            SubjFilt(f)=i;
            % we may be in singleton group mode, check for -1
            f=cfg.design(:,i)==-1;
            if sum(f)>0
                % assign straight
                SubjFilt(f)=i+1;
            end
        end
    end
    
    % load the values, assume that Matlab will handle CSVs
    mean_vals=cell(length(cfg.data),1);
    for i=1:length(cfg.data)
        mean_vals{i}=load(cfg.data{i});
    end
    
    % now subtract regressors (of no interest)
    for i=1:length(cfg.data)
        % find out which design colums are confounds
        for ci=size(cfg.design,2)+1:size(cfg.design,2)+length(cfg.reg_use)
            if any(strcmp(cfg.reg_use{ci-size(cfg.design,2)},cfg.global_viz_remove))
                % now find the matching contrast
                this_con=zeros([1,size(con_mat,2)]);
                this_con(ci)=1;
                fcon=find(all(con_mat==this_con,2));
                if ~isempty(fcon)
                    palmres_cope=fullfile(palmdir,['palm_out_dat_cope_m' num2str(i) '_c' num2str(fcon) '.csv']);
                    if exist(palmres_cope,'file')
                        cope=load(palmres_cope);
                        % now subtract this factor
                        mean_vals{i}=mean_vals{i}-(cope*design_mat(:,ci));
                    else
                        error(['Missing COPE file ' palmres_cope])
                    end
                else
                    error(['Could not find contrast for ' cfg.reg_use{ci-size(cfg.design,2)}])
                end
            end
        end
    end
    
    % estimate display range (across all datasets/frequencies)
    
    if cfg.global_vis_lock_scale || cfg.global_vis_log
        allmean=mean([mean_vals{:}],1);
        allstd=std([mean_vals{:}],1);
        lowTh=round(min(allmean-(2*allstd)),2);
        highTh=round(max(allmean+(3*allstd)),2);
    end
    
    % make sure low is plottable
    if cfg.global_vis_log
        logMin=round(min(allmean)/10,2); % 10-1 under lowest
        lowTh=max([lowTh,min(min([mean_vals{:}])),logMin]);
        % make sure we do not have negative data points
        for i=1:length(cfg.data)
            mean_vals{i}(mean_vals{i}<logMin)=logMin;
        end
    end
    %mcorr='m'; %use corrected p-vals between the freq bands / "modalities"
    mcorr=''; %use un-corrected p-vals between the freq bands / "modalities"
    
    mcmd=['montage -gravity East -mode concatenate -tile ' num2str(length(cfg.data)) 'x -geometry +15+10 -background white '];
    for i=1:length(cfg.data)
        % plot
        
        hFig=figure('Position',[0 0 400 800],'Visible','off','Color','w');
        %hFig=figure('Position',[0 0 400 800],'Color','w','Visible','on');
        hV=violinplot(mean_vals{i}, SubjFilt,'ViolinAlpha',0.4,'BoxColor',[0.2 0.2 0.2]);
        % make a grey background
        hFig.CurrentAxes.Color=[0.8 0.8 0.8];
        hFig.CurrentAxes.YGrid='on';
        hFig.CurrentAxes.XGrid='on';
        hFig.CurrentAxes.XTickLabel=[];
        hFig.CurrentAxes.GridColor='w';
        hFig.CurrentAxes.LineWidth=1;
        hFig.CurrentAxes.GridAlpha=0.5;
        
        % now style the Violins
        for ii=1:Ngrps
            if length(cfg.global_vis_colors)>=ii
                hV(ii).ViolinColor=cfg.global_vis_colors{ii};
            end
        end
        
        if cfg.global_vis_log
            hFig.CurrentAxes.MinorGridColor='w';
            hFig.CurrentAxes.MinorGridAlpha=0.8;
            % mark clipped points
            for ii=1:Ngrps
                if sum(mean_vals{i}(SubjFilt==ii)==logMin)>0
                    oldcol=repmat(hV(ii).ScatterPlot.MarkerFaceColor,[sum(SubjFilt==ii) 1]);
                    oldcol(mean_vals{i}(SubjFilt==ii)==logMin,:)=repmat([0 0 0],[sum(mean_vals{i}(SubjFilt==ii)==logMin) 1]); % make them black
                    hV(ii).ScatterPlot.MarkerFaceColor='flat';
                    hV(ii).ScatterPlot.CData=oldcol;
                end
            end
        end
        if cfg.global_vis_lock_scale
            ylim([lowTh highTh])
        end
        if cfg.global_vis_log
            hFig.CurrentAxes.YScale='log';
        end
        
        for gi=1:length(plot_contrasts)
            
            % check and print significance
            palmres_inc=fullfile(palmdir,['palm_out_dat_tstat_' mcorr 'fwep_m' num2str(i) '_c' num2str(plot_contrasts{gi}(2)) '.csv'] );
            sig_inc=csvread(palmres_inc); % inc in grp 2
            palmres_dec=fullfile(palmdir,['palm_out_dat_tstat_' mcorr 'fwep_m' num2str(i) '_c' num2str(plot_contrasts{gi}(1)) '.csv'] );
            sig_dec=csvread(palmres_dec); % decr in grp 2
            
            % take the bigger
            direction=0;
            if sig_inc>sig_dec
                sig=sig_inc;
                palmres_cohen=fullfile(palmdir,['palm_out_dat_cohen_m' num2str(i) '_c' num2str(plot_contrasts{gi}(2)) '.csv'] );
                direction=-1;
            else
                sig=sig_dec;
                palmres_cohen=fullfile(palmdir,['palm_out_dat_cohen_m' num2str(i) '_c' num2str(plot_contrasts{gi}(1)) '.csv'] );
                
                direction=1;
            end
            cohen=csvread(palmres_cohen);
            
            % we may need to flip
            tmp=con_mat(plot_contrasts{gi}(1),:);
            tmp(tmp==0)=[];
            if tmp(1)==-1
                direction=-direction;
            end
            
            %if direction==1
            % txt='>';
            %elseif direction==-1
            % txt='<';
            %else
            txt='';
            %end
            
            p=power(10,-sig);
            if p<0.001
                txt=['\bf' txt '\rm p<0.001, d=' num2str(abs(cohen),'%.2f')];
            elseif p<0.05
                txt=['\bf' txt '\rm p=' num2str(p,'%.3f') ', d=' num2str(abs(cohen),'%.2f')];
            elseif p<0.1
                txt=[txt ' (p=' num2str(p,'%.3f') ', d=' num2str(abs(cohen),'%.2f') ')'];
            else
                txt='';
            end
            if ~isempty(txt)
                sigstar({pot(gi,:)},{txt})
            end
        end
        
        % add a legend in the last
        if i==length(cfg.data)
            legend([hV(:).ViolinPlot],cfg.design_col_labels(plot_groups),'Color','w','FontSize',12,'Location','best');
        end
        
        if length(cfg.global_viz_remove)>0
            this_plot=fullfile(palmdir,['Global_fig_Violin_confounds_removed_m' num2str(i)  '.png']);
        else
            this_plot=fullfile(palmdir,['Global_fig_Violin_m' num2str(i)  '.png']);
        end
        export_fig(hFig,this_plot,'-nocrop','-r200')
        close(hFig)
        
        % add freq title with same font
        cmd=['convert ' this_plot ' -shave 0x64 -font FreeSans-Bold -pointsize 48 -fill black -background "' cfg.viz_data_back{i} '" -gravity North -splice 0x64 -annotate +0+2 "' cfg.data_labels{i} '" ' this_plot];
        system(cmd);
        
        mcmd=[mcmd ' ' this_plot ' '];
    end
    
    txt='';
    if length(cfg.global_viz_remove)>0
        txt='removed confounding factors: ';
        for ci=1:length(cfg.global_viz_remove)
            txt=[txt strrep(cfg.global_viz_remove{ci},'grp_','')];
            if ci<length(cfg.global_viz_remove)
                txt=[txt ', '];
            end
        end
        opf=fullfile(cfg.output,'Global_fig_Violin_confounds_removed.png');
    else
        opf=fullfile(cfg.output,'Global_fig_Violin.png');
    end
    system([mcmd opf]);
    
    if ~isempty(txt)
        cmd=['convert ' opf ' -font FreeSans -pointsize 48 -fill black -background white -gravity SouthEast -annotate +6+20 "' txt '" ' opf];
        system(cmd);
    end
end
end