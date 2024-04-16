% This script will load exported connectivity/power values at each vertex
% and average for each of the yeo network and store as .mat and -mgh files
% in the export folder (requires cs_nf_ft_report regional, which, by now, should be the same as nmri_report_regional)
% C. Stier 2022

% Modified by D. Garnica 10.2022

%% AVERAGE FC & POWER VALUES PER EVERY RSN YEO NETWORK
load('atlases directory.../Yeo2011_7Networks_N1000_suma-all-fsaverage-10.mat','suma_all') % change to your folder

metric = {'coh_img_', 'wpli_debiased_', 'power_'};  % change to your needs
% freq = {'Delta_','Theta_','Alpha_','Beta1_','Beta2_','Gamma_'}; % change to your needs
freq = {'Delta_','Theta_','Alpha_','Beta1_','Beta2_'};
scale = {'abs_not_scaled_'};
rootdir = ('export directory'); % change to your export folder
number_subjects = '_N120';      % change accordingly
space = 'nmriSUMA150_fs6_openmeeg';      % change accordingly

group = 'all';
yeo_all = {};

for var_freq = 1:length(freq)
    for var_metric = 1:length(metric)
        data = load([rootdir metric{var_metric} freq{var_freq} scale{1} space number_subjects '.mat']); % adapt here
        for i = 1:length(data.scale_metrics) % change here...
            m3 = data.scale_metrics{i};
            
            m3(2005:2338,1) = NaN; % I set subcortical nuclei to NaN and thus, these numbers/vertices won't be included in the analysis.
            regvalues_th = cs_nf_ft_report_regional(suma_all,m3);
            regvalues_th(1,:) = [];
            roi = regvalues_th(:,1);
            regvalues_th(:,1:3) = [];
            yeo_all{i,1} = cell2mat(regvalues_th);
        end
        
        filebase=fullfile(rootdir, [ metric{var_metric} freq{var_freq} scale{1} 'yeo_N' num2str(length(yeo_all))]);
        all_subjects = data.all_subjects;
        save([filebase '.mat'],'yeo_all','all_subjects')
        
        nmri_write_mgh([filebase '.mgh'],eye(4),yeo_all)
    end
end

%% RUN PALM FOR YEO PARCELLATED FC AND POWER
% Example code for PALM on yeo-networks
% change folders to your needs
analysis_dir = 'main dir';
results_dir = 'results dir'; % change

% Christina's example
% palm -i export/all_subjects/coh_img_Alpha_abs_not_scaled_yeo_N350.mgh -m mask_yeo7.csv -s 1 -d design_all_corr2_sex_manualdemean.csv -demean -t con_corr2_quadr2_all_sex.csv -logp -saveglm -savedof -accel tail -n 500 -o results/corr2_quadr2_sex_manualdemean_plusdemean_desikan/Alpha/coh_img

% CURRENT, TO USE:
% NOTE: a mask [0; 1; 1; 1; 1; 1; 1; 1; 0; 1; 1; 1; 1; 1; 1; 1] is included to ignore the Unknown and Freesurfer Medial Wall networks
load('main dir.../export_10/all_suma_msk.mat')
all_msk(17:2338) = [];
all_msk(9,1) = 0;
all_msk(12:13,1) = 1;
save_mgh(all_msk,'mask_yeo.mgh', eye(4)) % It is saved in the main directory, but better move it to the corresponding EXPORT folder

% ImCoh
palm -i /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/coh_img_Alpha_abs_not_scaled_yeo_N120.mgh -m /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/mask_yeo.mgh -s 1 -d /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/design.csv -demean -t /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/contrast.csv -logp -saveglm -savedof -accel tail -n 1000 -o /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/Dys_noIEDs_Controls/1_RSNs_dyslexia_controls/ImCoh/Alpha/
palm -i /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/coh_img_Beta1_abs_not_scaled_yeo_N120.mgh -m /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/mask_yeo.mgh -s 1 -d /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/design.csv -demean -t /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/contrast.csv -logp -saveglm -savedof -accel tail -n 1000 -o /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/Dys_noIEDs_Controls/1_RSNs_dyslexia_controls/ImCoh/Beta1/
palm -i /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/coh_img_Beta2_abs_not_scaled_yeo_N120.mgh -m /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/mask_yeo.mgh -s 1 -d /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/design.csv -demean -t /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/contrast.csv -logp -saveglm -savedof -accel tail -n 1000 -o /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/Dys_noIEDs_Controls/1_RSNs_dyslexia_controls/ImCoh/Beta2/
palm -i /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/coh_img_Delta_abs_not_scaled_yeo_N120.mgh -m /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/mask_yeo.mgh -s 1 -d /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/design.csv -demean -t /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/contrast.csv -logp -saveglm -savedof -accel tail -n 1000 -o /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/Dys_noIEDs_Controls/1_RSNs_dyslexia_controls/ImCoh/Delta/
palm -i /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/coh_img_Theta_abs_not_scaled_yeo_N120.mgh -m /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/mask_yeo.mgh -s 1 -d /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/design.csv -demean -t /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/contrast.csv -logp -saveglm -savedof -accel tail -n 1000 -o /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/Dys_noIEDs_Controls/1_RSNs_dyslexia_controls/ImCoh/Theta/

% wPLI
palm -i /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/wpli_debiased_Alpha_abs_not_scaled_yeo_N120.mgh -m /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/mask_yeo.mgh -s 1 -d /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/design.csv -demean -t /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/contrast.csv -logp -saveglm -savedof -accel tail -n 1000 -o /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/Dys_noIEDs_Controls/1_RSNs_dyslexia_controls/wPLI/Alpha/
palm -i /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/wpli_debiased_Beta1_abs_not_scaled_yeo_N120.mgh -m /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/mask_yeo.mgh -s 1 -d /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/design.csv -demean -t /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/contrast.csv -logp -saveglm -savedof -accel tail -n 1000 -o /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/Dys_noIEDs_Controls/1_RSNs_dyslexia_controls/wPLI/Beta1/
palm -i /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/wpli_debiased_Beta2_abs_not_scaled_yeo_N120.mgh -m /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/mask_yeo.mgh -s 1 -d /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/design.csv -demean -t /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/contrast.csv -logp -saveglm -savedof -accel tail -n 1000 -o /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/Dys_noIEDs_Controls/1_RSNs_dyslexia_controls/wPLI/Beta2/
palm -i /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/wpli_debiased_Delta_abs_not_scaled_yeo_N120.mgh -m /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/mask_yeo.mgh -s 1 -d /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/design.csv -demean -t /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/contrast.csv -logp -saveglm -savedof -accel tail -n 1000 -o /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/Dys_noIEDs_Controls/1_RSNs_dyslexia_controls/wPLI/Delta/
palm -i /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/wpli_debiased_Theta_abs_not_scaled_yeo_N120.mgh -m /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/mask_yeo.mgh -s 1 -d /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/design.csv -demean -t /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/contrast.csv -logp -saveglm -savedof -accel tail -n 1000 -o /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/Dys_noIEDs_Controls/1_RSNs_dyslexia_controls/wPLI/Theta/

% Power
palm -i /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/power_Alpha_abs_not_scaled_yeo_N120.mgh -m /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/mask_yeo.mgh -s 1 -d /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/design.csv -demean -t /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/contrast.csv -logp -saveglm -savedof -accel tail -n 1000 -o /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/Dys_noIEDs_Controls/1_RSNs_dyslexia_controls/Power/Alpha/
palm -i /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/power_Beta1_abs_not_scaled_yeo_N120.mgh -m /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/mask_yeo.mgh -s 1 -d /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/design.csv -demean -t /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/contrast.csv -logp -saveglm -savedof -accel tail -n 1000 -o /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/Dys_noIEDs_Controls/1_RSNs_dyslexia_controls/Power/Beta1/
palm -i /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/power_Beta2_abs_not_scaled_yeo_N120.mgh -m /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/mask_yeo.mgh -s 1 -d /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/design.csv -demean -t /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/contrast.csv -logp -saveglm -savedof -accel tail -n 1000 -o /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/Dys_noIEDs_Controls/1_RSNs_dyslexia_controls/Power/Beta2/
palm -i /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/power_Delta_abs_not_scaled_yeo_N120.mgh -m /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/mask_yeo.mgh -s 1 -d /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/design.csv -demean -t /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/contrast.csv -logp -saveglm -savedof -accel tail -n 1000 -o /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/Dys_noIEDs_Controls/1_RSNs_dyslexia_controls/Power/Delta/
palm -i /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/power_Theta_abs_not_scaled_yeo_N120.mgh -m /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/mask_yeo.mgh -s 1 -d /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/design.csv -demean -t /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/export_10/contrast.csv -logp -saveglm -savedof -accel tail -n 1000 -o /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/Dys_noIEDs_Controls/1_RSNs_dyslexia_controls/Power/Theta/

% metric = {'coh_img_', 'power_'};
% freq = {'Delta','Theta','Alpha','Beta1','Beta2','Gamma'};
% scale = {'abs_not_scaled_yeo_'};
%
% for var_freq = 1:length(freq)
%  for var_metric = 1:length(metric)
%   for i = 1:16
%    func_file = fullfile(export_func, ['7networks_yeo_N350_', metric{var_metric}, freq{var_freq}, '_', hemi{i}, '_', nw{i}, '.csv']);
%    design_file = fullfile(export_th, ['yeo7_N350_th_age_age2_sex_icv_', hemi{i}, '_', nw{i}, '.csv']);
%
%    results = fullfile(results_dir, freq{var_freq}, [metric{var_metric}, hemi{i}, '_', nw{i}]);
%
%    code_line = [{'-i', func_file, '-d', design_file, '-demean', '-t', 'con_yeo_th_age_age2_icv_sex.csv', '-saveglm', '-savedof', '-accel', 'tail', '-n', '500', '-o', results}];
%
%    palm(code_line{:})
%   end
%  end
% end
%
%   palm -i func_file -d design_file -demean -t con_yeo_th_age_age2_icv_sex.csv -saveglm -savedof -accel tail -n 500 -o results

%% Original CS visualization code
% Example code for visualization of yeo-results: example for one frequency band and contrast

% load([pwd 'atlases dir.../Yeo2011_7Networks_N1000_suma-all-fsaverage-10.mat'],'suma_all')
%
% opt=[];
% opt.per_hemi=1;
% opt.per_cortex=1;
% opt.rot=[90 0 ; -90 0];
% opt.thresh=1.3;
% opt.clim=[1.3 3.8];
% opt.colormap='hot'; % change color bars if needed
% opt.colorbar='hot';
% opt.scale=1;
%
% %%% Alpha C1
% freqname = {'Alpha_'};
% metric = {'_'}; % change how you like it
% analysis_type = {'dpx_'};
% analysis_type_2={'tstat_fwep_', 'tstat_uncp_'};% change to your needs
% contrasts = {'c1'}; % change to your needs
%
% rootdir = [results_dir '/ImCoh/Alpha/'];
% fig_dir = [rootdir 'fig_Alpha_'];
% if(~exist(fig_dir,'dir'))
%     mkdir(fig_dir)
% end
%
% for var_metric = 1:length(metric)
%     for var_analysis_type = 1:length(analysis_type)
%         for var_analysis_type_2 = 1:length(analysis_type_2)
%             for var_contrasts = 1:length(contrasts)
%
%                 logp=load_mgh([rootdir...
%                     metric{1, var_metric}...
%                     analysis_type{1,var_analysis_type} analysis_type_2{1, var_analysis_type_2}...
%                     contrasts{1, var_contrasts} '.mgz']);
%
%                    keydata = [logp (suma_all.annot_key{1,1}(2:17))]; % get rid of unknown index as subcortical nuclei have been set to NaN before and are not included
%                    vertex_data = zeros(2338,2);
%
%                    for v = 1:length(suma_all.annot)
%                        keynr = suma_all.annot(v);
%                        if keynr == 0
%                            vertex_data(v,1) = 0;
%                        else
%                            rownr = find(keydata(:,2) == keynr);
%                            log_p_v = keydata(rownr,1);
%                            vertex_data(v,1) = log_p_v;
%                            logp = vertex_data(:,1);
%                        end
%                    end
%
%                 opt.title=strcat(cell2mat(freqname),cell2mat(metric(1, var_metric)),cell2mat(analysis_type(1,var_analysis_type)),cell2mat(analysis_type_2(1, var_analysis_type_2)),cell2mat(contrasts(1, var_contrasts))) ;
%                 opt.output=strcat(fig_dir,cell2mat(metric(1, var_metric)),cell2mat(analysis_type(1,var_analysis_type)),cell2mat(analysis_type_2(1, var_analysis_type_2)),cell2mat(contrasts(1, var_contrasts)),'.png');
%
%                 hFig = nmri_plot_surface_suma(suma_all, logp, opt);
%             end
%         end
%     end
% end

%% Code for visualization of Yeo-results
cd /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls
load([pwd 'atlases dir.../Yeo2011_7Networks_N1000_suma-all-fsaverage-10.mat'],'suma_all')
% load('all_suma_msk_zerosubc.mat');

analysis_dir = 'main dir';
results_dir = 'results dir'; % change

freqname = {'Alpha', 'Beta1', 'Beta2', 'Delta', 'Theta'};
metric = {'coh_img', 'wpli_debiased', 'power'}; % change to your needs
analysis_type = {'dpx_'};
analysis_type_2={'tstat_fwep_', 'tstat_uncp_'};% change to your needs
contrasts = {'c1', 'c2', 'c3', 'c4', 'c5', 'c6'}; % change to your needs

for i = 1:length(freqname)
    for var_metric = 1:length(metric)
        
        if var_metric == 1
            folder = '/ImCoh/';
        elseif var_metric == 2
            folder = '/wPLI/';
        elseif var_metric == 3
            folder = '/Power/';
        end
        
        rootdir = [results_dir folder freqname{i} '/'];
        fig_dir = [rootdir 'fig_' freqname{i}];
        if(~exist(fig_dir,'dir'))
            mkdir(fig_dir)
        end
        
        for var_analysis_type = 1:length(analysis_type)
            for var_analysis_type_2 = 1:length(analysis_type_2)
                for var_contrasts = 1:length(contrasts)
                    
                    opt=[];
                    opt.per_hemi=1;
                    opt.per_cortex=1;
                    opt.rot=[90 0 ; -90 0];
                    opt.thresh=1.3;
                    opt.clim=[1.3 3.8];
                    opt.scale=1;
                    if var_contrasts == 1
                        opt.colormap='cool'; opt.colorbar='cool';
                    elseif var_contrasts == 2
                        opt.colormap='autumn'; opt.colorbar='autumn';
                    elseif var_contrasts == 3 || var_contrasts == 4 || var_contrasts == 5 || var_contrasts == 6
                        opt.colormap='parula'; opt.colorbar='parula';
                    end
                    
                    % logp = load_mgh([rootdir metric{1, var_metric} analysis_type{1,var_analysis_type} analysis_type_2{1, var_analysis_type_2} contrasts{1, var_contrasts} '.mgz']);
                    % logp = table2array(readtable([rootdir '_' analysis_type{1,var_analysis_type} analysis_type_2{1, var_analysis_type_2} contrasts{1, var_contrasts} '.csv']));
                    logp = load_mgh([rootdir '_' analysis_type{1,var_analysis_type} analysis_type_2{1, var_analysis_type_2} contrasts{1, var_contrasts} '.mgz']);
                    logp(1,1) = 0;
                    logp(9,1) = 0;
                    keydata = [logp (suma_all.annot_key{1,1}(2:17))]; % get rid of unknown index as subcortical nuclei have been set to NaN before and are not included
                    vertex_data = zeros(2338,2);
                    
                    for v = 1:length(suma_all.annot)
                        keynr = suma_all.annot(v);
                        if keynr == 0
                            vertex_data(v,1) = 0;
                        else
                            rownr = find(keydata(:,2) == keynr);
                            log_p_v = keydata(rownr,1);
                            vertex_data(v,1) = log_p_v;
                            logp = vertex_data(:,1);
                        end
                    end
                    
                    opt.title=strcat(cell2mat(freqname(i)),'_',cell2mat(metric(1, var_metric)),'_',cell2mat(analysis_type(1,var_analysis_type)),cell2mat(analysis_type_2(1, var_analysis_type_2)),cell2mat(contrasts(1, var_contrasts)));
                    opt.output=strcat(fig_dir,'/',cell2mat(metric(1, var_metric)),'_',cell2mat(analysis_type(1,var_analysis_type)),cell2mat(analysis_type_2(1, var_analysis_type_2)),cell2mat(contrasts(1, var_contrasts)),'.png');
                    
                    hFig = cs_nmri_plot_surface_suma_yeo(suma_all, logp, opt);
                end
            end
        end
    end
end

%% Concatenate the multiple images!
% you can also use this function, and the whole thing will be done for several contrasts (in this case 3) and
% frequency bands. In the end, you may concatenate results for each frequency band if needed if you do not use
% gamma, you have to exclude it in the script below. Also, change the amount of contrasts if needed ...

% cd /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls
% results_dir = '/home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/Dys_noIEDs_Controls/1_RSNs_dyslexia_controls'; % change
% load([pwd '/conf/atlas/Yeo2011_7Networks_N1000_suma-all-fsaverage-10.mat'],'suma_all')
% cs_visualize_single_palm_p_3contrasts_yeo_david(results_dir)
%
% % concatenate
% copyfile ([results_dir '/Alpha/fig_*'], [results_dir '/figures']);
% copyfile ([results_dir '/Beta1/fig_*'], [results_dir '/figures']);
% copyfile ([results_dir '/Beta2/fig_*'], [results_dir '/figures']);
% copyfile ([results_dir '/Delta/fig_*'], [results_dir '/figures']);
% copyfile ([results_dir '/Theta/fig_*'], [results_dir '/figures']);
% % copyfile ([results_dir '/Gamma/fig_*'], [results_dir '/figures']);
%
% rootdir = [results_dir '/figures'];
% fig_dir = [rootdir '/' 'concat_fig'];
% if(~exist(fig_dir,'dir'))
%     mkdir(fig_dir)
% end
%
% % frequencies = {'fig_Delta_','fig_Theta_','fig_Alpha_','fig_Beta1_','fig_Beta2_','fig_Gamma_'};
% frequencies = {'fig_Delta_','fig_Theta_','fig_Alpha_','fig_Beta1_','fig_Beta2_'};
% metric = {'coh_img', 'wpli_debiased', 'power'}; % change to what you need
% analysis_type = {'_dpx_'}; % change to what you need
% analysis_type_2={'tstat_fwep_', 'tstat_uncp_'}; % change to what you need
% contrasts = {'c1', 'c2', 'c3', 'c4', 'c5', 'c6'}; % change to what you need
%
% for var_frequencies = 1:length(frequencies)
%     for var_metric = 1:length(metric)
%         for var_analysis_type = 1:length(analysis_type)
%             for var_analysis_type_2 = 1:length(analysis_type_2)
%                 for var_contrasts = 1:length(contrasts)
%
%                     name =[fig_dir, '/fig_all_freq_' ...
%                         metric{1, var_metric}...
%                         analysis_type{1,var_analysis_type} analysis_type_2{1, var_analysis_type_2}...
%                         contrasts{1, var_contrasts} '.png'];
%                     img1 = imread([rootdir...
%                         frequencies{1}...
%                         metric{1, var_metric}...
%                         analysis_type{1,var_analysis_type}...
%                         analysis_type_2{1, var_analysis_type_2}...
%                         contrasts{1, var_contrasts} '.png']);
%                     img2 = imread([rootdir...
%                         frequencies{2}...
%                         metric{1, var_metric}...
%                         analysis_type{1,var_analysis_type}...
%                         analysis_type_2{1, var_analysis_type_2}...
%                         contrasts{1, var_contrasts} '.png']);
%                     img3 = imread([rootdir...
%                         frequencies{3}...
%                         metric{1, var_metric}...
%                         analysis_type{1,var_analysis_type}...
%                         analysis_type_2{1, var_analysis_type_2}...
%                         contrasts{1, var_contrasts} '.png']);
%                     img4 = imread([rootdir...
%                         frequencies{4}...
%                         metric{1, var_metric}...
%                         analysis_type{1,var_analysis_type}...
%                         analysis_type_2{1, var_analysis_type_2}...
%                         contrasts{1, var_contrasts} '.png']);
%                     img5 = imread([rootdir...
%                         frequencies{5}...
%                         metric{1, var_metric}...
%                         analysis_type{1,var_analysis_type}...
%                         analysis_type_2{1, var_analysis_type_2}...
%                         contrasts{1, var_contrasts} '.png']);
%                     % img6 = imread([rootdir...
%                     %     frequencies{6}...
%                     %     metric{1, var_metric}...
%                     %     analysis_type{1,var_analysis_type}...
%                     %     analysis_type_2{1, var_analysis_type_2}...
%                     %     contrasts{1, var_contrasts} '.png']);
%
%                     %                 img = [img1; img2; img3; img4; img5; img6]; % change if less than 6 frequency bands
%                     img = [img1; img2; img3; img4; img5]; % change if less than 6 frequency bands
%                     imwrite (img, name);
%                 end
%             end
%         end
%     end
% end

% base:
% name = 'C1_ImCoh_tstat_fwep.png';
% im1_delta_str = imread(string(fullfile(results_dir, 'ImCoh', 'Delta', 'fig_Delta', 'coh_img_dpx_tstat_fwep_c1.png')));
% im1_theta_str = imread(string(fullfile(results_dir, 'ImCoh', 'Theta', 'fig_Theta', 'coh_img_dpx_tstat_fwep_c1.png')));
% im1_alpha_str = imread(string(fullfile(results_dir, 'ImCoh', 'Alpha', 'fig_Alpha', 'coh_img_dpx_tstat_fwep_c1.png')));
% im1_beta1_str = imread(string(fullfile(results_dir, 'ImCoh', 'Beta1', 'fig_Beta1', 'coh_img_dpx_tstat_fwep_c1.png')));
% im1_beta2_str = imread(string(fullfile(results_dir, 'ImCoh', 'Beta2', 'fig_Beta2', 'coh_img_dpx_tstat_fwep_c1.png')));
% img = [im1_delta_str; im1_theta_str; im1_alpha_str; im1_beta1_str; im1_beta2_str];
% imwrite (img, name);

cd /main dir...
results_dir = 'results dir'; % change
load([pwd 'atlases dir.../Yeo2011_7Networks_N1000_suma-all-fsaverage-10.mat'],'suma_all')

freqname = {'Delta', 'Theta','Alpha', 'Beta1', 'Beta2'};
metric = {'coh_img', 'wpli_debiased', 'power'}; % change to your needs
analysis_type = {'dpx_'};
analysis_type_2={'tstat_fwep_', 'tstat_uncp_'}; % change to your needs
% contrasts = {'c1', 'c2', 'c3', 'c4', 'c5', 'c6'}; % change to your needs
contrasts = {'control_more_dyslexia', 'controls_less_dyslexia', 'Age-Positive', 'Age-Negative', 'Sex-Positive', 'Sex-Negative'};

for var_metric = 1:length(metric)
    if var_metric == 1
        folder = 'ImCoh';
    elseif var_metric == 2
        folder = 'wPLI';
    elseif var_metric == 3
        folder = 'Power';
    end
    
    for var_analysis_type = 1:length(analysis_type)
        for var_analysis_type_2 = 1:length(analysis_type_2)
            for var_contrasts = 1:length(contrasts)
                if var_contrasts == 1
                    contrast = 'c1';
                elseif var_contrasts == 2
                    contrast = 'c2';
                elseif var_contrasts == 3
                    contrast = 'c3';
                elseif var_contrasts == 4
                    contrast = 'c4';
                elseif var_contrasts == 5
                    contrast = 'c5';
                elseif var_contrasts == 6
                    contrast = 'c6';
                end
                
                image_string = [];
                for i = 1:length(freqname)
                    if i == 1
                        image_subfolder = 'fig_Delta';
                    elseif i == 2
                        image_subfolder = 'fig_Theta';
                    elseif i == 3
                        image_subfolder = 'fig_Alpha';
                    elseif i == 4
                        image_subfolder = 'fig_Beta1';
                    elseif i == 5
                        image_subfolder = 'fig_Beta2';
                    end
                    
                    name = [folder, '_', 'Yeo7_', analysis_type{var_analysis_type}, analysis_type_2{var_analysis_type_2}, contrasts{var_contrasts}, '.png'];
                    image_string{i} = imread(string([results_dir, '/', folder, '/', freqname{i}, '/', image_subfolder, '/', metric{var_metric}, ...
                        '_', analysis_type{var_analysis_type}, analysis_type_2{var_analysis_type_2}, contrast, '.png']));
                end
                image_string = image_string';
                img = cell2mat([image_string(1); image_string(2); image_string(3); image_string(4); image_string(5)]);
                imwrite (img, [results_dir, '/', name]);
            end
        end
    end
end
