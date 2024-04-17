% EXAMPLE (wiki)
palm -i export/17-Feb-2017/coh_img_Alpha_abs_not_scaled_N8_all.mgh -m export/17-Feb-2017/all_suma_msk.mgh -s conf/suma-all-fsaverage-10.gii -d design.csv -t contrast.csv -T -tfce2d -mv Wilks -demean -o results/palm_test/my_group_analysis/palm

% EXAMPLE (Christina)
palm -i export/all_subjects/coh_img_Alpha_abs_not_scaled_desikan_N350.mgh -s 1 -d design_all_corr2_sex_manualdemean.csv -demean -t con_corr2_quadr2_all_sex.csv -logp -saveglm -savedof -accel tail -n 5000 -o results/corr2_quadr2_sex_manualdemean_plusdemean_desikan/Alpha/coh_img

%% MY PALM RUN (runs PALM for connectivity-cognition analysis ---> 1-GROUP!)

% MY GUIDE: palm -i mghFreqBandFile -m sumaMaskFile -s mghSurfaceFile -d design.csv -demean -t contrast.csv -logp -saveglm -savedof -T -tfce2d -corrmod -corrcon -accel tail -n 5000 -o folder
% BUT I REMOVED THE DEMEAN OPTION!

% Power


% ImCoh
% done
palm -i /data dir.../coh_img_Alpha_abs_not_scaled_nmriSUMA150_fs6_openmeeg_N68_all.mgh -m /data dir.../export_dys_WISC_IQs/all_suma_msk.mgh -s /data dir.../conf/suma-all-fsaverage-10.gii -d /data dir.../group_filters/design_IQ_FS.csv -demean -t /data dir.../group_filters/contrast_IQ_FS.csv -logp -saveglm -savedof -T -tfce2d -corrmod -corrcon -accel tail -n 5000 -o /data dir.../8_dys_FSIQ_22_08_2022
% not done
palm -i /data dir.../export_dys_WISC_IQs/coh_img_Alpha_abs_not_scaled_nmriSUMA150_fs6_openmeeg_N68_all.mgh -m /data dir.../export_dys_WISC_IQs/all_suma_msk.mgh -s /data dir.../conf/suma-all-fsaverage-10.gii -d /data dir.../group_filters/design_IQ_FS.csv -demean -t /data dir.../group_filters/contrast_IQ_FS.csv -logp -saveglm -savedof -T -tfce2d -corrmod -corrcon -accel tail -n 5000 -o /data dir.../8_dys_FSIQ_22_08_2022
palm -i /data dir.../export_dys_WISC_IQs/coh_img_Alpha_abs_not_scaled_nmriSUMA150_fs6_openmeeg_N68_all.mgh -m /data dir.../export_dys_WISC_IQs/all_suma_msk.mgh -s /data dir.../conf/suma-all-fsaverage-10.gii -d /data dir.../group_filters/design_IQ_FS.csv -demean -t /data dir.../group_filters/contrast_IQ_FS.csv -logp -saveglm -savedof -T -tfce2d -corrmod -corrcon -accel tail -n 5000 -o /data dir.../8_dys_FSIQ_22_08_2022
palm -i /data dir.../export_dys_WISC_IQs/coh_img_Alpha_abs_not_scaled_nmriSUMA150_fs6_openmeeg_N68_all.mgh -m /data dir.../export_dys_WISC_IQs/all_suma_msk.mgh -s /data dir.../conf/suma-all-fsaverage-10.gii -d /data dir.../group_filters/design_IQ_FS.csv -demean -t /data dir.../group_filters/contrast_IQ_FS.csv -logp -saveglm -savedof -T -tfce2d -corrmod -corrcon -accel tail -n 5000 -o /data dir.../8_dys_FSIQ_22_08_2022
palm -i /data dir.../export_dys_WISC_IQs/coh_img_Alpha_abs_not_scaled_nmriSUMA150_fs6_openmeeg_N68_all.mgh -m /data dir.../export_dys_WISC_IQs/all_suma_msk.mgh -s /data dir.../conf/suma-all-fsaverage-10.gii -d /data dir.../group_filters/design_IQ_FS.csv -demean -t /data dir.../group_filters/contrast_IQ_FS.csv -logp -saveglm -savedof -T -tfce2d -corrmod -corrcon -accel tail -n 5000 -o /data dir.../8_dys_FSIQ_22_08_2022

% wPLI


% to open the mgz files later:
% uiopen('/data dir.../coh_img_Alpha/8_dys_FSIQ_22_08_2022_tfce_tstat_fwep_c1.mgz',1);
% load_mgh('data dir.../8_dys_FSIQ_22_08_2022/coh_img_Alpha/8_dys_FSIQ_22_08_2022_tfce_tstat_fwep_c1.mgz');

%% 