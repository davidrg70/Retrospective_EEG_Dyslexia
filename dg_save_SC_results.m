% dg_save_SC_results. DG April 2024

% This script saves the source-reconstructed results after (IN) export but
% before permutational analysis with PALM. It saves results of the 3
% studied metrics, and from processing+1st level analysis of children
% (Dys+Controls) at 1st-4th and 5th-8th school grades

clear all; clc;

analysis_dir = 'main dir...';
folders = {'export_12_1to4_pediatric_template','export_13_5to8_pediatric_template'};

metrics = {'coh_img','power','wpli_debiased'};
fbands = {'Delta','Theta','Alpha','Beta1','Beta2'};

%% loop over...
for f = 1:length(folders)
    if f == 1
        N = 96;
        Results_1to4 = struct;
    elseif f == 2
        N = 24;
        Results_5to8 = struct;
    end
    
    for m = 1:length(metrics)
        clear file metrics_cells scale_metrics;
        for ff = 1:length(fbands)
            clear scale_metrics;
            file = fullfile(analysis_dir, folders{f}, [metrics{m}, '_', fbands{ff}, ...
                '_abs_not_scaled_ANTS9Years3T_fs7_openmeeg_N', num2str(N), '.mat' ]);
            load(file);
            metrics_cells{:,ff} = scale_metrics;
        end
        metrics_cells = cat(1, fbands, metrics_cells);
        
        % save collection here into struct fields
        if f == 1
            Results_1to4.(metrics{m}) = metrics_cells;
        elseif f == 2
            Results_5to8.(metrics{m}) = metrics_cells;
        end
    end
end

save('SourceLevel_Results_96children_1st-4th_graders.mat', 'Results_1to4');
save('SourceLevel_Results_24children_5th-8th_graders.mat', 'Results_5to8');