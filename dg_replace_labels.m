function [long_labels_hemi, long_labels_list] = dg_replace_labels(short_labels_list)
% DavidG, 05.2023
% This function compares every short label and replaces it with its corresponding longer label.
% Short labels after Desikan-Killiany atlas (https://doi.org/10.1016/j.neuroimage.2006.01.021),
% and longer labels taken from Destrieux (https://doi.org/10.1016/j.neuroimage.2010.06.010)

load('Desikan_Destrieux_labels.mat');
new_label = {};

for i = 1:size(short_labels_list,1)
    tf_lh = contains(short_labels_list{i}, 'wm_lh_');
    tf_rh = contains(short_labels_list{i}, 'wm_rh_');
    
    % in case there are no "wm_lh_/wm_rh_" initials
    temp = lower(short_labels_list{i}); % convert the line to lowercase (in case there are uppercases)
    tf_lh2 = contains(temp, 'left');
    tf_rh2 = contains(temp, 'right');
    
    if tf_lh == 1 && tf_rh == 0
        lt = erase(short_labels_list{i}, 'wm_lh_');
    elseif tf_lh == 0 && tf_rh == 1
        lt = erase(short_labels_list{i}, 'wm_rh_');
    elseif tf_lh2 == 1 && tf_rh2 == 0
        lt = erase(temp, 'left ');
        % check if gyrus, sulcus, or both, and Replace!
        if contains(lt, 'gyrus and sulcus')
            temp = regexprep(lt, ' ', '_'); % replace spaces with underscore
            lt = regexprep(temp, 'gyrus_and_sulcus', 'G_and_S'); % replace with DK labelling
        elseif contains(lt, 'gyrus') && contains(lt, 'sulcus')
            temp = regexprep(lt, ' ', '_');
            lt = regexprep(temp, 'gyrus_', 'G_');
            lt = regexprep(lt, 'sulcus_', 'S_');
        elseif contains(lt, 'gyrus')
            temp = regexprep(lt, ' ', '_');
            lt = regexprep(temp, 'gyrus_', 'G_');
        elseif contains(lt, 'sulcus')
            temp = regexprep(lt, ' ', '_');
            lt = regexprep(temp, 'sulcus_', 'S_');
        elseif contains(lt, 'lat')
            temp = regexprep(lt, ' ', '_');
            lt = regexprep(temp, 'lat_', 'Lat_');
        elseif contains(lt, 'pole')
            temp = regexprep(lt, ' ', '_');
            lt = regexprep(temp, 'pole_', 'Pole_');
        end            
    elseif tf_lh2 == 0 && tf_rh2 == 1
        lt = erase(temp, 'right ');
                % check if gyrus, sulcus, or both, and Replace!
        if contains(lt, 'gyrus and sulcus')
            temp = regexprep(lt, ' ', '_'); % replace spaces with underscore
            lt = regexprep(temp, 'gyrus_and_sulcus', 'G_and_S'); % replace with DK labelling
        elseif contains(lt, 'gyrus') && contains(lt, 'sulcus')
            temp = regexprep(lt, ' ', '_');
            lt = regexprep(temp, 'gyrus_', 'G_');
            lt = regexprep(lt, 'sulcus_', 'S_');
        elseif contains(lt, 'gyrus')
            temp = regexprep(lt, ' ', '_');
            lt = regexprep(temp, 'gyrus_', 'G_');
        elseif contains(lt, 'sulcus')
            temp = regexprep(lt, ' ', '_');
            lt = regexprep(temp, 'sulcus_', 'S_');
        elseif contains(lt, 'lat')
            temp = regexprep(lt, ' ', '_');
            lt = regexprep(temp, 'lat_', 'Lat_');
        elseif contains(lt, 'pole')
            temp = regexprep(lt, ' ', '_');
            lt = regexprep(temp, 'pole_', 'Pole_');
        end
    end
    
    % tl = find(contains(Desikan_Destrieux_labels(:,2),lt)); % finds corresponding Destrieux label
    tl = strcmpi(Desikan_Destrieux_labels(:,2),lt); % % finds corresponding Destrieux label (IGNORING UPPERCASE)
    new_label{i} = Desikan_Destrieux_labels(tl,3); % replaces DK to Destrieux
    
    % check hemisphere of every label, and saves it
    if tf_lh == 1 || tf_lh2 == 1
        hem{i} = 'Left';
    elseif tf_rh == 1 || tf_rh2 == 1
        hem{i} = 'Right';
    end
end

long_labels_list = lower(string(new_label'));
long_labels_hemi = string(hem');
end