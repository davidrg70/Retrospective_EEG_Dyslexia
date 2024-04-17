# Retrospective_EEG_Dyslexia

** Scripts with "nmri" are modified but part of the pipeline for EEG processing, analysis, source reconstruction described by Marquetand et al (2019 - https://doi.org/10.1089/brain.2019.0662). The scripts available here were used for the retrospective Dyslexia study (Preprint: https://doi.org/10.21203/rs.3.rs-2895521/v1). Additionally, several scripts were used to bring the source reconstructed data to permutational analysis in PALM.

** Scripts with "dg" initials were additionally programmed by David Garnica for EEG data management, processing, analysis, and results visualization.

  ** Script "dg_fft_power_spectra.m" This script performs Power Spectrum Density (PSD) analysis of one or two groups. It implements usual   MATLAB functions to calculate relative and after absolute power. Finally, plots the absolute power according to the frequency bands       analyzed, and smoothes to the average line and the 95% confidence interval area.

  ** Script "dg_Dys_cognitiveData_SignVertices_averaged_wholeGroup.m" and "dg_Dys_cognitiveData_SignRegions_averaged_wholeGroup.m" These    scripts collect vertex-based (source-reconstructed) power or functional connectivity values and investigate a correlation between them    and Z-scores of cognitive tasks. They implement Pearson, Spearman, or Partial rank correlations, after the user's selection and           normality test (Shapiro-Wilk). Additionally, they apply FDR correction (https://www.mathworks.com/matlabcentral/fileexchange/27418-       fdr_bh) for multiple correlations. As controls did not have psychological testing in my study, the scripts create a mask to investigate   correlations exclusively with power or functional connectivity values of vertices that showed statistically significant differences       between patients and controls. The first script averages all vertices, but the second averages power or FC values of every region in      the Desikan-Kiliany atlas (Desikan et al., 2006 - https://doi.org/10.1016/j.neuroimage.2006.01.021).

  ** Script "dg_violinplots_globalData.m" This scripts plots power or functional connectivity values for 2 groups previously compared,      using violin plots. The script implements a violin plot function (https://github.com/bastibe/Violinplot-Matlab)

** Files called "SourceLevel_Results_96children_1st-4th_graders" and "SourceLevel_Results_24children_5th-8th_graders" are results after source reconstruction (3 metrics: Power ImCoh, and wPLI debiased) but obtained BEFORE permutational analysis.
