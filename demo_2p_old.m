%% Recording
% 1. find bloodvessel pattern
alignment
% 2. find picked planes
sbxshownplane
findmatchedplane
% 2. record PMT0+1
% 3. record running alone (no visual stim)
% 4. record the training stimulus
% 5. blood vessel pattern

%% Analysis
%1. use memmapdir for automatic cell extraction
memmapdir(4) % memmapdir(2);
% 2a. to sync with stimulus
caanalysisdir;

%% for binocular studies
calculate_ODI4; % 2c. to calculate ODI
% analysis_twoeye_data
align_two_days %compare_two_days
plot_twoday_data

%% for monocular studies
align_images;% 2b. to align two day images
analysis_twoday_data % 3. to match cells and analyze
calculate_MI
process_seperate_days_ori; % 4. use matched cell to analyze other data
MLvStimcorrelation

% for z-stack alignment and saveastif
sbxalignstackdir
