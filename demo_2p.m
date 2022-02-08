
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
%% add all the programs
run('/nas/data1/jsun/MATLAB/toaddpath.m')
%% 1. cell signal extraction
% suite2p

sbx2tifdir
% old way:use memmapdir for automatic cell extraction
% memmapdir(4) % memmapdir(2);
%% pick good cells and save to Fall.mat
% in terminal run python
conda activate suite2p
python -m suite2p
%% convert Fall.mat and put it together with .sbx files
printcells2
%% 1a. suite2P
% sbx2tif
% make_db
% % /nas/data1/mdadarla/CA2p/Suite2P-master/configFiles/make_db.m
% master_file
% % /nas/data1/mdadarla/CA2p/Suite2P-master/master_file.m
% printcells
% % /nas/data1/mdadarla/CA2p/Suite2pData/J1/printcell.m

%% 2. to sync with stimulus
caanalysisdir;
%caanalysis1;
caanalysis;
caanalysisplot
%% 3. analysis for binocular studies
calculate_ODI5;
        plot_two_eyes
        cal_SNR % calculate Signal to noise ratio and plot with velocity and stim
        plotODIanalysis
        %%% some criteria
        ManuelCellSelection;
        AutoCellSelection; %%% 
        autoselect
        autoselect2
        %calculate_ODI4; % 2c. to calculate ODI
        plotODIanalysis
        plotODIcomparison
        plotTC%plot tuning curve for each neuron
% analysis_twoday / two eye
align_two_days % ==align_images
    compare_two_days
    plot_twoday_data
    OKRplot_single  %OKR
    OKRplot %\\mps-zfs\data1\jsun\projects\Plasticity\optometry
    celldensityplot %imaging

fftalignimg % plot two images, return [u v];
allplanes % combine information for all planes
plotODIcomparison; % process paired cells for all planes in different days
finalODIcomparison;
catstr % concatenate two structures with same name fields
SNR_comparison
cal_SNR
casignaltrend % 

ExcelexportBen % this have the rumsum criteria for good cells
%% for monocular studies

% 1. run suite2p
% 2. generate figs using
printcells2    % generate cell pictures from Fall.mat
% 3. plane alignment and cell matching
align_two_session    

% 4. analyze calcium responses    
caanalysis
    

plotrunrep % plot repetitions of running trials
align_images;% 2b. to align two day images  ==align_two_days
analysis_twoday_data
plot_twoday_data
plotpeakselected_env3
plotpeakselected_env2
plotpeakselected_env

%\\mps-zfs\jsun\2pdata\Enhancement
summary_analysis
summary_analysis3_permouse

plotposition
plotrunrep
SNR_comparison

% 3. to match cells and analyze
calculate_MI
process_seperate_days_ori; % 4. use matched cell to analyze other data
MLvStimcorrelation %running and still data
cacorrelation

% for z-stack alignment and saveastif
sbxalignstackdir
tifalign
% for rename 
renamefiles

summary_analysis2
%% ICMS study

% Plane alignment and cell matching
    printcells2
    % generate cell pictures from Fall.mat
    align_two_session 
    % align planes between 2 sessions and match cells
    align_two_session_find_duplicates_aashish
    % auxillary file to aid cell matching 
    plot_calcium_minus_neuropil_aashish
    % plot calcium - (coeff * neuropil) for different coeffs 
    match_selected_cells_aashish
    % match cells across all sessions,create matched cell index for post session

% File creation and concatenation    
    update_Fall
    % add structures to Fall.mat for each plane of a mouse from another
    % set of suite2p recordings
    create_new_files_stim
    % convert Fall.mat to F_stim_*.mat for each plane
    create_new_files_pre_post 
    % convert Fall.mat to pre and post F_*.mat for each plane. add other data
    concatenate_planes_pre_post
    % concatenate all 3 planes of pre and post matfiles to 1 matfile F_*_allplanes.mat
    concatenate_planes_stim
    % concatenate all 3 planes of stim to 1 matfile F_stim*_allplanes.mat

% Analysis     
    pre_post_analysis
    % wrapper that calls all the following functions
        create_se
        % create se and epos
        create_pre_post
        % create structures pre, post, redcell, nframe, npair
        set_movie_frames
        % set movie frames for pre and post, create movie_frames_pre,
        % movie_frames_post, movie_sect_length
        plot_pca
        % perform and plot principal component analysis for each mouse pre and post
        plot_pca_corr
        % plot correlation of pca pre vs post
        plot_skewness
        % plot skewness pre vs post for all mice 
        plot_skewness_popn_corr
        % plot skewness change vs population correlation for all mice
        plot_skewness_rc
        % plot skewness change of redcell vs non redcell for all mice
        plot_skewness_dist_prepost
        % plot skewness change vs distance from electrode for all mice


% plot_two_session
% plotRG


%% making movies
sbxmovie
sbxmovie_all