clear all;close all;
%% load file
% function memmapmanuel(fn);             
fn='1793_329_000_1_memmap.mat';
data = matfile(fn);

tic;
%% Set parameters

K = 40;                                           % number of components to be found
tau = [4,6]*2;                                          % std of gaussian kernel (size of neuron) 
p = 2;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
magnification=data.magnification;
sizY=data.sizY;
merge_thr = 0.9;                                  % merging threshold
tsub=20;
options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'search_method','ellipse','dist',3,...      % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps
    'ssub',magnification,...
    'tsub',tsub,...
    'save_memory',1,...
    'max_size', ceil(6*data.magnification),...
    'merge_thr',merge_thr,...                    % merging threshold
    'gSig',tau...
    );
%% fast initialization of spatial components using greedyROI and HALS
P= preprocess_data(data.Y,p);
[Ain,Cin,bin,fin,center] = initialize_components(data.Y,K,tau,options);  % initialize

% display centers of found components
Cn =  reshape(P.sn,d1,d2); %correlation_image(Y); %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)
figure;imagesc(Cn);
    axis equal; axis tight; hold all;
    scatter(center(:,2),center(:,1),'mo');
    title('Center of ROIs found from initialization algorithm');
    drawnow;

%% manually refine components (optional)
refine_components = true;  % flag for manual refinement
if refine_components
    [Ain,Cin,center] = manually_refine_components(data.Y,Ain,Cin,center,Cn,tau,options);
end
    
%% update spatial components

[A,b,Cin] = update_spatial_components(data.Yr,Cin,fin,Ain,P,options);

%% update temporal components
[C,f,P,S] = update_temporal_components(data.Yr,A,b,Cin,fin,P,options);

%% merge found components
[Am,Cm,K_m,merged_ROIs,P,Sm] = merge_components(data.Yr,A,b,C,f,P,S,options);

%% repeat
[A2,b2,Cm] = update_spatial_components(data.Yr,Cm,f,Am,P,options);
[C2,f2,P,S2] = update_temporal_components(data.Yr,A2,b2,Cm,f,P,options);

%% do some plotting/output
[A_or,C_or,S_or,P] = order_ROIs(A2,C2,S2,P); % order components
K_m = size(C_or,1);

[C_df,~,S_df] = extract_DF_F(Yr,[A_or,b2],[C_or;f2],S_or,K_m+1); % extract DF/F values (optional)
elapsedtime=toc
%%
fname=sprintf('%s_%dseg_%d-%d',nam,K,sframe,sframe+num2read-1);
save([fname '.dF'],'C_df','S_df','elapsedtime');
contour_threshold = 0.95;                       % amount of energy used for each component to construct contour plot
h=figure;
[Coor,json_file] = plot_contours(A_or,reshape(P.sn,d1,d2),contour_threshold,1); % contour plot of spatial footprints
savefig(h,fname);
savejson('jmesh',json_file,[fname '.jmesh']);        % optional save json file with component coordinates (requires matlab json library)
%% display components
plot_components_GUI(Yr,A_or,C_or,b2,f2,Cn,options);
