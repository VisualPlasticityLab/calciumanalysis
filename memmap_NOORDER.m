function fname=memmap(nam,patchf);
%% load file
%path_to_package = 'C:\Users\kraken\Documents\MATLAB\ca_source_extraction-master\utilities';   % path to the folder that contains the package
%addpath(genpath(path_to_package));
%close all;
%nam='419_711_000';
if ~exist([nam '_memmap.mat'],'file')
	makememmap(nam);
end
data = matfile([nam '_memmap.mat'],'Writable',true);

%% Set parameters

tic
sbxread(nam,1,1);
global info;
sizY = size(data,'Y');                  % size of data matrix
if ~exist('patchf')
    mag=ceil(sqrt(sizY(3)/16000)); % if too many frames, chop the patch sizes to smaller ones
    if mod(sizY(1),128)==0      
        patchf=[ceil(sizY(1)*mag/128) ceil(sizY(2)*mag/100)];
    else
        patchf=[ceil(sizY(1)*mag/100) ceil(sizY(2)*mag/100)];
    end
end
patchf
patch_size = [ceil(sizY(1)/patchf(1)*1.2),ceil(sizY(2)/patchf(2)*1.2)];                   % size of each patch along each dimension (optional, default: [32,32])
overlap = [ceil(sizY(1)/patchf(1)*.1),ceil(sizY(2)/patchf(2)*.1)];                        % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1:end-1),patch_size,overlap);
K = ceil(sizY(1)*sizY(2)/10000/prod(patchf))    % number of components to
%be found per patch
K=ceil(100/prod(patchf));
tau = 2;                                          % std of gaussian kernel (size of neuron) 
p = 0;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.7;                                  % merging threshold
sizY = data.sizY;

options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'search_method','ellipse','dist',3,...      % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'ssub',info.config.magnification,...
    'tsub',2,...
    'fudge_factor',0.98,...                     % bias correction for AR coefficients
    'merge_thr',merge_thr,...                    % merging threshold
    'gSig',tau...
    );
display(sprintf('CNMFparms done %d min',round(toc/60)));

%% Run on patches
tic;

[A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(data,K,patches,tau,p,options);

display(sprintf('run patches in %d min',round(toc/60)));

%% DONNOT  order , just plot
%[A_or,C_or,S_or,P] = order_ROIs(A,C,S,P); % order components
K_m = size(C,1);
fname=sprintf('%s_%dcell',[nam '_memmap'],K_m);
tic
%[sig,S_df] = extract_DF_F(data,[A_or,b],[C_or;f],K_m+1); % extract DF/F values (sig) and spiking info
[sig,S_df] = extract_DF_F(data,[A,b],[C;f],K_m+1); % extract DF/F values (sig) and spiking info

display(sprintf('extracted sig in %d min',round(toc/60)));
sig=sig';
save([fname '.signals'],'sig','S_df');

%[Coor,json_file] = plot_contours(A_or,reshape(P.sn,sizY(1),sizY(2)),contour_threshold,1); % contour plot of spatial footprints
if ~isfield(info.aligned,'V')
    V=sbxVariancemap(nam);
else
    V=info.aligned.V;
end
contour_threshold = 0.95;                 % amount of energy used for each component to construct contour plot
%[Coor,json_file] = plot_contours(A_or,V,contour_threshold,1); 
[Coor,json_file] = plot_contours(A,V,contour_threshold,1); 
h=gcf;
h.Name=fname;
savefig(h,fname);
savejson('jmesh',json_file,[fname '.jmesh']);        % optional save json file with component coordinates (requires matlab json library)
