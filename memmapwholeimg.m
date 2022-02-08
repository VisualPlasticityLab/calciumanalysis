function fname=memmap(nam,magnification,picsiz)
%% load file
%path_to_package = 'C:\Users\kraken\Documents\MATLAB\ca_source_extraction\utilities';   % path to the folder that contains the package
%addpath(genpath(path_to_package));
%close all;

starttime=datetime;
data = matfile([nam '_memmap.mat'],'Writable',true);
if nargin<2
    try
        magnification=data.magnification;
    catch
        magnification=2;
    end
end
% if nargin<3
%     if numel(strfind(nam,'_'))<3
%         picsiz=[1 768 1 796] + [5 -5 5 -5];
%     else
%         picsiz=[100 760 10 786];
%     %two planes [100 760 10 786]
%     %bidirectional and two planes [150 760 134 786]
%     end
% end

%% Set parameters
%tic
% try
%     sbxread(nam,1,1);global info;magnification=info.config.magnification;
% catch
%     filenam=dir('*.sbx');
%     display(filenam(1).name);
%     sbxread(filenam(1).name(1:end-4),1,1);
% end

try
    sizY=data.sizY;
catch
    sizY = size(data,'Y');                  % size of data matrix
end
if nargin<3
    if numel(strfind(nam,'_'))<3
        picsiz=[1 sizY(1) 1 sizY(2)] + [5 -5 5 -5];
    else
        picsiz=[1 sizY(1) 10 sizY(2)]+ [5 -5 5 -5]*2;
    %picsiz=[100 sizY(1) 10 sizY(1)]
    %bidirectional and two planes [150 760 134 786]
    end
end

pcsX = picsiz(2)-picsiz(1) ;
pcsY = picsiz(4)-picsiz(3);
% sq=[180 300];
sq=[80 120];
overlap = [3,5]*ceil(magnification*1.5);                        % amount of overlap in each dimension (optional, default: [4,4])
patch_size(1) = ceil(mod(pcsX,sq(1))/floor(pcsX/sq(1)))+sq(1)+overlap(1); 
patch_size(2) = ceil(mod(pcsY,sq(2))/floor(pcsY/sq(2)))+sq(2)+overlap(2);
% size of each patch along each dimension (optional, default: [128,120])
patches = construct_partial_patches(picsiz,patch_size,overlap);

patch_size
display(sprintf('find %d cells /patch * %d patches',K,numel(patches)))
%K=ceil(15/magnification);
tau = ceil([3,5]*magnification);          % std of gaussian kernel (half width/height of neuron)
p = 0;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
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
    'max_size', ceil(6*magnification),...
    'merge_thr',merge_thr,...                    % merging threshold
    'gSig',tau...
    );
%     'fudge_factor',0.98,...                     % bias correction for AR coefficients

%% Run on whole picture
tic;
% fast initialization of spatial components using greedyROI and HALS
 % number of components to be found
K=100;
[Ain,Cin,bin,fin,center] = initialize_components(Y,K,tau,options);  % initialize

% display centers of found components
Cn =  reshape(P.sn,d1,d2); %correlation_image(Y); %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)
figure;imagesc(Cn);
    axis equal; axis tight; hold all;
    scatter(center(:,2),center(:,1),'mo');
    title('Center of ROIs found from initialization algorithm');
    drawnow;

% manually refine components (optional)
refine_components = false;  % flag for manual refinement
if refine_components
    [Ain,Cin,center] = manually_refine_components(Y,Ain,Cin,center,Cn,tau,options);
end
    
% update spatial components
Yr = reshape(Y,d,T);
clear Y;
[A,b,Cin] = update_spatial_components(Yr,Cin,fin,Ain,P,options);

% update temporal components
[C,f,P,S] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

%% merge found components
[Am,Cm,K_m,merged_ROIs,P,Sm] = merge_components(Yr,A,b,C,f,P,S,options);

%% repeat
[A2,b2,Cm] = update_spatial_components(Yr,Cm,f,Am,P,options);
[C2,f2,P,S2] = update_temporal_components(Yr,A2,b2,Cm,f,P,options);;

display(sprintf('run patches in %d min',round(toc/60)));

%% DONnotorder , just plot
contour_threshold = 0.95;                 % amount of energy used for each component to construct contour plot
K_m = size(C,1);
%[A_or,C_or,S_or,P] = order_ROIs(A,C,S,P); % order components
%     [sig,S_df] = extract_DF_F(data,[A_or,b],[C_or;f],K_m+1); % extract DF/F values (sig) and spiking info
[sig,S_df] = extract_DF_F(data,[A,b],[C;f],K_m+1); % extract DF/F values (sig) and spiking info
% try
%     V=data.V;
% catch
%     m=info.aligned.m;
%     V=zeros(sizY(1),sizY(2));
%     fac=sqrt(sizY(end));
%     for i=1:sizY(end)
%         V=V+((double(data.Y(:,:,i)-m))/fac).^2;
%     end
% end
V=reshape(P.sn,sizY(1),sizY(2));

%    [Coor,json_file] = plot_contours(A_or,V,contour_threshold,1);
    [Coor,json_file] = plot_contours(A,V,contour_threshold,1);
%% save signals into different files
resultname=sprintf('%s_%dcell_tau%d_%d_tsub%d',nam,K_m,tau(1),tau(2),tsub);
fname=strrep(resultname,'\\mps-zfs\data\jsun','\\mps-pc53\');
sig=sig';
title(strrep(resultname,'_','\_'));
save([fname '.signals'],'sig','S_df');

saveas(gcf,[fname '.fig'],'fig');
saveas(gcf,[fname '.png'],'png');
try
    T=0;
    for n=1:numel(data.eachsize)
        sig_chunk=sig(T+1:T+data.eachsize(1,n),:);
        sigs_chunk=sig(T+1:T+data.eachsize(1,n),:);
        sigsp_chunk=sig(T+1:T+data.eachsize(1,n),:);
        save([fname '_' num2str(n) '.signals'],'sig_chunk','S_df');
        T=T+data.eachsize(1,n);
    end
catch
end

endtime=datetime;
display(fname)
display(starttime);
display(endtime);
savejson('jmesh',json_file,[fname '.jmesh']);        % optional save json file with component coordinates (requires matlab json library)

% sendmail('j.suninchina@gmail.com','Finished project', ...
%          [fname 'finished']);