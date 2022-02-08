function fnames=memmap(nam,magnification,picsiz)
%% load file
%path_to_package = 'C:\Users\kraken\Documents\MATLAB\ca_source_extraction\utilities';   % path to the folder that contains the package
%addpath(genpath(path_to_package));
%close all;

starttime=datetime;
data = matfile([nam '_memmap.mat'],'Writable',true);
if nargin<2
    try
        load('magnificationlist');
        if mod(data.magnification,1)>0 || data.magnification==2
            data_magnification = data.magnification*2.5;
        else
            data_magnification = data.magnification;
        end
        magnification=maglist(data_magnification);
    catch
        magnification=1;
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

 sq=[240 300]/3;
%sq=ceil([20 30]*magnification);
% magnification = magnification/1.5;
% sq=ceil([50 78]*magnification);
overlap = ceil([3,5]*magnification);  % amount of overlap in each dimension

pcsX = picsiz(2)-picsiz(1) ;
pcsY = picsiz(4)-picsiz(3);
patch_size(1) = ceil(mod(pcsX,sq(1))/floor(pcsX/sq(1))+sq(1)+overlap(1)); 
patch_size(2) = ceil(mod(pcsY,sq(2))/floor(pcsY/sq(2))+sq(2)+overlap(2));
% size of each patch along each dimension (optional, default: [128,120])
patches = construct_partial_patches(picsiz,patch_size,overlap);
%K = ceil(sizY(1)*sizY(2)/10000/prod(patchf))    % number of components to be found
%K=round(prod(patch_size)/400/magnification^2);
K= 6;
patch_size
display(sprintf('find %d cells per patch * %d patches',K,numel(patches)))

tau = round([2,4]*magnification);          % std of gaussian kernel (half width/height of neuron)
p = 0;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.95;                                  % merging threshold
tsub=10;
options = CNMFSetParms(...
        'd1',sizY(1),'d2',sizY(2),...
    'search_method','ellipse','dist',3,...      % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps
    'ssub',magnification,...
    'tsub',tsub,...
    'save_memory',1,...
    'max_size', ceil(10*magnification),...
    'merge_thr',merge_thr,...                    % merging threshold
    'gSig',tau...
    );
%     'fudge_factor',0.98,...                     % bias correction for AR coefficients

%display(sprintf('CNMFparms done %d min',round(toc/60)));

%% Run on patches
tic;

[A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(data,K,patches,tau,p,options);

display(sprintf('run patches in %d min',round(toc/60)));

%% DONnotorder , just plot
contour_threshold = 0.95;                 % amount of energy used for each component to construct contour plot
K_m = size(C,1);
%[A_or,C_or,S_or,P] = order_ROIs(A,C,S,P); % order components
%     [sig,S_df] = extract_DF_F(data,[A_or,b],[C_or;f],K_m+1); % extract DF/F values (sig) and spiking info
[sig,~] = extract_DF_F(data,[A,b],[C;f],K_m+1); % extract DF/F values (sig) and spiking info
try
    V=data.V;
catch
    m=data.m;
    V=zeros(sizY(1),sizY(2));
    fac=max(m(:));
    for i=1:sizY(end)
        V=V+abs(double(data.Y(:,:,i)-m)/fac);
    end
    data.V=V;
end
V=reshape(P.sn,sizY(1),sizY(2));

%    [Coor,json_file] = plot_contours(A_or,V,contour_threshold,1);
    [Coor,json_file] = plot_contours(A,V,contour_threshold,1);
%% save signals and cell locations into different files
fname=sprintf('%s_%dcell_tau%d_%d_tsub%d',nam,K_m,tau(1),tau(2),tsub);
sig=sig';
% [sig, sigsp, options] = deconvolveCa(y, varargin)
title(strrep(fname,'_','\_'));

factor.picsiz = picsiz;
factor.K = K;
factor.sq =sq;
save([fname '.signals'],'sig','factor','A','V','Coor','-v7.3');
savejson('jmesh',json_file,[fname '.jmesh']);        % optional save json file with component coordinates (requires matlab json library)

try
    saveas(gcf,[fname '.png'],'png');
    saveas(gcf,[fname '.fig'],'fig');
    close all;
end



if contains(nam,'&') % derive spike signals for multiple files
    fnames = splitfile(fname);
    sigs = nan(size(sig));
    sigsp = nan(size(sig));
%     for i=1:K_m
%         try
%             [sigs(:,i),sigsp(:,i)]=deconvolveCa(sig(:,i));
%         catch
%             sprintf('Warning! Cell %d did not process',i)
%         end
%     end
        eachsize=data.eachsize;
        T=0;
        sigALL=sig;
        sigsALL=sigs;
        sigspALL=sigsp;
        for n=1:numel(eachsize)
            sig=sigALL(T+1:T+eachsize(1,n),:);
            sigs=sigsALL(T+1:T+eachsize(1,n),:);
            sigsp=sigspALL(T+1:T+eachsize(1,n),:);
            save([fnames{n} '.signals'],'sig','sigs','sigsp','factor','A','V','Coor');
            T=T+eachsize(1,n);
        end
else
            fnames = fname;
end

endtime=datetime;
display(fname)
display(starttime);
display(endtime);

setpref('Internet','SMTP_Server','keck.ucsf.edu')
setpref('Internet','E_mail','jsun@phy.ucsf.edu')
sendmail({'j.suninchina@gmail.com'},'Finished project', ...
         [fname 'finished']); %,'mdadarla@phy.ucsf.edu'