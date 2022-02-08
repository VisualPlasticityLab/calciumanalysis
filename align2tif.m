function align2tif(fname,varargin)

% align2tif
% Generates tif file from .align files
% Argument is the number of frames to convert
% If no argument is passed the whole file is written

if nargin ==0
    [fname,p]=uigetfile('.align','load aligned data');
else
    [p,fnam,ext]=fileparts(fn);
    fname=[fnam ext];
end

load(fullfile(p,fname),'m','-mat');
tifname = [strtok(fname,'.') '_average.tif'];
nplane=size(m,3);
for k=1:nplane
    q = squeeze(m(:,:,k));
    if(k==1)
        imwrite(q,tifname,'tif');
    else
        imwrite(q,tifname,'tif','writemode','append');
    end
end
