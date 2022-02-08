function sbxnplane2tif(fn)

if nargin ==0
    [fname,p]=uigetfile('.align','load aligned data');
    fn = strtok(fname,'.');
else
    [p,fnam,ext]=fileparts(fn);
    fname=[fnam ext];
end

load(fullfile(p,fname),'m','-mat');

N=size(m,3);
q = squeeze(m(:,:,1));
imwrite(q,[fn '.tif'],'tif');

for i=2:N
    q = squeeze(m(:,:,i));
    imwrite(q,[fn '.tif'],'tif','writemode','append');
end