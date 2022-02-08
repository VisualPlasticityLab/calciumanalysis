function [m,T] = sbxalign_nfiles(varargin)

% Aligns images for multiple files
% 
% m - mean image after the alignment
% T - optimal translation for each frame


sbxaligndir;
ref=varargin{1};
[m0,T0] = load([ref '.align']);

for i=2:nargin
    fn=varargin{i};
    [m1,T1] = load([fn '.align']);
    [u v] = fftalign(m1,m0);     
    m= circshift(m1,[u, v]);
    T = ones(size(T1,1),1)*[u v] + T1;
    save([fn 'to' ref(ref~=fn) '.align'],'m','T')
end
