function [m,mg,mr,T] = sbxalignx2(fname,idx,chan)

% Aligns images in fname for all indices in idx
% 
% m - mean image after the alignment
% T - optimal translation for each frame
if nargin<3 
    chan =1; %default to use green channel
end

if(length(idx)==1)
    
    A0 = sbxread(fname,idx(1),1);
    m = squeeze(A0(chan,:,:));
    mg = squeeze(A0(1,:,:));
    mr = squeeze(A0(2,:,:));
    T = [0 0];
    
elseif (length(idx)==2)
    
    A0 = sbxread(fname,idx(1),1);
    B0 = sbxread(fname,idx(2),1);
    A = squeeze(A0(chan,:,:));
    B = squeeze(B0(chan,:,:));
    gA = squeeze(A0(1,:,:));
    gB = squeeze(B0(1,:,:));
    rA = squeeze(A0(2,:,:));
    rB = squeeze(B0(2,:,:));
    [u v] = fftalign(A,B);
    
    m = (circshift(A,[u,v])+B)/2;
    mg = (circshift(gA,[u,v])+gB)/2;
        mr = (circshift(rA,[u,v])+rB)/2;

    T = [[u v] ; [0 0]];
    
else
    
    idx0 = idx(1:floor(end/2));
    idx1 = idx(floor(end/2)+1 : end);
    
    [A,gA,rA,T0] = sbxalignx2(fname,idx0,chan);
    [B,gB,rB,T1] = sbxalignx2(fname,idx1,chan);
   
    [u v] = fftalign(A,B);
     
    m = (circshift(A,[u, v])+B)/2;
    mg = (circshift(gA,[u,v])+gB)/2;
        mr = (circshift(rA,[u,v])+rB)/2;

    T = [(ones(size(T0,1),1)*[u v] + T0) ; T1];
    
end
