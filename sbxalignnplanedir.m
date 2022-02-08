function sbxalignnplanedir(chan,varargin)
% chan: 1 for green, 2 for red
%

if(nargin>=2) % cell with filenames to be aligned
    for(i=1:numel(varargin))
        d(i).name = varargin{i};
    end
elseif nargin>=1
    %specify the green/red channel for alignment
    d = dir('*.sbx');
else
    chan =0;  % default to use green channel(1) for alignment
    d = dir('*.sbx');
end

% Align all *.sbx files in the list

for(i=1:length(d))
    try
        fn = strtok(d(i).name,'.');
         if ~ exist([fn '.align'])
            sbxread(fn,1,1);            % read one frame to read the header of the image sequence
            global info;                % this contains the information about the structure of the image
            tic
            if info.volscan ==1
                nplanes=info.otparam(3);
            else
                nplanes=1;
            end
            T=zeros(info.max_idx+1,2);%%%%CAN'T be uint16!!!
            display(sprintf('%s:Aligning %d images',fn,info.max_idx+1));
            if  info.nchan >1 
                if  chan ==0
                chan = menu('Use which channel to align?','Green','Red')
                end
                for i=1:nplanes
                    [m(:,:,i),mg(:,:,i),mr(:,:,i),T(i:nplanes:info.max_idx+1,:)] = sbxalignx2(fn,i-1:nplanes:info.max_idx,chan);
                end
                save([fn '.align'],'m','mg','mr','T');
                clear m mg mr T;
            else
                for i=1:nplanes
                    [m(:,:,i),T(i:nplanes:info.max_idx+1,:)] = sbxalignx(fn,i-1:nplanes:info.max_idx);
                end
                save([fn '.align'],'m','T');
                clear m T;
                display(sprintf('Done %s: Aligned %d images in %d min',fn,info.max_idx+1,round(toc/60)));
            end
        else
            sprintf('File %s is already aligned',fn)
        end
    catch
        sprintf('Could not align %s',fn)
    end
end

function [m,T] = sbxalignx(fname,idx)

% Aligns images in fname for all indices in idx
% 
% m - mean image after the alignment
% T - optimal translation for each frame

if(length(idx)==1)
    
    A = sbxread(fname,idx(1),1);
    A = squeeze(A(1,:,:));
    m = A;
    T = [0 0];
    
elseif (length(idx)==2)
    
    A = sbxread(fname,idx(1),1);
    B = sbxread(fname,idx(2),1);
    A = squeeze(A(1,:,:));
    B = squeeze(B(1,:,:));
    
    [u v] = fftalign(A,B);
    
    Ar = circshift(A,[u,v]);
    m = (Ar+B)/2;
    
    if sqrt(u^2+v^2)>25*2
        disp(idx);u=0;v=0;
    end  
    T = [[u v] ; [0 0]];
else
    
    idx0 = idx(1:floor(end/2));
    idx1 = idx(floor(end/2)+1 : end);
    
    [A,T0] = sbxalignx(fname,idx0);
    [B,T1] = sbxalignx(fname,idx1);
   
    [u v] = fftalign(A,B);
    if sqrt(u^2+v^2)>25*2
        disp(idx);u=0;v=0;
    end    
    Ar = circshift(A,[u, v]);
    m = (Ar+B)/2;
    T = [(ones(size(T0,1),1)*[u v] + T0) ; T1];
    
end

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

% function [m,T] = sbxalignx(fname,idx)
% 
% Aligns images in fname for all indices in idx
% 
% m - mean image after the alignment
% T - optimal translation for each frame
% 
% if(length(idx)==1)
%     
%     A = sbxread(fname,idx(1),1);
%     A = squeeze(A(1,:,:));
%     m = A;
%     T = [0 0];
%     
% elseif (length(idx)==2)
%     
%     A = sbxread(fname,idx(1),1);
%     B = sbxread(fname,idx(2),1);
%     A = squeeze(A(1,:,:));
%     B = squeeze(B(1,:,:));
%     
%     [u v] = fftalign(A,B);
%     
%     Ar = circshift(A,[u,v]);
%     m = (Ar+B)/2;
%     
%     if sqrt(u^2+v^2)>25*2
%         disp(idx);u=0;v=0;
%     end  
%     T = [[u v] ; [0 0]];
% else
%     
%     idx0 = idx(1:floor(end/2));
%     idx1 = idx(floor(end/2)+1 : end);
%     
%     [A,T0] = sbxalignx(fname,idx0);
%     [B,T1] = sbxalignx(fname,idx1);
%    
%     [u v] = fftalign(A,B);
%     if sqrt(u^2+v^2)>25*2
%         disp(idx);u=0;v=0;
%     end    
%     Ar = circshift(A,[u, v]);
%     m = (Ar+B)/2;
%     T = [(ones(size(T0,1),1)*[u v] + T0) ; T1];
%     
% end

