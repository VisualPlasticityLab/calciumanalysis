function sbx2tif(fname,varargin)
% (fname,nbinning = 1, N = info.max_idx);
% sbx2tif
% Generates tif file from sbx files
% Argument is the number of frames to convert
% If no argument is passed the whole file is written

z = sbxread(fname,1,1);
global info;


if(nargin>2)
    N = min(varargin{1},info.max_idx);
else
    N = info.max_idx;
end

% Make sure newN+1 is divisible by 2 or 3 planes!
if info.volscan 
    nplane = info.otparam(3);
else
    nplane = 1;
end
newN = N - mod(N+1,nplane) ; %make sure equal number for each plane
eachimage = whos('z');
maxk = ceil(4.252*10e8/eachimage.bytes) - 1;  % max tiff file size is 5216 frames for 512*796 or 3477 frames 768 *796 so find first number smaller than that that is divisible by 3:
% (5219+1)/6 = 1740/2 = 870
% maxklist = maxk*[1:1:10];

k = 0;
numfiles = 1;
while k <= newN
    tifname = [fname,'_',num2str(numfiles), '.tif'];
    startfile = 1;
    eachk = 0;
    while (eachk<=maxk && k<=newN) %within one tif file
        q = sbxread(fname,k,1);
        q = squeeze(q(1,:,:));
        if startfile == 1
            imwrite(q,tifname,'tif');
            startfile = 0;
        else
            imwrite(q,tifname,'tif','writemode','append');
        end
        k = k+1;
        eachk = eachk+1;
        if mod(k+1,500)==0
            fprintf('%d/%d complete,file%s\n',k+1,newN+1,tifname)
        end
    end
    fprintf('close file %d at k = %d\n',numfiles,k)
    numfiles = numfiles + 1;
end
