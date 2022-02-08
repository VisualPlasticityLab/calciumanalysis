function sbx2tif_green_and_red(fname,varargin)

% sbx2tif
% Generates tif file from sbx files
% Argument is the number of frames to convert
% If no argument is passed the whole file is written

z = sbxread(fname,2,1);
global info;

if(nargin>1)
    N = min(varargin{1},info.max_idx);
else
    N = info.max_idx;
end
if info.volscan 
    nplane = info.otparam(3);
else
    nplane = 1;
end
newN = N - mod(N+1,nplane) ; %make sure equal number for each plane
maxk = ceil(4.252*10e8/info.nsamples) - 1;  % max tiff file size is 2608 frames for 2*512*796 or 1738 frames 768 *796 so find first number smaller than that:


k = 0;
numfiles = 1;
done = 0;
while k <= newN && ~done
    try 
    eachK = 0;
    while(eachK<maxk && k<=newN)
        q = sbxread(fname,k,1);
        qg = squeeze(q(1,:,:));
        qr = squeeze(q(2,:,:));
        if(k==0)
            imwrite(qg,[fname,'_',num2str(numfiles),'gr', '.tif'],'tif');
            imwrite(qr,[fname,'_',num2str(numfiles),'gr', '.tif'],'tif','writemode','append');
        else
            imwrite(qg,[fname,'_',num2str(numfiles),'gr', '.tif'],'tif','writemode','append');
            imwrite(qr,[fname,'_',num2str(numfiles),'gr', '.tif'],'tif','writemode','append');
        end
        if mod(k,500)==0
            fprintf('%d/%d complete...\n',k,N)
        end
        k = k+1;
        eachK = eachK+1;
    end
    fprintf('close file %d at k = %d\n',numfiles,k)
    eachk=0;
    numfiles = numfiles + 1;
    catch
        fprintf('finished at k= %d\n',k);
        done = 1;
    end
end
