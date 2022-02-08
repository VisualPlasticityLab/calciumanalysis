function sbx2tifnbinning(fname,varargin)
% (fname,nbinning = 1, N = info.max_idx);
% sbx2tif
% Generates tif file from sbx files
% Argument is the number of frames to convert
% If no argument is passed the whole file is written

z = sbxread(fname,1,1);
global info;

if(nargin>1)
    Nbinning = varargin{2};
else
    Nbinning = 1;
end

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
newN = N - mod(N+1,nplane) ; 
% maxk = 5219;  % max tiff file size is 5224 frames so find first number smaller than that that is divisible by 3:
% (5219+1)/6 = 1740/2 = 870
% maxklist = maxk*[1:1:10];

k = 0;
numfiles = 1;
while k <= newN
    done = 0;
    tifname = [fname,'_',num2str(numfiles), '_ aligned.tif'];
    startfile = 1;
    while(~done && k<=newN)
        try
            q = sbxread(fname,k,Nbinning*nplane);
%             q = mean(q(:,:,:,1:nplane:Nbinning*nplane),4); %binning them to get better signal
            q = squeeze(q(1,:,:));
%             if ~isempty(info.aligned)
%                 q = circshift(q,info.aligned.T(k+1,:));%+[9 53]); % align the image
%             end
            if startfile == 1
                imwrite(q,tifname,'tif');
                startfile = 0;
            else
                imwrite(q,tifname,'tif','writemode','append');
            end
            k = k+Nbinning;
        catch
            fprintf('Had to close file %d at k = %d\n',numfiles,k)
            numfiles = numfiles + 1;
            done = 1;
        end
        if mod(k,100*Nbinning)==1
            fprintf('%d/%d complete...\n',k+1,newN+1)
        end
%          if ismember(k,maxklist)
%             fprintf('Had to close file %d at k = %d\n',numfiles,k)
%             numfiles = numfiles + 1;
%             done = 1;
%          end
    end
end
