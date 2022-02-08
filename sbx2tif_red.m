function sbx2tif_red(fname,varargin)

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

% Make sure N+1 is divisible by 3 planes!
newN = N - mod(N+1,3);
maxk = 5219; % max tiff file size is 5224 frames (k = 0:2611) so find first number smaller than that that is divisible by 3:
% (2609+1)*2 = 5220; 5220/3 = 1740
maxklist = maxk*[1:10];

k = 0;

numfiles = 1;
while k <= newN
    done = 0;
    while(~done && k<=newN)
        try
            q = sbxread(fname,k,1);
            %qg = squeeze(q(1,:,:));
            qr = squeeze(q(2,:,:));
            if(k==0)
                imwrite(qr,[fname,'_',num2str(numfiles),'_r', '.tif'],'tif');
            else
                imwrite(qr,[fname,'_',num2str(numfiles),'_r', '.tif'],'tif','writemode','append');
            end
            
        catch
            fprintf('Had to close file %d at k = %d\n',numfiles,k)
            numfiles = numfiles + 1;
            done = 1;
        end
        if mod(k,1000)==0
            fprintf('%d/%d complete...\n',k,N)
        end
        
        if ismember(k,maxklist)
            fprintf('Had to close file %d at k = %d\n',numfiles,k)
            numfiles = numfiles + 1;
            done = 1;
        end
        k = k+1;
    end
end
