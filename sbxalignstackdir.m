function sbxalignstackdir(varargin)
% chan: 1 for green, 2 for red
%

if(nargin>=1) % cell with filenames to be aligned
    for files=1:numel(varargin)
        d(files).name = varargin{files};
    end
else
    d = dir('*.sbx');
end

% Align all *.sbx files in the list

for(files=1:length(d))
    try
        fn = strtok(d(files).name,'.');
        sbxread(fn,1,1);            % read one frame to read the header of the image sequence
        global info;                % this contains the information about the structure of the image
        if ~info.volscan
            knobby = info.config.knobby.schedule(:,5);
            knobby = [ 0;knobby ; info.max_idx+1];
            nplanes=numel(knobby)-1;
            T=zeros(info.max_idx+1,2);%%%%CAN'T be uint16!!!
            fnametif = sprintf('%s_average.tif',fn);
            for n=1:nplanes
                [m(:,:,n),T(knobby(n)+1:knobby(n+1),:)] = sbxalignx(fn,knobby(n):knobby(n+1)-1);
                q = squeeze(m(:,:,n));
                if(n==1)
                    imwrite(q, fnametif,'tif');
                else
                    imwrite(q, fnametif,'tif','writemode','append');
                end
            end
            save([fn '.align'],'m','T');
            
            clear m T;
            display(sprintf('Done %s: Aligned %d images in %d min',fn,info.max_idx+1,round(toc/60)));
        else
            sprintf('Aborted. detected volscan, please use sbxalignplanedir.m')
        end
    catch
        sprintf('Could not align %s',fn)
    end
end
