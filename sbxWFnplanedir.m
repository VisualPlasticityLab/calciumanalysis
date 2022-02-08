function sbxWFnplanedir(varargin)
%

if(nargin>=1) % cell with filenames to be aligned
    for(i=1:numel(varargin))
        d(i).name = varargin{i};
    end
else
    d = dir('*.sbx');
end

% Align all *.sbx files in the list

for(i=1:length(d))
    try
        fn = strtok(d(i).name,'.');
            sbxread(fn,1,1);            % read one frame to read the header of the image sequence
            global info;                % this contains the information about the structure of the image
            tic
           0-vb                      00
            T=zeros(info.max_idx+1,2);%%%%CAN'T be uint16!!!
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
