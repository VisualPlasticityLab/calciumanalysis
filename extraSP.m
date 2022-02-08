function signew=extraSP(sig,fname,fast)

ncell = size(sig,2);
switch fast
    case 0 %using deconv and MLspike to extract Calcium signal
        load(fname,'sigML','-mat');
            if ~exist('sigML','var') %using deconvolve to extract Calcium signal
            fprintf('MLspiking,size %d %d,expecting %d sec\n',size(sig),round(prod(size(sig))/600));tic
            sigML = nan(size(sig));
            sigs = nan(size(sig));
            sigsp = nan(size(sig));
            for ii=1:ncell
                try
                    sigML(:,ii)=MLspike(sig(:,ii));
                catch
                         fprintf('cell %d did not process',ii);
                end               
                if mod(ii,50)==1
                    fprintf('MLspiking %d out of %d cells in %d secs\n',ii,ncell,round(toc))
                end
            end
            fprintf('Done deconvolveCa and MLspiking, took %d secs\n',round(toc))
%              save(fname,'sigML','-append');
        end
        signew = sigML;
    case .5
        load(fname,'sigsp','-mat');
        if ~exist('sigsp','var') %using deconvolve to extract Calcium signal
            fprintf('deconvlve...,size %d %d\n',size(sig));tic
            sigML = nan(size(sig));
            sigs = nan(size(sig));
            sigsp = nan(size(sig));
            for ii=1:ncell
                try
                    [sigs(:,ii),sigsp(:,ii)]=deconvolveCa(sig(:,ii));
                catch
                    fprintf('cell %d did not process',ii);
                end
                if mod(ii,50)==10
                    fprintf('deConvolveCa %d out of %d cells in %d secs\n',ii,ncell,round(toc))
                end
            end
            fprintf('Done deConvolveCa, took %d secs\n',round(toc))
%             save(fname,'sigsp','sigs','-append');
        end
        signew = sigsp;
    case 1
        fprintf('skip spiking extraction\n')
        signew = sig;
end