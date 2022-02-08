function makememmapdir(varargin)

if(nargin>=1) % cell with filenames to be aligned
    for(i=1:length(varargin{1}))
        d(i).name = varargin{1}{i};
    end
else
    d = dir('*.sbx');
end

path=pwd;
zfs_path=strrep(path,'C:','\\mps-zfs\data\jsun');

for(i=1:numel(d))
    fn = strtok(d(i).name,'.')
    if ~ exist(fullfile( zfs_path,[fn '_memmap.mat']))
        tic
        makememmap(fn);
        display(sprintf('saved memmap file in %d min',round(toc/60)));
    else
        sprintf('%s_memmap.mat is already made',fn)
    end
end