function newname=makememmap(fn)
startime=datetime;
path=pwd;
zfs_path=strrep(path,'c:','\\mps-zfs\data\jsun');
newname=fullfile(zfs_path,[fn '_memmap.mat']);
if exist(newname,'file')
    disp([' memmap file already existing']);return;
elseif ~exist(zfs_path,'dir')
    mkdir(zfs_path);
end

global info_loaded info
if(~isempty(info_loaded))   % try closing previous...
    try
        fclose(info.fid);
    catch
    end
end
sbxread(fn,1,1);        % read one frame to read the header of the image sequence
m=info.aligned.m;
Y=zeros([info.sz info.max_idx+1],'uint16');
V=zeros(info.sz);
fac=sqrt(info.max_idx);
tic
for i=0:info.max_idx
    z =sbxread(fn,i,1);
    z= squeeze(z);
    Y(:,:,i+1) = circshift(z,info.aligned.T(i+1,:)); % align the image
    temp=((double(Y(:,:,i+1)-m))/fac).^2;
    V=V+temp;
    if mod(i,1000)==0
        fprintf('Aligning frame #%d/%d\n',i,info.max_idx+1);toc;
    end
end

sizY = size(Y);
Yr = reshape(Y,prod(sizY(1:end-1)),[]);
nY = min(Yr(:));
save([fn '.align'],'V','-append');

if mod(info.config.magnification*2,5) == 0
    magnification = info.config.magnification*2/5;
else
    magnification = info.config.magnification;
end
disp(' try saving memmap file');
save(newname,'sizY','Yr','Y','nY','V','m','magnification','-v7.3');
display(startime);
display(datetime);

