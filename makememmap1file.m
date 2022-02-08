function newname=makememmap1file(fn)
%% define the name and save to MPS-ZFS directly
path=pwd;
if strcmpi(path(1),'C')==1
    zfs_path=strrep(path,path(1:2),'\\mps-zfs\data1\jsun');
else
    zfs_path=strrep(path,'\\MPS-PC53','\\mps-zfs\data1\jsun');
end
newname=[fn '_memmap.mat'];
newname=fullfile(zfs_path,newname);
% 
% if exist(newname,'file')
%     disp([newname ' memmap file already existing']);return;
% elseif ~exist(zfs_path,'dir')
%     mkdir(zfs_path);
% end
data = matfile(newname,'Writable',true);
display(sprintf(' try saving memmapfile: %s',newname));
starttime=datestr(now);
%% save 1 piece
tic
global info;            % this contains the information about the structure of the image
sbxread(fn,1,1);% read one frame to read the header of the image sequence
T0=info.max_idx+1;
m=info.aligned.m;
V=zeros(info.sz);
fac=sqrt(info.max_idx);

data.m=zeros(info.sz);
data.Y=zeros([info.sz 2],'uint16');
data.Yr = zeros([prod(info.sz) 2],'uint16');
%data.V=zeros(info.sz);
data.sizY =[info.sz T0];
data.nY = inf ;
data.m=m;
chunksize=500;
for i=0:chunksize:T0-1
    for j=i:min(chunksize-1+i,T0-1)
    z =sbxread(fn,j,1);    z= squeeze(z(1,:,:));
    img = circshift(z,info.aligned.T(j+1,:)); % align the image
    data.Y(:,:,j+1) = img;
    V=V+((double(img-m))/fac).^2;
    end
    data.Yr(:,i+1:j+1) = reshape(data.Y(:,:,i+1:j+1),prod(info.sz),j-i+1);
    fprintf('File %d Frame %d/%d for %.2f seconds\n ',i,j,T0-1,toc);
    nY= min(min(data.Yr(:,i+1:j+1)));
    if nY< data.nY
        data.nY = nY;
    end
end
save([fn '.align'],'V','-append');
data.V=V/max(V(:))*double(max(m(:)));

%data.Yr = reshape(data.Y,prod(info.sz),T0);
%data.nY = min(reshape(data.Yr,prod(info.sz)*T0,1));
if info.config.magnification>=5 
    data.magnification = info.config.magnification*2/5;
else
    data.magnification = info.config.magnification;
end
display(sprintf('Finished making memmap %s',newname));
display(starttime)
display(datestr(now))

