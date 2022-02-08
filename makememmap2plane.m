function newnames=makememmap2plane(fn)
%% define the name
path=pwd;
if strcmpi(path(1),'C')==1
    zfs_path=strrep(path,path(1:2),'\\mps-zfs\data1\jsun');
else
    zfs_path=strrep(path,'\\mps-pc53','\\mps-zfs\data1\jsun');
end
newnames{1}=fullfile(zfs_path,[fn '_1_memmap.mat']);
newnames{2}=fullfile(zfs_path,[fn '_2_memmap.mat']);

if exist(newnames{1},'file')
    display('memmap file already exist!'); return;
elseif ~exist(zfs_path,'dir')
    mkdir(zfs_path);
end

data1 = matfile(newnames{1},'Writable',true);
data2 = matfile(newnames{2},'Writable',true);

disp(sprintf('Saving memmapfiles for 2 planes:%s,%s',newnames{1},newnames{2}));

%% save each piece
sbxread(fn,1,1);% read one frame to read the header of the image sequence
global info;

nY = Inf;
m1=info.aligned.m(:,:,1);
m2=info.aligned.m(:,:,2);
V1=zeros(info.sz);
V2=zeros(info.sz);
T2=floor((info.max_idx+1)/2);
T1=info.max_idx+1-T2;
fac=sqrt(info.max_idx);
data1.Y=zeros([info.sz 2],'uint16');
data2.Y=zeros([info.sz 2],'uint16');
data1.m = m1;
data2.m = m2;
data1.sizY =[info.sz T1];
data2.sizY =[info.sz T2];

% data1.Yr = zeros([prod(info.sz) T1],'uint16');
% data2.Yr = zeros([prod(info.sz) T2],'uint16');
% 
% data1.nY = inf;
% data2.nY = inf;
tic
for j=0:info.max_idx
    z =sbxread(fn,j,1);
    z= squeeze(z);
    if mod(j,2)==0
        img = circshift(z,info.aligned.T(j+1,:)); % align the image
        data1.Y(:,:,j/2+1) = img;
%         data1.Yr(:,j/2+1) = reshape(img,:,1);
%         if data1.nY > min(data1.Yr(:,j/2+1))
%             data1.nY = min(data1.Yr(:,j/2+1));
%         end
        V1=V1+((double(img-m1))/fac).^2;
    else
        img = circshift(z,info.aligned.T(j+1,:)); % align the image
        data2.Y(:,:,(j+1)/2) = img;
%         data2.Yr(:,(j+1)/2) = reshape(img,:,1);
%         if data2.nY > min(data2.Yr(:,(j+1)/2))
%             data2.nY = min(data2.Yr(:,(j+1)/2));
%         end
        V2=V2+((double(img-m2))/fac).^2;
    end
    
    if mod(j,500)==0
        fprintf('Saved Frame %d/%d for %.2f seconds\n ',j,info.max_idx,toc);
    end
end
%V=cat(3,V1,V2); save([fn '.align'],'V','-append');


data1.V=V1;
data2.V=V2;
% assert(T==(info.max_idx+1)/2,'the total image size is not matched!');


data1.Yr = reshape(data1.Y,prod(info.sz),T1);
data2.Yr = reshape(data2.Y,prod(info.sz),T2);

data1.nY = min(reshape(data1.Yr,prod(info.sz)*T1,1));
data2.nY = min(reshape(data2.Yr,prod(info.sz)*T2,1));



if info.config.magnification~=2 && info.config.magnification~=1 
    magnification = info.config.magnification*2/5;
else
    magnification = info.config.magnification;
end

data1.magnification = magnification;
data2.magnification = magnification;
