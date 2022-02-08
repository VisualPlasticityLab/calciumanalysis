function newnames=makememmapnplane(fn)
%% define the name
sbxread(fn,1,1);% read one frame to read the header of the image sequence
global info;

path=pwd;
[~,remain] = strtok(path,'2'); %2pdata
zfs_path=fullfile('\\mps-zfs\data1\jsun',remain)
if ~exist(zfs_path,'dir')
    mkdir(zfs_path);
end

nplanes=info.otparam(3);
for i=1:nplanes
    newnames{i}=fullfile(zfs_path,[fn '_' num2str(i) '_memmap.mat']);
    data{i} = matfile(newnames{i},'Writable',true);
    disp(sprintf('Saving memmapfiles for %d planes:%s',nplanes,newnames{i}));
end

new = 1;

% if exist(newnames{1},'file')
%     %     new = menu('memmap file already exist!','Rewrite','Continue','exit');
%     %     if new == 3
%     return;
%     %     end
% end
%% save each piece
k=1;
tic
j0=0;
if new == 1
    nY = Inf;
    fac=sqrt(info.max_idx);
    for i =1:nplanes
        try
            m(:,:,i)=info.aligned.mr(:,:,k);
        catch
            m(:,:,i)=info.aligned.m(:,:,k);
        end
        try
            mg(:,:,i)=info.aligned.mg(:,:,k);
        catch
            mg(:,:,i)=info.aligned.m(:,:,k);
        end
        V{i}=zeros(info.sz);
        T{i}=numel(i:nplanes:info.max_idx+1);
        data{i}.Y=zeros([info.sz 2],'uint16');
        data{i}.m = m{i};
        data{i}.sizY =[info.sz T{i}];
    end
else
    for i =1:nplanes
        try
            m(:,:,i)=info.aligned.mr(:,:,k);
        catch
            m(:,:,i)=info.aligned.m(:,:,k);
        end
        try
            mg(:,:,i)=info.aligned.mg(:,:,k);
        catch
            mg(:,:,i)=info.aligned.m(:,:,k);
        end
        
        V{i}=zeros(info.sz);
        T{i}=numel(i:nplanes:info.max_idx+1);
        temp = size(data{i},'Y');
        j0 = j0+temp(end);
    end
end

for j=j0:info.max_idx
    z =sbxread(fn,j,1);
    z= squeeze(z(1,:,:));
    i = mod(j,nplanes)+1;
    img = circshift(z,info.aligned.T(j+1,:)); % align the image
    data{i}.Y(:,:,ceil((j+1)/nplanes)) = img;
    V{i}=V{i}+((double(img-m{i}))/fac).^2;
    if mod(j,500)==0
        fprintf('Saved Frame %d/%d for %.2f seconds\n ',j,info.max_idx,toc);
    end
end

if info.config.magnification~=2 && info.config.magnification~=1
    magnification = info.config.magnification*2/5;
else
    magnification = info.config.magnification;
end

for i =1:nplanes
    data{i}.magnification = magnification;
    data{i}.V=V{i}/max(max(V{i}));
    data{i}.Yr = reshape(data{i}.Y,prod(info.sz),T{i});
    data{i}.nY = min(reshape(data{i}.Yr,prod(info.sz)*T{i},1));
end
toc
