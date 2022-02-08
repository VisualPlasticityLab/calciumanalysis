
function newname=makememmapnfiles
%% define the name and save to MPS-ZFS directly
path=pwd;
if strcmpi(path(1),'C')==1
    zfs_path=strrep(path,path(1:2),'\\mps-zfs\data1\idl');
else
    zfs_path=strrep(path,'\\mps-pc53','\\mps-zfs\data1\idl');
end
files = dir('*.sbx');
newname=[files(1).name(1:end-5) 'x_memmap.mat'];
newname=fullfile(zfs_path,newname);

if exist(newname,'file')
    disp([newname ' memmap file already existing']);return;
elseif ~exist(zfs_path,'dir')
    mkdir(zfs_path);
end
data = matfile(newname,'Writable',true);
disp(sprintf(' try saving memmap file: %s',newname));

%% predefine total length
global info;            % this contains the information about the structure of the image
%calculate the total size of the image
for i=1:numel(files)
    fn=files(i).name(1:end-4);
    sbxread(fn,1,1);% read one frame to read the header of the image sequence
    eachsize(i)=info.max_idx+1;
end
T0= sum(eachsize);
data.eachsize=eachsize;
data.sizY =[info.sz T0];

%% save each piece
T=0;
for i=1:numel(files)
    fn=files(i).name(1:end-4);
    sbxread(fn,1,1);% read one frame to read the header of the image sequence
    m(:,:,i)=info.aligned.m;
    V=zeros(info.sz);
    fac=max(max(m(:,:,i)));
    tic
    if i==1
        data.m=zeros(info.sz);
        data.Y=zeros([info.sz 2],'uint16');
        data.V=zeros(info.sz);
        u=0;v=0;
    else
%         [u v] = fftalign(m(:,:,i),m(:,:,1));
            [optimizer, metric] = imregconfig('multimodal');
            moving= m(:,:,i);
            fixed = m(:,:,1);
            [moving_reg,R_reg] =  imregister(moving, fixed, 'rigid', optimizer, metric);
            tform = imregtform(moving, fixed, 'translation', optimizer, metric);
            v=round(tform.T(3,1));
            u=round(tform.T(3,2));

        figure('Name',['aligned file#' num2str(i)]);
        imshowpair(circshift(squeeze(m(:,:,i)),[u v]),squeeze(m(:,:,1)))
    end
    for j=0:info.max_idx
        z =sbxread(fn,j,1);
        z= squeeze(z);
        img = circshift(z,info.aligned.T(j+1,:)+[u v]); % align the image
        data.Y(:,:,T+j+1) = img;
        V=V+abs(double(img-m(:,:,i))/fac);
        if mod(j,1000)==0
            fprintf('File %d Frame %d/%d for %.2f seconds\n ',i,j,info.max_idx,toc);
        end
    end
    save([fn '.align'],'V','-append');
    ratio=T/(T+eachsize(i));
    data.V=data.V*ratio+circshift(V,[u v])*(1-ratio);
    data.m=data.m*ratio+circshift(double(m(:,:,i)),[u v])*(1-ratio);
    T=T+eachsize(i);
end
%assert(T==info.max_idx,'the total image size is not matched!')
data.Yr = reshape(data.Y,prod(info.sz),T0);
% data.nY = min(reshape(data.Yr,prod(info.sz)*T0,1));
data.nY =3;
if (info.config.magnification>=5)
    data.magnification = info.config.magnification/2.5;
else
    data.magnification = info.config.magnification;
end