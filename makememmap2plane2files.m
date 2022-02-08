function newnames=makememmap2plane2files
%% define the name and save to MPS-ZFS directly
path=pwd;
if strcmpi(path(1),'C')==1
    zfs_path=strrep(path,path(1:2),'\\mps-zfs\data\jsun');
else
    zfs_path=strrep(path,'\\mps-pc53','\\mps-zfs\data\jsun');
end
files = dir('*.sbx');
newnames{1}=fullfile(zfs_path,[files(1).name(1:end-5) 'x_1_memmap.mat']);
newnames{2}=fullfile(zfs_path,[files(1).name(1:end-5) 'x_2_memmap.mat']);
newnames{3}=fullfile(zfs_path,[files(1).name(1:end-5) 'x_3_memmap.mat']);

if exist(newname{1},'file')
   disp([' memmap file already existing']);return;
elseif ~exist(zfs_path,'dir')
    mkdir(zfs_path);
end

data1 = matfile(newnames{1},'Writable',true);
data2 = matfile(newnames{2},'Writable',true);
disp(sprintf('Saving memmapfiles for 3 planes:%s,%s',newnames{1},newnames{2}));

%% save each piece
T=0;
nY = Inf;
for i=1:numel(files)
    fn=files(i).name(1:end-4);
    sbxread(fn,1,1);% read one frame to read the header of the image sequence
    global info;            % this contains the information about the structure of the image
    m1(:,:,i)=info.aligned.m(:,:,1);
    m2(:,:,i)=info.aligned.m(:,:,2);
    m3(:,:,i)=info.aligned.m(:,:,3);
    
    V1=zeros(info.sz);
    V2=zeros(info.sz);
    V3=zeros(info.sz)
    fac=sqrt(info.max_idx+1);
    
    eachsize(i)=(info.max_idx+1)/2;
    tic
    if i==1
        data1.Y=zeros([info.sz (info.max_idx+1)/2],'uint16');
        data2.Y=zeros([info.sz (info.max_idx+1)/2],'uint16');
        data3.Y=zeros([info.sz (info.max_idx+1)/2],'uint16');
        data1.V=zeros(info.sz);
        data2.V=zeros(info.sz);
        data3.V=zeros(info.sz);

        data1.m=zeros(info.sz);
        data2.m=zeros(info.sz);
        data3.m=zeros(info.sz);

        u1=0;v1=0;
    else
        [u1 v1] = fftalign(m1(:,:,i),m1(:,:,1));
        [u2 v2] = fftalign(m2(:,:,i),m2(:,:,1));
        [u3 v3] = fftalign(m3(:,:,i),m3(:,:,1));

    end
   
    for j=0:info.max_idx
        z =sbxread(fn,j,1);
        z= squeeze(z);
        if mod(j,2)==0
            data1.Y(:,:,T+j/2+1) = circshift(z,info.aligned.T(j+1,:)+[u1 v1]); % align the image
            V1=V1+((double(data1.Y(:,:,i+1)-m1(:,:,i)))/fac).^2;

        else 
            data2.Y(:,:,T+(j+1)/2) = circshift(z,info.aligned.T(j+1,:)+[u2 v2]); % align the image
            V2=V2+((double(data2.Y(:,:,i+1)-m2(:,:,i)))/fac).^2;
        end
        
        if mod(j,500)==0
            fprintf('File %d Frame %d/%d for %.2f seconds\n ',i,j,info.max_idx,toc);
        end
    end
    
    V=cat(3,V1,V2);
    save([fn '.align'],'V','-append');
    
    ratio=T/(T+eachsize(i));
    data1.V=data1.V*ratio+V1*(1-ratio);
    data2.V=data2.V*ratio+V2*(1-ratio);
    data1.m=data1.m*ratio+circshift(double(m1(:,:,i)),[u1 v1])*(1-ratio);
    data2.m=data2.m*ratio+circshift(double(m2(:,:,i)),[u2 v2])*(1-ratio);
    T=T+eachsize(i);

end

data1.eachsize=eachsize;
data2.eachsize=eachsize;

data2.Yr = reshape(data2.Y,prod(info.sz),T);
data1.Yr = reshape(data1.Y,prod(info.sz),T);

data1.nY = min(reshape(data1.Yr,prod(info.sz)*T,1));
data2.nY = min(reshape(data2.Yr,prod(info.sz)*T,1));

data1.sizY =[info.sz T];
data2.sizY =[info.sz T];

if mod(info.config.magnification*2,5) == 0
    magnification = info.config.magnification*2/5;
else
    magnification = info.config.magnification;
end

data1.magnification = magnification;
data2.magnification = magnification;

