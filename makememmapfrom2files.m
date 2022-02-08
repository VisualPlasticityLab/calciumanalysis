% function newname=makememmapfrom2files
%% define the name

files = dir('*_memmap.mat');
files=files(end:-1:1);
newname='';
for i=1:numel(files)
    pos=strfind(files(i).name,'_memmap.mat');
    newname=[newname files(i).name(1:pos-1) '&'];
end
newname=[newname(1:end-1) '_memmap.mat'];

if exist(newname,'file')
    disp([newname 'already existing']);
    return;
end

data = matfile(newname,'Writable',true);
disp(sprintf(' try saving memmap file: %s',newname));

%% save each piece
T = 0;

for i=1:numel(files)
    fn= files(i).name;
    load(fn,'Y');
    eachsize(i)=size(Y,3);
    load(fn,'m');
    imgm(:,:,i)=m;
    tic
    if i==1
        data.Y=Y;
    else
        [u v] = fftalign(imgm(:,:,i),imgm(:,:,1));
        figure;
        imshowpair(circshift(imgm(:,:,i),[u v]),imgm(:,:,1));
        title(fn)
        for j=1:eachsize(i)
            Y(:,:,j)=circshift(Y(:,:,j),[u v]);
        end
        data.Y(:,:,T+1:T+eachsize(i))=Y;
    end
    T=T+eachsize(i);      
end
data.eachsize=eachsize;

load(fn,'magnification','nY');
data.magnification = magnification;
data.nY = nY;
data.sizY =size(data,'Y');
X = data.sizY(1,1);
Y = data.sizY(1,2);
A = X*Y;
data.Yr = reshape(data.Y,A,T);

