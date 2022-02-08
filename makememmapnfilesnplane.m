function newnames=makememmapnfilesnplane(fns)
%% define the path, save to MPS-ZFS\data\jsun directly
path=pwd;
% [~,remain] = strtok(path,'1'); %data
% zfs_path=fullfile('\\mps-zfs\data1\jsun\',remain)
% 
% if ~exist(zfs_path,'dir')
%     mkdir(zfs_path);
% end

nfiles=numel(fns);
newname='';
for i=1:nfiles
    newname=[newname fns{i} '&'];
end
newname=newname(1:end-1);
sbxread(fns{1},1,1);
global info;% read one frame to read the header of the image sequence

if info.volscan==1
    nplanes=info.otparam(3);
    realplane = info.otparam(3);
else
    nplanes=1;
    realplane =1;
end

for k=1:realplane
    newnames{k}=fullfile(path,sprintf('%s_%d_memmap.mat',newname,k));
00000000000000000000    data{k} = matfile(newnames{k},'Writable',true);
    disp(sprintf('Saving memmapfiles for %d planes:%s',nplanes,newnames{k}));
end

if exist(newnames{k},'file')
    disp([newnames{k} ' memmap file already existing']);
    return;
end

setpref('Internet','SMTP_Server','keck.ucsf.edu')
setpref('Internet','E_mail','jsun@phy.ucsf.edu')

%% predefine total length
% this contains the information about the structure of the image
%calculate the total size of the image
eachsize=zeros(nplanes,nfiles);
for i=1:nfiles
    sbxread(fns{i},1,1);
   for  k=1:realplane
        eachsize(k,i)=numel(k:nplanes:info.max_idx+1);
    end
end
T0= sum(eachsize,2);

for  k=1:realplane
    data{k}.eachsize = eachsize(k,:);
    data{k}.sizY =[info.sz T0(k)];
    data{k}.magnification = info.config.magnification/2.5;
end
%% align different files
for k=1:realplane  % for k=1:nplanes
    % First align all files
    for i=1:nfiles
        sbxread(fns{i},1,1);% read one frame to read the header of the image sequence
        
        % for single files wasn't loading alignment data -- MCD
        if(exist([fns{i} ,'.align'])) % aligned?
            info.aligned = load([fns{i} ,'.align'],'-mat');
        else
            info.aligned = [];
        end   
        
%         try
%             m(:,:,i)=info.aligned.mr(:,:,k);
%         catch
            m(:,:,i)=info.aligned.m(:,:,k);
%         end
        if i==1
            u=0; v=0;
        else
%             [u1(i) v1(i)] = fftalign(m(Center(1)-150:Center(1)+150,Center(2)-250:Center(2)+250,i),...
%                 m(Center(1)-150:Center(1)+150,Center(2)-250:Center(2)+250,1));
%             h0 = figure('Name','fftalign')
%             imshowpair(circshift(squeeze(m(:,:,i)),[u1(i) v1(i)]),(m(:,:,1)))
%             title(sprintf('u=%d,v=%d',u1(i),v1(i)));
%             [~,pos]=max(reshape(m(:,:,1),prod(info.sz),1));
%             [x,y]=ind2sub(info.sz,pos);
%             moving= m(x-100:x+100,y-100:y+100,i);
%             fixed = m(x-100:x+100,y-100:y+100,1);
%             [moving_reg,R_reg] =  imregister(moving, fixed, 'rigid', optimizer, metric);
            CenterL = info.sz/3;
            CenterR = info.sz/3*2;
            [~,I] = sort(reshape(m(:,:,1),prod(info.sz),1),'descend');
            disttoCenterL = (mod(I,info.sz(1))-CenterL(1)).^2+(ceil(I/info.sz(1))-CenterL(2)).^2;
            [~, pos1] = min(disttoCenterL(1:1000));
            disttoCenterR = (mod(I,info.sz(1))-CenterR(1)).^2+(ceil(I/info.sz(1))-CenterR(2)).^2;
            [~, pos2] = min(disttoCenterR(1:1000));
%             Area = mod(abs(I-I(pos1)),info.sz(1)).* ceil(abs(I-I(pos1))/info.sz(1));
%             pos2 = find(abs(mod(abs(I-I(pos1)),info.sz(1)))>200 ...
%                 & abs(mod(abs(I-I(pos1)),info.sz(1)))<300 ...
%                 & ceil(abs(I-I(pos1))/info.sz(1))>300 ...
%                 & ceil(abs(I-I(pos1))/info.sz(1))<500,1);
            [x1,y1] = ind2sub(info.sz,I(pos1));
            [x2,y2] = ind2sub(info.sz,I(pos2));
            moving = m(min(x1,x2):max(x1,x2),min(y1,y2)-10:max(y1,y2)+10,i);
            fixed = m(min(x1,x2):max(x1,x2),min(y1,y2)-10:max(y1,y2)+10,1);
            % fftalign
%             [u1(i) v1(i)] = fftalign(moving,fixed);
%              [u1(i) v1(i)]=fftalign(m(50:end-50,100:end-100,i),m(50:end-50,100:end-100,1))
            h0 = figure('Name','fftalign')
%              [u1(i) v1(i)]=fftalign(m(100:450,100:500,i),m(100:450,100:500,1))
%             imshowpair(circshift(squeeze(m(:,:,i)),[u1(i) v1(i)]),(m(:,:,1)))
%             title(sprintf('u=%d,v=%d',u1(i),v1(i)));

            
            while menu('Select areas for image alignment?','Yes','Done')==1
                align1img = squeeze(m(:,:,i));
                align2img = squeeze(m(:,:,1));
                imshowpair(align1img,align2img);
                rect=round(getrect(h0));
                img1=align1img(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3));
                img2=align2img(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3));
                [u1(i) v1(i)] = fftalign(img1,img2);
                align1imgnew = circshift(align1img,[u1(i) v1(i)]);
                imshowpair(align1imgnew,align2img);
                title(sprintf('u=%d,v=%d',u1(i),v1(i)))
            end
            figname= sprintf('plane%d aligned file#%d',k,i);
            saveas(h0,[figname '_fft.fig']);
            close(h0)
            u(i) = u1(i);v(i) = v1(i);
            % multimodal
%             [optimizer, metric] = imregconfig('multimodal');
%             optimizer.MaximumIterations = 1000;
%             tform = imregtform(moving, fixed, 'translation', optimizer, metric);
%             v(i)=round(tform.T(3,1));
%             u(i)=round(tform.T(3,2));
%             h1 = figure('Name','imregtform');
%             imshowpair(circshift(squeeze(m(:,:,i)),[u(i) v(i)]),(m(:,:,1)))            
%             title(sprintf('u=%d,v=%d',u(i),v(i)))
%             figname= sprintf('plane%d aligned file#%d',k,i);
%             if (u1(i)-u(i))^2+(v1(i)-v(i))^2>4
%                 sendmail({'j.suninchina@gmail.com','mdadarla@phy.ucsf.edu'},'Check alignment', ...
%                     [newnames{k} 'needs attention']);
%                 if menu('Select which alignment method?','fft','multimodal')==1
%                     u(i) = u1(i);v(i) = v1(i);saveas(h0,[figname '_fft.fig']);
%                 end
%             else
%                 saveas(h1,[figname '_imreg.fig']);
%             end
        end
    end
    data{k}.uv=[u;v];
    data{k}.nY = 3;
end

%% save each 
for k=1:realplane
    T=0;
    V=zeros(info.sz);
    data{k}.m=zeros(info.sz);
    data{k}.V=zeros(info.sz);
    data{k}.Y=zeros([info.sz 2],'uint16');
    uv = data{k}.uv; 
    u = uv(1,:); 
    v = uv(2,:);
    for i=1:nfiles
        tic
        sbxread(fns{i},1,1);% read one frame to read the header of the image sequence
        try
            mg=info.aligned.mg(:,:,k);
        catch
            mg=info.aligned.m(:,:,k);
        end
        mg = circshift(squeeze(mg),[u(i) v(i)]);
        fac=double(max(mg(:)));
        for j=k:nplanes:info.max_idx+1
            z =sbxread(fns{i},j-1,1);
            z= squeeze(z(1,:,:));
            nth=ceil(j/nplanes);
%             img(:,:,nth) = circshift(z,info.aligned.T(j,:)+[u v]); % align the image
            %             temp = min(img(:));
            %             if temp < nY
            %                 nY = temp;
            %             end
            img = circshift(z,info.aligned.T(j,:)+[u(i) v(i)]); % align the image
            V=V+abs(double(img-mg))/fac;
            data{k}.Y(1:info.sz(1),1:info.sz(2),nth+T)=img;
            if mod(nth,500)==0
                fprintf('alinged Plane%d File%d Frame%d in %.1f seconds\n',k,i,nth,toc);
            end
        end
%         data{k}.Y(1:info.sz(1),1:info.sz(2),T+1:T+eachsize(k,i)) = img;
        data{k}.Yr(1:prod(info.sz),T+1:T+eachsize(k,i)) = reshape(...
        data{k}.Y(1:info.sz(1),1:info.sz(2),T+1:T+eachsize(k,i)),...
        prod(info.sz),eachsize(k,i));
        ratio=T/(T+eachsize(k,i));
        data{k}.V=data{k}.V*ratio+V*(1-ratio);
        data{k}.m=data{k}.m*ratio+double(mg)*(1-ratio);
        T=T+eachsize(k,i);
        fprintf('reshaped Plane%d File%d Frame%d in %.1f seconds\n',k,i,nth,toc);
    end

end
