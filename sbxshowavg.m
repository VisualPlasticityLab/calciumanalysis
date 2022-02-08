function sbxshowavg(fname,idx,nbinning,channel)
%specify the numbers of frames for binniing
global info;
sbxread(fname,1,1);

if nargin<4
    channel =1;
end
if nargin<3
    nbinning=1;
end
if nargin<2 || isempty('idx')
    idx = 0:nbinning:info.max_idx;
end

h=figure

for(i=1:length(idx))
    y = sbxread(fname,idx(i),nbinning);
    y = mean(y,4); % make average from the nframes for binning
    if size(y,1)>1 && ~exist('channel','var')
        z = squeeze(y(2,:,:))*5;
        z(:,:,2)=y(1,:,:)*10;
        z(:,:,3)=0;
    else
        z = squeeze(y(channel,:,:));
    end
    if ~isempty(info.aligned)
        z = circshift(z,info.aligned.T(idx(i)+1,:)); % align the image
    end
    imshow(z,[]) %0 max(z(:))*1.25 display image from 5-95% range
    colormap gray
    if i==1
        temp=double(z(:,:,1));
    else
    temp=temp+double(z(:,:,1));
    end
    %imshow(z,'InitialMagnification',100);
    %colormap gray
    %imagesc(z); 
    title(num2str(idx(i)+1));
    pause(.05);
    drawnow;
    axis on;
end
figure;imagesc(temp);colormap gray