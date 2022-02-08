function z=sbxplotline(fname,idx)
if nargin<2
    sbxread(fname,0,1);
    global info;
    disp(info.max_idx);
    idx=1:info.max_idx;
end
z=[];
for(i=1:length(idx))
    y = sbxread(fname,idx(i),1);
%     if ~isempty(info.aligned)
%         y = circshift(y(:,:,i),info.aligned.T(idx(i)+1,:)); % align the image
%     end
if info.area_line %area or line scanmode
    z(:,i,:)=y(:,end/2,:);% all channels
else
    z=cat(2,z,y);
end
end
h=figure; %('Position', [100, 100, 2000, 2000])
temp=permute(z,[2 3 1]);
if size(temp,3)>1 & size(temp,3)<3
     temp(:,:,3)=0;
     temp(:,:,2)=0;
     
end
imshow(temp);
save([fname '_line.mat'],'z');