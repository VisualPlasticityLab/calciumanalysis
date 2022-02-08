function h=oricolormap(x,y,mag,ori,cond1,cond2)

if ~exist('x')
    x=0;y=0;
    h=figure;hold on
    axis([-1 1 -1 1])
end
if ~exist('mag')
    mag =1;
end

if ~exist('ori')
    ori=4;
end

de=0:180/ori:360;
colortype = ori;
h = (1:colortype)'/colortype;
s = ones(colortype,1);
v = ones(colortype,1);
clut =squeeze(hsv2rgb(h,s,v));
clut(ori+1:ori*2,:) = clut/2;

%clut(:,:)= .5;


for i=1:numel(de)-1 % plot 0degree not 360degree
arrow([x y],[x+cosd(de(i))*mag y+sind(de(i))*mag],30,'Color',clut(i,:));
text(x+cosd(de(i))*mag,y+sind(de(i))*mag,num2str(de(i)),'HorizontalAlignment','center');
end
scatter(x,y,mag*30,[1 1 1]/2,'filled')
scatter(x,y,mag*100,[1 1 1]/2)
text(x,y,cond1,'HorizontalAlignment','center')
text(x+mag,y,cond2)

