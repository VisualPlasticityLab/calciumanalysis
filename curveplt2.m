function curveplt2(data1,data2,ori)
%% define colorscheme
if ~exist('ori')
    ori=4;
end
colortype = ori;
h = (1:colortype)'/colortype;
s = ones(colortype,1);
v = ones(colortype,1);
clut =squeeze(hsv2rgb(h,s,v));
clut(ori+1:ori*2,:) = clut/2;
clut(end+1,:) = [1 1 1]/2;

%% calculate curve but avoid NaN data on either day to contanimate
d1 = data1;
d2 = data2;

for i=1:ori*2+1
    temp1 = find(d1 ==i);
    d1(temp1) = d1(temp1) +rand(size(temp1))/4-1/8;
    temp2 = find(d2 ==i);
    d2(temp2) = d2(temp2) +rand(size(temp2))/4-1/8;
end

x0 = -20:20;
y0 = (d2-d1)/2*atan(x0)/atan(20)+(d2+d1)/2*ones(size(x0));
for i=1:numel(d2)
	plot(x0,y0(i,:),'Color',[clut(data1(i),:) .5],'LineWidth',4);
end

% dim = (size(data1,1)>1)+1; %the dimension that's 1
% avg = nanmean(cat(dim,data1,data2),dim);
% x0 = -20:0;
% y0 = (avg-data1)/2*atan(x0)/atan(20)+(avg+data1)/2*ones(size(x0));
% x2 = 0:20;
% y2 = (data2-avg)/2*atan(x2)/atan(20)+(data2+avg)/2*ones(size(x0));
% for i=1:numel(data1)
% 	plot(x0,y0(i,:),'Color',[clut(data1(i),:) .3],'LineWidth',2);
%     plot(x2,y2(i,:),'Color',clut(data2(i),:),'LineWidth',2);
% end

h1.FaceColor=[0 0 0.5];
h2.FaceColor=[0.5 0 0];
% legend('still','run','Location','north','Orientation','horizontal')
% legend('boxoff')
% set(gca,'XTickLabel','')
% set(gca,'XLim', [h2.BinEdges(1) h2.BinEdges(end)])
linkaxes(gca,'x')
% axis tight