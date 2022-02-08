function hh = histplt3(data1,data2,color)

if size(data1,1)>1
    data1 = data1';
end
if size(data2,1)>1
    data2 = data2';
end
if ~exist('color')
    color = [.5 0.5 0.5];
end
pairnum = numel(data1);

figure; 
subplot('Position',[.35 .1 .35 .8])
hold on;
plot([1;2]*ones(1,pairnum),cat(1,data1,data2),'color',[.8 .8 .8])
plot([1;2],cat(1,mean(data1),mean(data2)),'ko-')
xlim([0.5 2.5])
hh(2) = gca;
bins = hh(2).YTick(1):diff(hh(2).YTick([1 end]))/10:hh(2).YTick(end);
set(hh(2),'YLim',bins([1 end]))
set(hh(2),'XTick',1:2)
set(hh(2),'XTickLabel',{'Pre-','Post-'});
hh(2).YAxis.Visible = 'off';   % remove y-axis

subplot('Position',[.1 .1 .2 .8])
hist(data1,bins)
hh(1) = gca;
hpatch = findobj(gca,'Type','patch');
hpatch.FaceColor =color;
hh(1).XTick = bins(1:round(sqrt(numel(hh(2).YTick))):end);
box off
camroll(90)
% view([180 90])

subplot('Position',[.7 .1 .2 .8])
hist(data2,bins)
hpatch = findobj(gca,'Type','patch');
hpatch.FaceColor =color/2;
hh(3) = gca;
hh(3).XTick = bins(1:round(sqrt(numel(hh(2).YTick))):end);
set(hh(3),'XDir','reverse')   
hh(3).YAxisLocation = 'right';
box off
camroll(-90)
linkaxes(hh([1 3]))
set(gca,'XLim',bins([1 end]))
