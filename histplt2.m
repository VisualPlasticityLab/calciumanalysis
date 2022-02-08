function histplt2(data1,data2,sc)

hold on
if exist('sc','var')
    h1 = histogram(data1(:),'BinEdges', sc);
    h2 = histogram(data2(:),'BinEdges', sc);
else
    blim = prctile(cat(1,data1(:),data2(:)),[5 95]);
    h1 = histogram(data1(:),'BinEdges',blim(1):diff(blim)/10:blim(end));
    h2 = histogram(data2(:),'BinEdges',blim(1):diff(blim)/10:blim(end));
end

h1.FaceColor=[0 0 0.5]; %blue
h2.FaceColor=[0.5 0 0];
% legend('still','run','Location','north','Orientation','horizontal')
% legend('boxoff')
% set(gca,'XTickLabel','')
% set(gca,'XLim', [h2.BinEdges(1) h2.BinEdges(end)])
linkaxes(gca,'x')
% axis tight