function h=histplt2(figH,data1,data2);

hist(figH,data1(:));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75);
hold on
hist(figH,data2(:));
h = findobj(gac,'Type','patch');
set(h,'facealpha',0.75);
legend('running','still');