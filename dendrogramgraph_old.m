function neworder = dendrogramgraph_old(M,N,figtitle1,figtitle2)
n=size(M,1);
Z = linkage(M);
step = round(sqrt(n));

figure('Position',[1000 400 1000 700])
subplot('Position',[0.1 .75 .8 .15])
dendrogram(Z,0);
g=gca;
neworder=str2num(g.XTickLabel);
set(g,'XTick',1:step:n)
set(g,'XTickLabel',1:step:n)

subplot('Position',[0.5 0.1 .35 .5])
imagesc(M(neworder,neworder)); % plot the matrix
set(gca, 'XTick', 1:step:n); % center x-axis ticks on bins
set(gca, 'YTick', 1:step:n); % center y-axis ticks on bins
% set(gca, 'XTickLabel',neworder(1:step:n));
% set(gca, 'YTickLabel',neworder(1:step:n));
colormap('jet'); % set the colorscheme
title(figtitle1, 'FontSize', 14); % set title

subplot('Position',[0.1 0.1 .35 .5])
imagesc(N(neworder,neworder)); % plot the matrix
set(gca, 'XTick', 1:step:n); % center x-axis ticks on bins
set(gca, 'YTick', 1:step:n); % center y-axis ticks on bins
colormap('jet'); % set the colorscheme
title(figtitle2, 'FontSize', 14); % set title