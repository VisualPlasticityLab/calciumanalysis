function correlationgraph(M,figtitle,L)
n=size(M,1);
imagesc(M,[-.2 1]); % plot the matrix
set(gca, 'XTick', 1:round(sqrt(n)):n); % center x-axis ticks on bins
set(gca, 'YTick', 1:round(sqrt(n)):n); % center y-axis ticks on bins
colormap('jet'); % set the colorscheme
%colorbar
title(figtitle, 'FontSize', 14); % set title
