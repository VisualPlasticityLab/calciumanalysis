function newfig = overlayfigs(fig1,fig2,fig3)

normed1 = normfig(fig1);
normed2 = normfig(fig2);
if nargin<3
    normed3 = normed1;
else
    normed3 = normfig(fig3);
end

newfig(:,:,1) = normed1;
newfig(:,:,2) = normed2;
newfig(:,:,3) = normed3;



function normed = normfig(fig)
fig = double(fig);
figrange = prctile(fig(:),[.1 99.9]);
normed = (fig-figrange(1))/(figrange(2)-figrange(1));