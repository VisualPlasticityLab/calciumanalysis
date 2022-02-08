function [xpos,ypos,xwidth,yheight]=figurepara(xcol,yrow,xmin,xmax,ymin,ymax)
margin=0.05;
xwidth=subsiz(margin,xcol);
yheight=subsiz(margin,yrow);
xpos=subpos(margin,xcol);
ypos=subpos(margin,yrow);

if nargin>2
% xscale=(xmax-xmin)*xcol/50;
% yscale=(ymax-ymin)*yrow/2;
% subplot('position',[xpos(1) ypos(1) xwidth yheight]);
% plot([xmin;xscale+xmin],[ymax;ymax],'-k','Linewidth',2);
% plot([xmin;xmin],[ymax-yscale;ymax],'-k','Linewidth',2);
% text(xmin,ymax-yscale, num2str(yscale), 'HorizontalAlignment','left','VerticalAlignment','top');
% text((xmin+xscale)/2,ymax, num2str(xscale), 'HorizontalAlignment','center','VerticalAlignment','top');
end

function pos=subpos(margin,div)
pos=margin+(1-margin*2)/div*((1:div)-1);

function siz=subsiz(margin,div)
siz=(1-margin*2)/div*.75;





