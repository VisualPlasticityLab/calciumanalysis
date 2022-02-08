function [xpos,ypos,xwidth,yheight]=figurepos(xcol,yrow)
margin=0.1;
xwidth=subsiz(margin,xcol);
yheight=subsiz(margin,yrow);
xpos=subpos(margin,xcol);
ypos=subpos(margin,yrow);

function pos=subpos(margin,div)
pos=margin+(1-margin*2)/div*((1:div)-1);

function siz=subsiz(margin,div)
siz=(1-margin*2)/div*0.8;





