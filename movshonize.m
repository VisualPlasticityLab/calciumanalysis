function [] = movshonize()
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%Pascal Wallisch
%02/25/2011
%This function does the Movshon procedure on a figure. Makes it more
%palatable to the man. 
box off
set(gca,'TickDir','out')
set(gca,'XScale', 'linear')
% set(gca,'fontangle','oblique')
set(gca,'fontsize',16)
set(gca,'ticklength',[0.01 0.025])
axis square

end

