function [a1,ori,peak2]=fit2peakgaussion(pref0,peak0)

% s = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[0,0],...
%                'Upper',[Inf,max(cdate)],...
%                'Startpoint',[1 1]);
% f = fittype('a*(x-b)^n','problem','n','options',s);
% Fit the data using the fit options and a value of n = 2:
% [c2,gof2] = fit(cdate,pop,f,'problem',2)


%% prepross data
x0=1:16/numel(peak0):16;
x=1:16;
peak = interp1( x0,peak0,x,'pchip','extrap');
pref = (pref0-1) /numel(peak0)*16+1;
y=circshift(peak,[1,5-pref]);
ymin = min(peak);
[ymax2,ypos2]=max(peak(9:16));
ypos2=ypos2+8;

% xdata=1:.1:16;
% ydata=interp1(x,y,xdata,'pchip');
% 
% [yupper,ylower]=envelope(xdata,ydata,'pchip');
% ydata=yupper;
% ydata2= (yupper+ylower)/2;
% figure;plot(xdata,ydata,'o',xdata,yupper,xdata,ylower,xdata,ydata2)
%% fit two-peak guassaian
s = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[   0,  0,  1, 9,0,0,-y(5)],...
               'Upper',[y(5),ymax2,9,16,8,8,y(5)],...
               'Startpoint',[y(5)-ymin ymax2-ymin 5 ypos2 2 2 ymin]);
f = fittype('a1*exp(-((x-b1)/c1)^2)+a2*exp(-((x-b2)/c2)^2)+d','options',s);
[fit1,gof,fitinfo] = fit(x',y',f);
peak1=circshift(fit1(x)',[1,pref-5]);
peak2=peak1(x0);
% plot(x,y,'o',xdata,ydata,'--',x,fit1(x),'-');
%axis off;
%title(gof.rsquare);
%title(gof.rmse)
% plot(peak2,'-x','Linewidth',2);
[~,ori] =max(peak2);
% ori = mod(ori-1, floor(numel(x0)/2))+1;
a1=fit1.a1;

