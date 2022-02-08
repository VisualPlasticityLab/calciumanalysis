function peak=Fit2Gaussian(pref,xR,yR)

% s = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[0,0],...
%                'Upper',[Inf,max(cdate)],...
%                'Startpoint',[1 1]);
% f = fittype('a*(x-b)^n','problem','n','options',s);
% Fit the data using the fit options and a value of n = 2:
% [c2,gof2] = fit(cdate,pop,f,'problem',2)


%% prepross data
x=mod(xR+(4-pref),16)+1;
y=yR;
ymin = min(y);
ymax = max(y);
[ymax2,pos]=max(y(x>=9));
center1 = 5;
center2 =x(pos);

% xdata=1:.1:16;
% ydata=interp1(x,y,xdata,'pchip');
% 
% [yupper,ylower]=envelope(xdata,ydata,'pchip');
% ydata=yupper;
% ydata2= (yupper+ylower)/2;
% figure;plot(xdata,ydata,'o',xdata,yupper,xdata,ylower,xdata,ydata2)
%% fit two-peak guassaian
s = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[   0,  0,  1, 9,0,0,-ymax/2],...
               'Upper',[ymax,ymax, 9,16,8,8,ymax],...
               'Startpoint',[y(5)-ymin ymax2-ymin center1 center2 2 2 ymin]);
f = fittype('a1*exp(-((x-mod(b1,16)-1)/c1)^2)+a2*exp(-((x-mod(b2,16)-1)/c2)^2)+d','options',s);
[fit1,gof,fitinfo] = fit(x',y',f);
%[fit1,gof,fitinfo] = fit(xdata',ydata',f);
fit1
% plot(x,y,'o',xdata,ydata,'--',x,fit1(x),'-');
plot(x,y,'o--',x,fit1(x),'-');
axis tight;
%axis off
title(gof.rsquare);
%title(gof.rmse)
peak=fit1(center1);