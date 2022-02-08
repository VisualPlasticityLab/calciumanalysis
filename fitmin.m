function [Fred_corrected,p] = fitmin(Fgreen,Fred,debug)
if size(Fgreen,2)>1
    SIZ=size(Fgreen);
    Fgreen = Fgreen(:);
    Fred =Fred(:);
end
step = floor(sum(Fgreen<=prctile(Fgreen,99.9))/100);
xf=linspace(prctile(Fgreen,0.1),prctile(Fgreen,99.9),step);
yf=nan(size(xf)-1);
for k=1:step-1
    try
        yf(k)=min(Fred(Fgreen>=xf(k)&Fgreen<=xf(k+1)));
    end
end
xf = xf(1:end-1);

xData = xf(~isnan(yf));
yData = yf(~isnan(yf));
xData = xData(yData~=0);
yData = yData(yData~=0);

% fitresult = fit( xData', yData', 'poly1' );
% ybaseline =fitresult(green_mean);
p = polyfit(xData',yData',1);
% red_baseline = polyval(p,green_mean);
Fred_corrected = Fred-polyval(p,Fgreen);

%% plot
if ~exist('debug','var')
debug =0;
end
if debug
    linex= [0,max(Fgreen)];
    hold on;
    scatter(xData,yData,'bo');
    scatter(Fgreen,Fred,'b.');
    scatter(Fgreen,Fred_corrected,'r.');
    plot(linex,polyval(p,linex),'k-');
    plot(mean(Fgreen),mean(Fred_corrected),'kx','MarkerSize',16*debug,'Linewidth',3);
    axis tight
end
%% reshape
if exist('SIZ','var')
    Fred_corrected =reshape(Fred_corrected,SIZ(1),SIZ(2));
end