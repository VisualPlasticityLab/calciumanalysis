function peakRSplt(peakR,peakS);

Var=size(peakR,1);
ncell=size(peakR,2);


Ori=Var;  %
x=[1:(Ori+1)]';
xp=[0:2*pi/Ori:2*pi]';
xx=[1:.1:(Ori+1)]';
yR=[peakR(:,:); peakR(1,:)];
yS=[peakS(:,:); peakS(1,:)];

for j=1:ceil(ncell/10)
    num=min(10,ncell-10*(j-1)); %plot 10 or mod(ncell,10) cells per page
    figure('Name',['reorganized_ page#' num2str(j)], 'position',[ 200 100 1600 100*num]);
    [xpos,ypos,xwidth,yheight]=figurepos(2,num);
    
    for i=1:num    %%%%%cell#=i+10*(j-1)
        try
        ymax=max(peakR(:,i+10*(j-1)));
        catch
        ymax=max(peakS(:,i+10*(j-1)));
        end

        subplot('position',[ ypos(i) xpos(1)  yheight xwidth]);
        plot(x,yR(:,i+10*(j-1)),'r*--',x,yS(:,i+10*(j-1)),'b*-');      
        ylabel(['cell' num2str(i)])
        axis([ 1 Ori+1 0 ymax])
        axis off
        
        subplot('position',[ ypos(i) xpos(2)  yheight xwidth]);
        P = polar(xp, ymax * ones(Ori+1,1));
        set(P, 'Visible', 'off');
        hold on
        polar(xp,yR(:,i+10*(j-1)),'--r');
        polar(xp,yS(:,i+10*(j-1)),'-b');
        axis off
    end
    
end