function h=peakplt(peak);

Var=size(peak,1);
ncell=size(peak,2);

Ori=Var;  %
x=[1:(Ori+1)]';
xp=[0:2*pi/Ori:2*pi]';
xx=[1:.1:(Ori+1)]';
y=[peak(:,:); peak(1,:)];

for j=1:ceil(ncell/10)
    num=min(10,ncell-10*(j-1)); %plot 10 or mod(ncell,10) cells per page
    h(j)=figure('Name',['reorganized_ page#' num2str(j)], 'position',[ 200 100 1600 100*num]);
    [xpos,ypos,xwidth,yheight]=figurepos(2,num);
    
    for i=1:num    %%%%%cell#=i+10*(j-1)
        ymax=max(peak(:,i+10*(j-1)));
        subplot('position',[ ypos(i) xpos(1)  yheight xwidth]);
        plot(x,y(:,i+10*(j-1)),'--o');      
        ylabel(['C#' num2str(i+10*(j-1))])
        axis([ 1 Ori+1 0 ymax])
        axis off
        
        subplot('position',[ ypos(i) xpos(2)  yheight xwidth]);
        P = polar(xp, ymax * ones(Ori+1,1));
        set(P, 'Visible', 'off');
        hold on
        polar(xp,y(:,i+10*(j-1)),'--o');
        axis off
    end
    
end