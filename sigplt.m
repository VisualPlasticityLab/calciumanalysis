function hsig=sigplt(sig,stim,Cor,sigs,sigsp)

ncell=size(sig,2);
nframe=size(sig,1);
page=ceil(ncell/10);

if ~exist('stim')
    stim=[];
end
if ~exist('Cor')
    Cor=1:ncell;
end
if ~exist('sigsp')
    for i=1:ncell
%         [sigs(:,i),sigsp(:,i)]=deconvolveCa(sig(:,i));
    end
end


for j=1:ceil(ncell/10)
    num=min(10,ncell-10*(j-1)); %plot 10 or mod(ncell,10) cells per page            
    hsig(j)=figure('Name',['Ca sig w stim_ page#' num2str(j)], 'position',[50 50 200*num 50*num]);
    [xpos,ypos,xwidth,yheight]=figurepara(1,num);

    for k=1:num
        nth=k+10*(j-1);
        ymax=max(sig(:,nth));

        subplot('position',[xpos ypos(k) xwidth yheight]);hold on;
%         plot(1:nframe,sig(:,nth),'linewidth',1); 
        plot(1:nframe,sig(:,nth),'Linewidth',1,'Color',[1 0 0 0.3]); 
%         plot(1:nframe,sigs(:,nth),'Linewidth',1,'Color',[1 0 0 0.7]);                
%         plot(1:nframe,sigsp(:,nth),'b');        
        plot(1:nframe,stim'*ymax,'b','Linewidth',1);
        text(nframe*1.00,0,['Cell#' num2str(Cor(nth))],'HorizontalAlignment','left');
        axis tight
        %xlim([1 nframe]);
        subplot('position',[xpos+xwidth*1.05 ypos(k) .1 yheight]);hold on;
        try
%             histogram(sig(:,nth),[0:ymax/50:ymax],'Normalization','probability')
%             [~,pos]=findpeaks(sig(:,nth));
%             histogram(diff(pos))
Fr=pwrplt(sigs(:,nth));
        catch
            title('noise');
        end
    end
end