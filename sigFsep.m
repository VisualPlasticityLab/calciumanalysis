function [sigR,errorR,sigS,errorS]= sigFcmp(sigF,matrix,prestim);
%sigR,sigS seg*Var*ncell; sigF  seg,rep,Var,ncell;  rep*Var
seg=size(sigF,1);
Var=size(sigF,3);
ncell=size(sigF,4);
page=ceil(ncell/10);
sigR=zeros(seg,Var,ncell);
sigS=sigR;

ymin=0;
ymax=2;


for k=1:Var
    tempR=sigF(:,matrix(:,k),k,:);        %tempR  seg,some rep,1,ncell
    tempS=sigF(:,~matrix(:,k),k,:);
    sigR(:,k,:)=squeeze(mean(tempR,2)); % sigR seg,1,1,ncell;
    errorR(:,k,:) = std(tempR,0,2);     % errorT seg,1,Var,ncell;
    sigS(:,k,:)=squeeze(mean(tempS,2));
    errorS(:,k,:) = std(tempS,0,2);
end

ymax=max(max(sigF(:,:,k,nth)));

plot(tempR,'k','linewidth',.5);
plot(tempS,'g','linewidth',.5);
errorbar(1:seg,sigR(:,k,nth),errorR(:,k,nth),'r','linewidth',1);
errorbar(1:seg,sigS(:,k,nth),errorS(:,k,nth),'b','linewidth',1);

plot([prestim prestim],[0 ymax],'--','linewidth',.5);% stimulus on time
axis([1 seg 0 ymax]);
a=gca;
set(a,'xticklabel',[]);
text(seg/2,0,['ori' num2str(k) 'cell#' num2str(i+10*(j-1))],'VerticalAlignment','top');
hold off;
%drawnow;
end
end
end