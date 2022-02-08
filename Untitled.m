nSteps = size(eye1.finalvalueR,2);

h1=figure('Position',[100 200 1500 500],'Name',['alltraces_cell#' num2str(j)])
ncol =nSteps*2+2;
Cor =[1 2 4 5 8 9 10];
ncell=numel(Cor);
%%
for j= 1:ncell
    nth1= pair1(Cor(j));
    nth2= pair2(Cor(j));
    
    
    temp_1st=squeeze(eye1.sigF(nth1,:,:,:));
    temp_2nd=squeeze(eye2.sigF(nth2,:,:,:));
    ymax=max([temp_1st(:);temp_2nd(:)]);
    ymin=min([temp_1st(:);temp_2nd(:)]);
    buffer=(ymax-ymin)*.1;
    for ori=1:nSteps
        subplot(ncell,ncol,ori+ncol*(j-1)); hold on;
        plot(squeeze(temp_1st(ori,eye1.matrix(:,ori),:))','b','linewidth',1);
        plot(squeeze(mean(temp_1st(ori,eye1.matrix(:,ori),:))),'r','linewidth',3);
        axis([win_sig1(1)-5 win_sig1(end)+5 ymin-buffer ymax+buffer])
        axis off
        
        subplot(ncell,ncol,ori+nSteps+ncol*(j-1)); hold on;
        plot(squeeze(temp_1st(ori,eye1.matrix(:,ori),:))','b','linewidth',1);
        plot(squeeze(mean(temp_1st(ori,eye1.matrix(:,ori),:))),'r','linewidth',3);
        axis([win_sig1(1)-5 win_sig1(end)+5 ymin-buffer ymax+buffer])
        axis off
    end  
    
    subplot(ncell,ncol,nSteps*2+1+ncol*(j-1));    hold on;
    errorbar(1:nSteps,eye1.finalvalueR(nth1,:),eye1.stdofeachvalueR(nth1,:),'r*-');
    errorbar(1:nSteps,eye2.finalvalueR(nth2,:),eye2.stdofeachvalueR(nth2,:),'b+-');
    
    subplot(ncell,ncol,nSteps*2+2+ncol*(j-1));    hold on;
    polarplt(cat(1,eye1.finalvalueR(nth1,:),eye2.finalvalueR(nth2,:)));
end