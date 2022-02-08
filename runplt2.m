function [matrix,h,run_thr]=runplt2(run,win_sig) 
% run:rep,Var,seg(prestim+stimON+stimOFF-prestim),
% run_adjusted:rep,Var,seg(prestim+stimON+stimOFF-prestim),
global info;
rep = size(run,1);
Var = size(run,2);
seg = size(run,3);
run_adjusted=run(:,:,win_sig(1):win_sig(2));
stimON = size(run_adjusted,3); 

%% define threshold for running/ still status
h(1)=figure('Name','running speed','position',[ 200 200 600 600]);
%run_diff=max(diff(run,1,3),[],3); % run_diff:rep,Var
% run_max=squeeze(max(run,[],3));
subplot(3,1,1);
% %plot(1:rep,run_max','*');
% histogram(run_max',0:.25:10);
% title('max run')
% legend('show')

maxspeed = min(prctile(run(:),95),10);
if maxspeed >0
histogram(run(:),0:maxspeed/50:maxspeed);
xlabel('speed')
drawnow;

subplot(3,1,2)
run_avg=squeeze(mean(run_adjusted,3));% run_avg:rep,Var
run_var = squeeze(var(run_adjusted,0,3));
%plot(1:rep,run_avg','*');
histogram(run_avg,0:maxspeed/25:maxspeed);
title('average adjusted speed')
xlabel('speed')

% cutoff1=inputdlg('define states','',1,{'1'});
% acc_thr=str2num(cutoff1{1});
% subplot(1,3,1);hold on;
% plot([ 1 rep],[acc_thr acc_thr],'--');
% cutoff2=inputdlg('define states','',1,{'.5'});
% run_thr=str2num(cutoff2{1});
run_thr=1.5;
subplot(3,1,2);hold on;
plot([run_thr run_thr],[ 1 rep],'--');
drawnow;

% subplot(3,2,5)
% ACC=run_max>acc_thr;
% imagesc(ACC);

subplot(3,1,3);
RUN=run_avg>run_thr;
imagesc(RUN);
%matrix=RUN&ACC;  %%%%matrix=rep*stim    %%%get two methonds, pick one to determine different states
matrix=RUN;
else
    matrix=[];
    run_thr=0;
end

%% to show the running/still status for each recorded segment
if run_thr>0
    h(2)=figure('Name','running segment','position',[ 200 200 1200 200]);

if length(info.steps)==1
    info.steps(2)=1;
    info.steps(1)=Var;
end
[xpos,ypos,xwidth,yheight]=figurepara(info.steps(1),info.steps(2));
ymax=prctile(run(:),99);
for i=1:info.steps(1)
    
    for j=1:info.steps(2)
        nth=j+(i-1)*info.steps(2);
        subplot('Position',[xpos(i) ypos(j) xwidth yheight]);hold on      
        plot(squeeze(run(matrix(:,nth),nth,:))','r');
        plot(squeeze(run(~matrix(:,nth),nth,:))','b');
        plot([win_sig(1) win_sig(1)],[0 ymax],'k--');
        plot([win_sig(2) win_sig(2)],[0 ymax],'k-');
        xlim([1 seg]);
        ylim([0 ymax]);
        axis off;
        title(['ori' num2str(i) 'ctr' num2str(j)])
        if i==1&j==1
            line([1. 1],[0 run_thr],'Linewidth',3,'color','r')
            text(1,run_thr,'1.5cm/s','FontSize',8,'HorizontalAlignment','Right')
        end
    end
end

% if any(all(matrix)==1) || any(all(~matrix)==1)    
%     choice = questdlg([' status are not covered for every type of stimulus' num2str(find(all(matrix)==1)) num2str(find(all(~matrix)==1))],'','no running','ignore','rerun','rerun' );    
%     switch choice
%         case 'rerun'
%             [matrix,h,run_thr]=runplt2(run,win_sig);%runplt(run,run_adjusted); 
%         case 'no running'
%             matrix=[];
%             run_thr=0;
%         case 'ignore'
%     end
% end

end
