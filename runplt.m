function [matrix,h,run_thr]=runplt(run); % run:rep,Var,stimON
%%%%%% new tryout 07-28-16, to show histgram instead of the scatter points 
%%%%%% for the running/still trials and then select threshold
global info;
rep=size(run,1);
Var=size(run,2);

%% define threshold for running/ still status
h(1)=figure('Name','running speed','position',[ 0 0 900 800]);
run_max=squeeze(max(run,[],3));
subplot(3,1,1);
hist(run_max(:),20);
title('max run')
xlabel('speed cm/s')
drawnow;

subplot(3,1,2)
run_avg=squeeze(mean(run,3));% run_avg:rep,Var
hist(run_avg(:),20);
title('average speed')
xlabel('speed cm/s')
drawnow;

cutoff1=inputdlg('define states','',1,{'4'});
acc_thr=str2num(cutoff1{1});
subplot(3,1,1);hold on;
plot([acc_thr acc_thr],[ 1 rep],'--');
subplot(3,1,2);hold on;

cutoff2=inputdlg('define states','',1,{'1'});
run_thr=str2num(cutoff2{1});
plot([run_thr run_thr],[ 1 rep],'--');

subplot(3,2,5)
ACC=run_max>acc_thr;
imagesc(ACC);

subplot(3,2,6)
RUN=run_avg>run_thr;
imagesc(RUN);
matrix=RUN&ACC;  %%%%matrix=rep*stim=rep*    %%%get two methonds, pick one to determine different states
%% show results and decide if redo
h(2)=figure('Name','running segment','position',[ 200 200 1200 800]);
[xpos,ypos,xwidth,yheight]=figurepara(Var,info.steps(2));
ymax=max(run(:));
for i=1:Var
    for j=1:info.steps(2)
        nth=j+(i-1)*info.steps(2);
        subplot('Position',[xpos(i) ypos(j) xwidth yheight]);hold on
        
        plot(squeeze(run(matrix(:,nth),nth,:))','r');
        plot(squeeze(run(~matrix(:,nth),nth,:))','b');
        
        ylim([0 ymax]);
        title(['ori' num2str(i) 'ctr' num2str(j)])
    end
end

if any(all(matrix)==1) || any(all(~matrix)==1)
    
    choice = questdlg([' status are not covered for every type of stimulus' num2str(find(all(matrix)==1)) num2str(find(all(~matrix)==1))],'','no running','ignore','rerun','rerun' );
    
    switch choice
        case 'rerun'
            [matrix,h]=runplt(run);
        case 'no running'
            matrix=[];
        case 'ignore'
    end
end


