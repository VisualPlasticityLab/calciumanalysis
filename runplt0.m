function [matrix,h,run_thr]=runplt0(run,win_sig)
% run:numberoftrials,seg(prestim+stimON+stimOFF), or
% seg = size(run,3);
% run_adjusted:numberoftrials,seg(stimON),
global info;
if ndims(run)>2
    runsizes = size(run);
    run = reshape(run,prod(runsizes(1:end-1)),runsizes(end));
end
seg = size(run,2);
numberoftrials = size(run,1);
run_adjusted=run(:,win_sig(1):win_sig(2));
stimON = numel(win_sig);
%% define threshold for running/ still status

maxspeed = min(prctile(run(:),95),10);
if maxspeed >0
    h(1)=figure('Name','running speed','position',[ 200 200 600 600]);
    subplot(3,1,1);
    histogram(run(:),0:maxspeed/50:maxspeed);
    xlabel('speed')
    drawnow;
    run_avg = nanmean(run_adjusted);
    run_var = nanvar(run_adjusted);
    run_thr=1;
    RUN=run_avg>run_thr;
    RUN(run_var>run_avg) = RUN(run_var>run_avg)/2;
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
    ymax=max(run(:));
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
