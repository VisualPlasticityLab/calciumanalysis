function [matrix,hfig,run_thr]=runplt3(run,win_sig)
% run:rep,Var,seg(prestim+stimON+stimOFF-prestim),
% run_adjusted:rep,Var,seg(prestim+stimON+stimOFF-prestim),
% win_sig [prestim+1 prestim+stimON seg];

global info;
rep = size(run,1);
Var = size(run,2);
seg = size(run,3);

run_adjusted=run(:,:,win_sig(1):win_sig(2));

%% Define threshold for running still status

hfig(1)=figure('Name','running speed','position',[ 200 200 600 600]);

maxspeed = min(prctile(run(:),95),10);
if maxspeed >0
    subplot(3,1,1);
    histogram(run(:),0:maxspeed/50:maxspeed);
    xlabel('speed')
    drawnow;
    
    subplot(3,1,3);
    runavg=reshape(run_adjusted,rep*Var,[]);
    [idx,C] = kmeans(runavg,2);
    [run_thr,status] = sort(mean(C,2));
    for i=1:numel(run_thr)
        matrix(idx==status(i)) = i;
    end

    matrix = reshape(matrix,rep,Var);
    imagesc(matrix);

else
    matrix=[];
    run_thr=0;
end

%% to show the running/still status for each recorded segment
    hfig(2)=figure('Name','running segment','position',[ 200 200 1200 200]);
    if length(info.steps)==1
        info.steps(2)=1;
    end
    [xpos,ypos,xwidth,yheight]=figurepara(Var,info.steps(2));
    ymax=max(run(:));
    for i=1:Var
        for j=1:info.steps(2)
            nth=j+(i-1)*info.steps(2);
            subplot('Position',[xpos(i) ypos(j) xwidth yheight]);hold on
            for k=1:numel(run_thr)
                plot(squeeze(run(matrix(:,nth)==k,nth,:))','Color',hsv2rgb(k/numel(run_thr),1,1));
            end
            plot([win_sig(1) win_sig(1)],[0 ymax],'k--');
            plot([win_sig(2) win_sig(2)],[0 ymax],'k-');
            
            xlim([1 seg]);
            ylim([0 ymax]);
            %         axis off;
            title(['ori' num2str(i) 'ctr' num2str(j)])
        end
    end
    
    %% convert to logical run(1)/still(0) for consistency
        matrix = logical(matrix-1);

    % if any(all(matrix)==1) || any(all(~matrix)==1)
    %     choice = questdlg([' status are not covered for every type of stimulus' num2str(find(all(matrix)==1)) num2str(find(all(~matrix)==1))],'','no running','ignore','rerun','rerun' );
    %     switch choice
    %         case 'rerun'
    %             [matrix,h,run_thr]=runplt(run,run_adjusted);
    %         case 'no running'
    %             matrix=[];
    %             run_thr=0;
    %         case 'ignore'
    %     end
    % end
    
