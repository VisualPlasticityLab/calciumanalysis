function h=sigONEplt(SI,sigF,matrix,win_sig,Cor)
% sigF seg*rep*Var*ncell
% maxtrix rep*Var

seg=size(sigF,1);
rep=size(sigF,2);
Var=size(sigF,3);
if ~exist('Cor')
    ncell = size(sigF,4);
    Cor=1:ncell;
else
    ncell = numel(Cor);
end


[b,a] = butter(12,.2,'low');           % IIR filter design
sigF = filtfilt(b,a,sigF);
baseline=nanmean(sigF(win_sig(1)-5:win_sig(1),:,:,:));
sigF=sigF-repmat(baseline,[seg 1 1 1 ]);
% Gd=[];

if isempty(matrix)
    matrix=logical(ones(rep,Var)); skip=1;
else
    skip=0
end
skip =1;

cellperpage=20;
page=ceil(ncell/cellperpage);
tic
for j=1:page
    num=min(cellperpage,ncell-cellperpage*(j-1)); %plot 10 or mod(ncell,10) cells per page
    
    h(j)=figure('Name',['reorganized_ page#' num2str(j)], 'position',[ 50 50 50*Var 50*num]);
    [xpos,ypos,xwidth,yheight]=figurepara(Var,num,1,seg);
    for i=1:num    %%%%%cell#=i+10*(j-1)
        nth=Cor(i+cellperpage*(j-1));
        ymax=max(reshape(sigF(:,:,:,nth),[],1))/1.5;  %let's just try with max for now
        ymin=min(reshape(sigF(:,:,:,nth),[],1))/5;
        subplot('position',[xpos(1) ypos(i) xwidth yheight]);
        for k=1:Var
            try
                subplot('position',[xpos(k) ypos(i) xwidth yheight]);hold on;
                tempR=sigF(:,matrix(:,k),k,nth);        %tempR  seg,some rep,1,1;
                plot(1:seg,tempR,'k','linewidth',.25);
                plot(1:seg,nanmean(tempR,2),'r','linewidth',2);
                if ~skip
                tempS=sigF(:,~matrix(:,k),k,nth);
                plot(1:seg,tempS,'g','linewidth',.25);
                plot(1:seg,nanmean(tempS,2),'b','linewidth',2);
                end
                %	k.color = [.25 .25 .25];   %grey
                %	g.color = [0 .25 0];
                %	patchline(1:seg,tempR,'edgecolor','k','linewidth',.5,'edgealpha',0.5)
                %	patchline(1:seg,tempS,'edgecolor','g','linewidth',.5,'edgealpha',0.5)
                % plot(1:seg,nanmean(tempS,2),'b','linewidth',2);
                %	ymax = max(max(max(nanmean(tempR,2)),max(nanmean(tempS,2))));
                plot([win_sig(1) win_sig(1)],[0 ymax],'k--',[win_sig(end) win_sig(end)], [0 ymax],'k--','linewidth',.25);% stimulus on time
                axis([win_sig(1)-5 seg ymin ymax]);
                %                 a=gca;
                %                 set(a,'xticklabel',[]);
                axis off;
                %peakvalues=sprintf('R=%.2f,S=%.2f',SI.peakR(:,k,nth),SI.peakS(:,k,nth));
                %text(0,0,peakvalues,'VerticalAlignment','top');
            end
        end
        %ODIvalues=sprintf('#%dR:%.2fS:%.2f',Cor(nth),SI.gOSI_R(nth),SI.gOSI_S(nth));
        %text(seg*1.05,ymax/2,ODIvalues,'HorizontalAlignment','left');
        
    end
    drawnow;
    %         answer = inputdlg('Enter gd cells','Pick gd cells',[1 80],{num2str(Cor(nth))});
    %         Gd=[Gd str2num(answer{1})];
end
disp(toc)
