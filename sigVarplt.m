function h=sigVarplt(sigF,Var,win_sig,matrix,Idx,folder)
%sigF seg*rep*Var*ncell
%maxtrix rep*Var

if prod(Var)~=size(sigF,3)
    disp('resize not matched, retry');return
end

seg=size(sigF,1);
rep=size(sigF,2);
ncell=numel(Idx);

if isempty(matrix)
    matrix=logical(zeros(rep,prod(Var)));
end

cellperpage=1;
page=ceil(ncell/cellperpage);
for j=1:page
    num=min(cellperpage,ncell-cellperpage*(j-1)); %plot 10 or mod(ncell,10) cells per page
%     h(j)=figure('Name',['reorganized_ page#' num2str(j)], 'position',[ 200 100 100*Var(1) 100*num]);
    h(j)=figure('Name',['reorganized_ page#' num2str(Idx(j))], 'position',[ 200 100 100*Var(1) 100*num]);
    [xpos,ypos,xwidth,yheight]=figurepara(Var(1),num,1,seg);
    for i=1:num    %%%%%cell#=i+10*(j-1)
            nth=Idx(i+cellperpage*(j-1));
            ymax=prctile(reshape(sigF(:,:,:,nth),[],1),99);
%             ymax = min(ymax,1000);
            subplot('position',[xpos(1) ypos(i) xwidth yheight/Var(2)]);
            text(0,0,['cell#' num2str(nth)],'HorizontalAlignment','right','VerticalAlignment','bottom');
            
        for k=1:Var(1)   %%%ori
            for m=1:Var(2)   %%%%%var#=m+Var(2)*(k-1)
                kth=m+(k-1)*Var(2);
                try
                    tempR=sigF(:,matrix(:,kth),kth,nth);        %tempR  seg,some rep,1,1;
                    rhoR = corr(tempR);
                    cor(2)=mean(rhoR(~triu(ones(size(rhoR,1)))));
                end
                try
                    tempS=sigF(:,~matrix(:,kth),kth,nth);
                    rhoS = corr(tempS);
                    cor(1)=mean(rhoS(~triu(ones(size(rhoS,1)))));
                    
                end
                subplot('position',[xpos(k) ypos(i)+yheight/Var(2)*(m-1) xwidth yheight/Var(2)]);hold on;
                plot(tempR,'Color',[1 0 0 .3],'linewidth',.5);
                plot(tempS,'Color',[0 0 1 .3],'linewidth',.5);
                plot(1:seg,mean(tempR,2),'r','linewidth',1);
                plot(1:seg,mean(tempS,2),'b','linewidth',1);
                plot([win_sig(1) win_sig(1)] ,[0 ymax*.8] ,[win_sig(end) win_sig(end)], [0 ymax],'--','linewidth',1);% stimulus on time
                axis([1 seg 0 ymax]);
                a=gca;
                set(a,'xticklabel',[]);
                axis off;
                text(0,ymax*.9,sprintf('S%.2fR%.2f',cor),'HorizontalAlignment','left','Fontsize',6);
                
            end
            text(seg*.9,ymax*.9,[num2str(k)],'VerticalAlignment','bottom','Fontsize',6);
        end
        scalebar('unit',{'frame','df/F'})
%        answer = inputdlg('Enter gd cells','Pick gd cells',[1 80],{num2str(Cor(nth))});
%        Gd=[Gd str2num(answer{1})];

    end
     drawnow;
    saveas(h(j),fullfile(folder,h(j).Name),'png');
     close(h(j))
end
