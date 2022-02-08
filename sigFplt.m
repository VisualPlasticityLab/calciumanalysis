function h=sigFplt(sigF,matrix,win_sig,Cor,folder)
% sigF seg*rep*Var*ncell
% maxtrix rep*Var
%%
ncell = numel(Cor);
seg=size(sigF,1);
rep=size(sigF,2);
Var=size(sigF,3);  % Var=1;

if isempty(matrix)
    matrix=ones(rep,Var);
end
matrix = double(matrix);
if min(matrix(:))==0
    matrix = matrix+1;
end
status = numel(unique(matrix));
%% adjust baseline for calcium raw traces
baseline= nanmean(sigF(1:win_sig(1)-1,:,:,:),1);
sigF =sigF - repmat(baseline,[seg 1 1 1]);
%%
cellperpage=1;
page=ceil(ncell/cellperpage);
tic
for j=1:page
    num=min(cellperpage,ncell-cellperpage*(j-1)); %plot 10 or mod(ncell,10) cells per page
    h(j)=figure('Name',['reorganized_ page#' num2str(Cor(j))], 'position',[ 50 50 100*Var 200*num]);
    
%     h(j)=figure('Name',['reorganized_ page#' num2str(j)], 'position',[ 50 50 100*Var 200*num]);
    [xpos,ypos,xwidth,yheight]=figurepara(Var,num,1,seg);
    for i=1:num    %%%%%cell#=i+10*(j-1)
        nth=Cor(i+cellperpage*(j-1));
            
        ymax=max(reshape(sigF(:,:,:,nth),[],1));  %let's just try with max for now
        ymin=min(reshape(sigF(:,:,:,nth),[],1));  %let's just try with max for now
        subplot('position',[xpos(1) ypos(i) xwidth yheight]);
        
        for k=1:Var
                subplot('position',[xpos(k) ypos(i) xwidth yheight]);hold on;
                for n=status:-1:1
                    temp{n}=sigF(:,matrix(:,k)==n,k,nth);        %tempR  seg,some rep,1,1;
                    if ~isempty(temp{n})
                       plot(1:seg,temp{n},'Color',hsv2rgb(n/status,.5,1),'linewidth',.25);
                        %plot(1:seg,temp{n},'Color',hsv2rgb(0.8,0.5,1),'linewidth',.5); %changed the colours for clarity,
                        plot(1:seg,nanmean(temp{n},2),'Color',hsv2rgb(n/status,0.2,1),'linewidth',2);
                    end
                try    
                    rho = corr(temp{n});
                    cor(n)=mean(rho(~triu(ones(size(rho,1)))));
                catch
                    cor(n)=NaN;
                end
                end
                plot([win_sig(1) win_sig(1)] ,[ymin ymax] ,'b-','linewidth',1);% stimulus on time
                plot([win_sig(2) win_sig(2)], [ymin ymax],'b--','linewidth',1);% stimulus off time
                axis([1 seg 0 ymax]);
                axis off;
                title(sprintf('S%.2fR%.2f',cor),'Fontsize',6)
            
        end
        text(seg*1.05,0,['N' num2str(nth)],'HorizontalAlignment','left');
        plot([seg seg], [0 ymax/2],'k-','linewidth',1);% stimulus off time
        text(seg*1.05,ymax/2,[num2str(round(ymax/2))],'HorizontalAlignment','left');

    end
    drawnow;
%     saveas(h(j),fullfile(folder,h(j).Name),'png');
%     close(h(j))
end

