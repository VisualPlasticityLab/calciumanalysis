function plotTC(fign,selected1)
load([fign 'all.mat'],'PEAK0','ORI0','ORI1','ORI2','ODI0','ODI1','eye1','eye2','SNR');
% eye1.win_sig = [ 7 30];
c_ =[0.5 0.5 0.5];
c_s = [0 0 1];
c_r = [1 0 0];

%% Initialize parameters

    
sigtype = fign(end-2:end-1);
nSteps = size(eye1.peak,2);
nSteps1 = nSteps-mod(nSteps,2);
   if ~isfield(eye1,'sig')
       eye1.sig=eye1.sigF;
       eye2.sig=eye2.sigF;
   end
framestocapture= min(size(eye1.sig,1), size(eye2.sig,1));
[names,names_dir] = num2ori(nSteps);
%%
for jj=selected1
            pngname = ['alltraces_cell#' num2str(jj) '_' sigtype '.png'];
            figname = ['alltraces_cell#' num2str(jj) '_' sigtype '.fig'];

            if ~exist(pngname,'file')
                
%             pngname = [h1.Name '_' sigtype '_2.png'];
    nth1= jj;
    nth2= jj;
        a = PEAK0(jj);
        odi = ODI1(:,jj) ;
        switch sigtype
            case 'ca'
                temp_1st=squeeze(eye1.sigF(:,:,:,nth1));
                temp_2nd=squeeze(eye2.sigF(:,:,:,nth2));
%                 temp_1st= temp_1st-repmat(prctile(temp_1st,5),[framestocapture 1 1]);
%                 temp_2nd= temp_2nd-repmat(prctile(temp_2nd,5),[framestocapture 1 1]);
                ymax=prctile([temp_1st(:);temp_2nd(:)],99);
%                 ymax = max(1,ymax);
%                 ymin=prctile([temp_1st(:);temp_2nd(:)],1);%min([temp_1st(:);temp_2nd(:)]);
%                 ymax = 4;
                 ymin = 0;
                buffer=(ymax-ymin)*.05;
           case 'sp'
                temp_1st=squeeze(eye1.sigF(:,:,:,nth1));%eye1.sigspF, in old fiels
                temp_2nd=squeeze(eye2.sigF(:,:,:,nth2));
                ymax=prctile([temp_1st(:);temp_2nd(:)],99);
                ymin=0;%min([temp_1st(:);temp_2nd(:)]);
                buffer=(ymax-ymin)*.025;
            case 'ML'
                temp_1st=squeeze(eye1.sigMLF(:,:,:,nth1));
                temp_2nd=squeeze(eye2.sigMLF(:,:,:,nth2));
                ymax=prctile([temp_1st(:);temp_2nd(:)],99.9);
                ymin=0;%min([temp_1st(:);temp_2nd(:)]);
                buffer=(ymax-ymin)*.025;
        end
        
        h1=figure('Position',[100 200 1200 400],'Name',['alltraces_cell#' num2str(jj)]);
        ax(1)=subplot(2,2,1);    hold on;
        if nSteps1<nSteps
        errorbar(1:nSteps1,eye1.peak(:,1:end-1,nth1),eye1.error(:,1:end-1,nth1),'o-','Linewidth',2,'Color',c_);
        errorbar(1:nSteps1,eye1.peakS(:,1:end-1,nth1),eye1.errorS(:,1:end-1,nth1),'o-','Linewidth',1,'Color',c_s);
        errorbar(1:nSteps1,eye1.peakR(:,1:end-1,nth1),eye1.errorR(:,1:end-1,nth1),'o-','Linewidth',1,'Color',c_r);
        plot([1 nSteps1],[eye1.peak(:,end,nth1) eye1.peak(:,end,nth1)],'Linewidth',2,'Color',c_);
        plot([1 nSteps1],[eye1.peakS(:,end,nth1) eye1.peakS(:,end,nth1)],'Linewidth',1,'Color',c_s);
        plot([1 nSteps1],[eye1.peakR(:,end,nth1) eye1.peakR(:,end,nth1)],'Linewidth',1,'Color',c_r);
        else
        errorbar(1:nSteps,eye1.peak(:,:,nth1),eye1.error(:,:,nth1),'o-','Linewidth',2,'Color',c_);
        errorbar(1:nSteps,eye1.peakS(:,:,nth1),eye1.errorS(:,:,nth1),'o-','Linewidth',1,'Color',c_s);
        errorbar(1:nSteps,eye1.peakR(:,:,nth1),eye1.errorR(:,:,nth1),'o-','Linewidth',1,'Color',c_r);
        end  
        %         title(sprintf('ipsi eye, ori(all,s,r)=(%d,%d,%d)',ORI1(1,jj),ORI2([2 5],jj)));
        title('Ipsi')
        legend({'all','still','run'},'orientation','horizontal');
        legend('boxoff','Location','top') ;
        ax(1).XAxis.TickValues = 1:2:nSteps;
        ax(1).XAxis.TickLabels = names_dir(1:2:end);
          
        ax(2)=subplot(2,2,2);    hold on;
        if nSteps1<nSteps

        errorbar(1:nSteps1,eye2.peak(:,1:end-1,nth2),eye2.error(:,1:end-1,nth2),'o-','Linewidth',2,'Color',c_/2);
        errorbar(1:nSteps1,eye2.peakS(:,1:end-1,nth2),eye2.errorS(:,1:end-1,nth2),'o-','Linewidth',1,'Color',c_s/2);
        errorbar(1:nSteps1,eye2.peakR(:,1:end-1,nth2),eye2.errorR(:,1:end-1,nth2),'o-','Linewidth',1,'Color',c_r/2);
        plot([1 nSteps1],[eye2.peak(:,end,nth1) eye2.peak(:,end,nth1)],'Linewidth',2,'Color',c_/2);
        plot([1 nSteps1],[eye2.peakS(:,end,nth1) eye2.peakS(:,end,nth1)],'Linewidth',1,'Color',c_s/2);
        plot([1 nSteps1],[eye2.peakR(:,end,nth1) eye2.peakR(:,end,nth1)],'Linewidth',1,'Color',c_r/2);
        else
            errorbar(1:nSteps,eye2.peak(:,:,nth2),eye2.error(:,:,nth2),'o-','Linewidth',2,'Color',c_/2);
        errorbar(1:nSteps,eye2.peakS(:,:,nth2),eye2.errorS(:,:,nth2),'o-','Linewidth',1,'Color',c_s/2);
        errorbar(1:nSteps,eye2.peakR(:,:,nth2),eye2.errorR(:,:,nth2),'o-','Linewidth',1,'Color',c_r/2);
        end
                    %         title(sprintf('contra eye, ori(all,s,r)=(%d,%d,%d)',ORI1(2,jj),ORI2([3 6],jj)));
        title('Contra')
        legend({'all','still','run'},'orientation','horizontal');
        legend('boxoff','Location','top') ;            
        axis tight
        ax(2).XAxis.TickValues = 1:2:nSteps;
        ax(2).XAxis.TickLabels = names_dir(1:2:end);
       
%         ax(3)=subplot(2,3,3);    hold on;
%         plot(1:nSteps1,eye1.peak(:,1:end-1,nth1),'o-','Linewidth',2,'Color',c_);
%         plot(1:nSteps1,eye2.peak(:,1:end-1,nth2),'o-','Linewidth',2,'Color',c_/2);
%         plot([1 nSteps1],[eye1.peak(:,end,nth1) eye1.peak(:,end,nth1)],'Linewidth',2,'Color',c_);
%         plot([1 nSteps1],[eye2.peak(:,end,nth1) eye2.peak(:,end,nth1)],'Linewidth',2,'Color',c_/2);
% %         plot(1:nSteps,eye1.peakR(:,:,nth1),'Linewidth',1,'Color',c_r);
% %         plot(1:nSteps,eye1.peakS(:,:,nth1),'Linewidth',1,'Color',c_s);
% %         plot(1:nSteps,eye2.peakR(:,:,nth2),'Linewidth',1,'Color',c_r/2);
% %         plot(1:nSteps,eye2.peakS(:,:,nth2),'Linewidth',1,'Color',c_s/2);

%         ax(3).XAxis.TickLabels = names_dir(1:end-1);
%         legend('ispi','contra');
%         legend('boxoff','Location','top') ;
        try
%            title(sprintf('SNR=%.1f, a=%.2f,ori(all,s,r)=(%d,%d,%d),odi=%.1f %.1f',SNR.skew(jj),a,ORI0(jj),ORI2([1 4],jj),odi));
           title(sprintf('SNR=%.1f, a=%.2f,ori(s,r)=%d,%d,odi(s,r)=%.1f,%.1f',SNR.skew(jj),a,ORI2([1 4],jj),odi));
        catch
            title(sprintf('a=%.2f,ori(s,r)=%d,%d,odi(s,r)=%.1f,%.1f',a,ORI2([1 4],jj),odi));
        end
        xlim([1 nSteps1]);
        ylim([ymin-buffer ymax+buffer])
        axis tight
        linkaxes(ax)
        
        for ori=1:nSteps
            traces1r = squeeze(temp_1st(:,eye1.matrix(:,ori),ori));
            traces1s = squeeze(temp_1st(:,~eye1.matrix(:,ori),ori));
            traces1r_avg=nanmean(traces1r,2);
            traces1s_avg=nanmean(traces1s,2);
            traces1_avg=nanmean(temp_1st(:,:,ori),2);
            try
                rho = corr(squeeze(temp_1st(eye1.win_sig(1):eye1.win_sig(end)...
                    ,~eye1.matrix(:,ori),ori)));
                cor1(1)=mean(rho(~triu(ones(size(rho,1)))));
            end
            try
                rho = corr(squeeze(temp_1st(eye1.win_sig(1):eye1.win_sig(end)...
                    ,eye1.matrix(:,ori),ori)));
                cor1(2)=mean(rho(~triu(ones(size(rho,1)))));
            end
            traces2r= squeeze(temp_2nd(:,eye2.matrix(:,ori),ori));
            traces2s= squeeze(temp_2nd(:,~eye2.matrix(:,ori),ori));
            traces2r_avg=nanmean(traces2r,2);
            traces2s_avg=nanmean(traces2s,2);
            traces2_avg=nanmean(temp_2nd(:,:,ori),2);
            try
                rho = corr(squeeze(temp_2nd(eye1.win_sig(1):eye1.win_sig(end)...
                             ,~eye1.matrix(:,ori),ori)));
                cor2(1)=mean(rho(~triu(ones(size(rho,1)))));
            end
            try
                rho = corr(squeeze(temp_2nd(eye1.win_sig(1):eye1.win_sig(end)...
                    ,eye1.matrix(:,ori),ori)));
                cor2(2)=mean(rho(~triu(ones(size(rho,1)))));
            end
            try
                Fr1(ori)=pwrplt(squeeze(nanmean(temp_1st(:,:,ori),2)),1);
            catch
                Fr1(ori) = nan;
            end
            try
                Fr2(ori)=pwrplt(squeeze(nanmean(temp_2nd(:,:,ori),2)),1);
            catch
                Fr2(ori)= nan;
            end
            
            bx(1)=subplot(2,nSteps*2+1,nSteps*2+ori+1); hold on;
                title(sprintf('S%.1f\nR%.1f',cor1),'Fontsize',8)

            plot(traces1r,'Color',[c_r 0.3],'Linewidth',.5);
            plot(traces1s,'Color',[c_s 0.3],'Linewidth',.5);
%             plot(traces1_avg,'Color',c_,'Linewidth',2);
            plot(traces1r_avg,'Color',[c_r 0.7],'Linewidth',1);
            plot(traces1s_avg,'Color',[c_s 0.7],'Linewidth',1);
%             title(sprintf('%.2f',eye1.peak(:,ori,nth1)));   %  title(sprintf('%.2f',Fr1(ori))
            axis off
            rectangle('Position',[eye1.win_sig(1) ymin-2*buffer eye1.win_sig(end)-eye1.win_sig(1) 2*buffer],...
                'FaceColor',[.5 .5 .5 1],'EdgeColor',[.5 .5 .5 1]);% stimulus on&off time
            
            bx(2)=subplot(2,nSteps*2+1,nSteps*3+ori+2); hold on;
                title(sprintf('S%.1f\nR%.1f',cor2),'Fontsize',8)
            plot(traces2s,'Color',[c_s/2 0.3],'Linewidth',.5);           
            plot(traces2r,'Color',[c_r/2 0.3],'Linewidth',.5);
%             plot(traces2_avg,'Color',c_/2,'Linewidth',2);
            plot(traces2s_avg,'Linewidth',1,'Color',[c_s/2 0.7]);
            plot(traces2r_avg,'Linewidth',1,'Color',[c_r/2 0.7]);
%             title(sprintf('%.2f',eye2.peak(:,ori,nth2)));   % title(sprintf('%.2f',Fr2(ori)))
%             bx(2).XLabel.String=names_dir{ori};
            axis off
            rectangle('Position',[eye1.win_sig(1) ymin-2*buffer eye1.win_sig(end)-eye1.win_sig(1) 2*buffer],...
                'FaceColor',[0.5 .5 .5 1],'EdgeColor',[0.5 .5 .5 1]);% stimulus on&off time

            bx(1).XTickLabel='';
            bx(2).XTickLabel='';
%             scalebar('unit',{'',''})
            if ori==1
                bx(1).YLabel.String=sprintf('ipsi # %d',nth1);
                bx(2).YLabel.String=sprintf('contra # %d',nth2);
                scalebar('unit',{'','df/F'},'location','northwest')
            else
                bx(1).YTickLabel='';
                bx(2).YTickLabel='';
            end
            linkaxes(bx);
            axis([1 framestocapture ymin-buffer ymax+buffer]);
        end
        saveas(h1,pngname)
%         saveas(h1,figname)
        %close(h1)
            end
end
