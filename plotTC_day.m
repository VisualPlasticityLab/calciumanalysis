function plotTC_day(fign,selected1)

load([fign 'all.mat'],'day1','day2','SNR');

sigtype = fign(end-2:end-1);
c_ =[0.5 0.5 0.5];
c_s = [0 0 1];
c_r = [1 0 0];

% Initialize parameters
nSteps = size(eye1.sig,2);
nSteps1 = nSteps-mod(nSteps,2);
[names,names_dir] = num2ori(nSteps);
framestocapture= min(size(eye1.sig,1), size(eye2.sig,1));

for jj=selected1
            pngname = ['alltraces_cell#' num2str(jj) '_' sigtype '.png'];

            if ~exist(pngname,'file')
                
%             pngname = [h1.Name '_' sigtype '_2.png'];
    nth1= jj;
    nth2= jj;
        a = PEAK0(jj);
        odi = ODI0(jj) ;
        switch fign(end-2:end-1)
            case 'ca'
                temp_1st=squeeze(eye1.sigF(:,:,:,nth1));
                temp_2nd=squeeze(eye2.sigF(:,:,:,nth2));
%                 temp_1st= temp_1st-repmat(prctile(temp_1st,5),[framestocapture 1 1]);
%                 temp_2nd= temp_2nd-repmat(prctile(temp_2nd,5),[framestocapture 1 1]);
                ymax=prctile([temp_1st(:);temp_2nd(:)],99.9);
                ymin=prctile([temp_1st(:);temp_2nd(:)],2.5);%min([temp_1st(:);temp_2nd(:)]);
                buffer=(ymax-ymin)*.05;
           case 'sp'
                temp_1st=squeeze(eye1.sigspF(:,:,:,nth1));
                temp_2nd=squeeze(eye2.sigspF(:,:,:,nth2));
                ymax=prctile([temp_1st(:);temp_2nd(:)],99.9);
                ymin=0;%min([temp_1st(:);temp_2nd(:)]);
                buffer=0;%(ymax-ymin)*.05;
            case 'ML'
                temp_1st=squeeze(eye1.sigMLF(:,:,:,nth1));
                temp_2nd=squeeze(eye2.sigMLF(:,:,:,nth2));
                ymax=prctile([temp_1st(:);temp_2nd(:)],99.9);
                ymin=0;%min([temp_1st(:);temp_2nd(:)]);
                buffer=0;%(ymax-ymin)*.05;
        end
        
        h1=figure('Position',[100 200 1200 400],'Name',['alltraces_cell#' num2str(jj)]);
        ax(1)=subplot(3,3,1);    hold on;
        errorbar(1:nSteps1,eye1.peak(:,1:end-1,nth1),eye1.error(:,1:end-1,nth1),'o-','Linewidth',2,'Color',c_);
        errorbar(1:nSteps1,eye1.peakS(:,1:end-1,nth1),eye1.errorS(:,1:end-1,nth1),'o-','Linewidth',1,'Color',c_s);
        errorbar(1:nSteps1,eye1.peakR(:,1:end-1,nth1),eye1.errorR(:,1:end-1,nth1),'o-','Linewidth',1,'Color',c_r);
        plot([1 nSteps1],[eye1.peak(:,end,nth1) eye1.peak(:,end,nth1)],'Linewidth',2,'Color',c_);
        plot([1 nSteps1],[eye1.peakS(:,end,nth1) eye1.peakS(:,end,nth1)],'Linewidth',1,'Color',c_s);
        plot([1 nSteps1],[eye1.peakR(:,end,nth1) eye1.peakR(:,end,nth1)],'Linewidth',1,'Color',c_r);
        title(sprintf('ipsi eye, ori(all,s,r)=(%d,%d,%d)',ORI1(1,jj),ORI2([2 5],jj)));
        axis tight
           
        ax(2)=subplot(3,3,2);    hold on;
        errorbar(1:nSteps1,eye2.peak(:,1:end-1,nth2),eye2.error(:,1:end-1,nth2),'o-','Linewidth',2,'Color',c_/2);
        errorbar(1:nSteps1,eye2.peakS(:,1:end-1,nth2),eye2.errorS(:,1:end-1,nth2),'o-','Linewidth',1,'Color',c_s/2);
        errorbar(1:nSteps1,eye2.peakR(:,1:end-1,nth2),eye2.errorR(:,1:end-1,nth2),'o-','Linewidth',1,'Color',c_r/2);
        plot([1 nSteps1],[eye2.peak(:,end,nth1) eye2.peak(:,end,nth1)],'Linewidth',2,'Color',c_/2);
        plot([1 nSteps1],[eye2.peakS(:,end,nth1) eye2.peakS(:,end,nth1)],'Linewidth',1,'Color',c_s/2);
        plot([1 nSteps1],[eye2.peakR(:,end,nth1) eye2.peakR(:,end,nth1)],'Linewidth',1,'Color',c_r/2);
        title(sprintf('contra eye, ori(all,s,r)=(%d,%d,%d)',ORI1(2,jj),ORI2([3 6],jj)));
        axis tight
        
        ax(3)=subplot(3,3,3);    hold on;
        plot(1:nSteps1,eye1.peak(:,1:end-1,nth1),'o-','Linewidth',2,'Color',c_);
        plot(1:nSteps1,eye2.peak(:,1:end-1,nth2),'o-','Linewidth',2,'Color',c_/2);
        plot([1 nSteps1],[eye1.peak(:,end,nth1) eye1.peak(:,end,nth1)],'Linewidth',2,'Color',c_);
        plot([1 nSteps1],[eye2.peak(:,end,nth1) eye2.peak(:,end,nth1)],'Linewidth',2,'Color',c_/2);
%         plot(1:nSteps,eye1.peakR(:,:,nth1),'Linewidth',1,'Color',c_r);
%         plot(1:nSteps,eye1.peakS(:,:,nth1),'Linewidth',1,'Color',c_s);
%         plot(1:nSteps,eye2.peakR(:,:,nth2),'Linewidth',1,'Color',c_r/2);
%         plot(1:nSteps,eye2.peakS(:,:,nth2),'Linewidth',1,'Color',c_s/2);
        xlim([1 nSteps1]);
        legend('ispi','contra');
        legend('boxoff','Location','top') ;
        try
            title(sprintf('SNR=%.1f, a=%.2f,ori(all,s,r)=(%d,%d,%d),odi=%.2f',SNR.skew(jj),a,ORI0(jj),ORI2([1 4],jj),odi));
        catch
            title(sprintf('a=%.2f,ori(all,s,r)=(%d,%d,%d),odi=%.2f',a,ORI0(jj),ORI2([1 4],jj),odi));
        end
        linkaxes(ax)
        
        for ori=1:nSteps
            traces1r = squeeze(temp_1st(:,eye1.matrix(:,ori),ori));
            traces1s = squeeze(temp_1st(:,~eye1.matrix(:,ori),ori));
            traces1r_avg=eye1.sigR(:,ori,nth1);
            traces1s_avg=eye1.sigS(:,ori,nth1);
            traces1_avg=eye1.sig(:,ori,nth1);
            
            traces2r= squeeze(temp_2nd(:,eye2.matrix(:,ori),ori));
            traces2s= squeeze(temp_2nd(:,~eye2.matrix(:,ori),ori));
            traces2r_avg=eye2.sigR(:,ori,nth2);
            traces2s_avg=eye2.sigS(:,ori,nth2);
            traces2_avg=eye2.sig(:,ori,nth2);
            
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
            
            bx(1)=subplot(3,nSteps,nSteps+ori); hold on;
            rectangle('Position',[eye1.win_sig(1) ymin-buffer eye1.win_sig(end)-eye1.win_sig(1) ymax-ymin+2*buffer],...
                'FaceColor',[0.5 .5 .5 .3],'EdgeColor',[0.5 .5 .5 .3]);% stimulus on&off time
            plot(traces1r,'Color',[c_r 0.3],'Linewidth',.5);
            plot(traces1s,'Color',[c_s 0.3],'Linewidth',.5);
            plot(traces1_avg,'Color',c_,'Linewidth',2);
%             plot(traces1r_avg,'Color',[0 1 0 0.7],'Linewidth',1);
%             plot(traces1s_avg,'Color',[0 0 0 0.7],'Linewidth',1);
            title(sprintf('%.2f',eye1.peak(:,ori,nth1)));   %  title(sprintf('%.2f',Fr1(ori))
            
            bx(2)=subplot(3,nSteps,2*nSteps+ori); hold on;
            rectangle('Position',[eye1.win_sig(1) ymin-buffer eye1.win_sig(end)-eye1.win_sig(1) ymax-ymin+2*buffer],...
                'FaceColor',[0.5 .5 .5 .3],'EdgeColor',[0.5 .5 .5 .3]);% stimulus on&off time
            plot(traces2r,'Color',[c_r/2 0.3],'Linewidth',.5);
            plot(traces2s,'Color',[c_s/2 0.3],'Linewidth',.5);           
            plot(traces2_avg,'Color',c_/2,'Linewidth',2);
%             plot(traces2s_avg,'Linewidth',1,'Color',[0 0 0 0.7]);
%             plot(traces2r_avg,'Linewidth',2,'Color',[1 0 0 0.7]);
%             title(sprintf('%.2f',eye2.peak(:,ori,nth2)));   % title(sprintf('%.2f',Fr2(ori)))
            bx(2).XLabel.String=names_dir{ori};

            bx(1).XTickLabel='';
            bx(2).XTickLabel='';
            if ori==1
                bx(1).YLabel.String=sprintf('ipsieye,cell# %d',nth1);
                bx(2).YLabel.String=sprintf('contraeye,cell# %d',nth2);
            else
                bx(1).YTickLabel='';
                bx(2).YTickLabel='';
            end
            linkaxes(bx);
            axis([1 framestocapture ymin-buffer ymax+buffer]);

        end
        saveas(h1,pngname)
        close(h1)
            end
end
