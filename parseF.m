function [sigP,variant]=parseF(sigF,window,matrix)
global info;
seg=size(sigF,1);
rep=size(sigF,2);
ncell=size(sigF,4);
% sigF:seg,rep,Var,ncell


% %%  calculate best direction and best Contrast
% [~,peakR,errorR,~,peakS,errorS]= sigFcmp(sigF,window,matrix);  % calculate peak from prestim:prestim+stimON-1
% % use peakS to calculate bestOri and bestContrast
% % peak=max(mean(sigF(,2),[],1);  % peak 1*1*Var*ncell
% 
% peak_s=reshape(peakS,info.steps(2),info.steps(1),ncell);% peak:'contrast','Orientation',ncell
% peak_r=reshape(peakR,info.steps(2),info.steps(1),ncell);% peak:'contrast','Orientation',ncell
% 
% [~,bestDir_s]=nanmax(nanmean(peak_s,1),[],2);%find best ori for all different contrast for each cell
% bestDir_s=reshape(bestDir_s,1,ncell);
% [~,bestCtr_s]=nanmax(nanmean(peak_s,2),[],1);%find best contrast for all different ori for each cell
% bestCtr_s=reshape(bestCtr_s,1,ncell);
% 
% [~,bestDir_r]=nanmax(nanmean(peak_r,1),[],2);%find best ori for all different contrast for each cell
% bestDir_r=reshape(bestDir_r,1,ncell);
% [~,bestCtr_r]=nanmax(nanmean(peak_r,2),[],1);%find best contrast for all different ori for each cell
% bestCtr_r=reshape(bestCtr_r,1,ncell);
% 
% hparse(1)=figure('Name','best Direction and best Contrast','position',[200 400 500 500]);
% sub1=subplot(2,2,1);roseplt2(sub1,bestDir_r,bestDir_s);title('best Direction');
% sub2=subplot(2,2,2);histplt2(bestCtr_r,bestCtr_s);title('best Contrast');
% 
% subplot(2,2,3);hold on;
% scatter(bestDir_s,bestDir_r,'jitter','on', 'jitterAmount', 0.15);
% plot([0 max(bestDir_s)],[0 max(bestDir_s)],'--')
% xlabel('best Direction still');ylabel('best Direction running');
% 
% subplot(2,2,4); hold on;
% scatter(bestCtr_s,bestCtr_r,'jitter','on', 'jitterAmount', 0.15);
% plot([0 max(bestCtr_s)],[0 max(bestCtr_s) ],'--')
% xlabel('best Contrast still');ylabel('best Contrast running');
% 
% %% bin the direction to Orientation, calculate best Orietation and best Contrast
% peak_s=(peak_s(:,1:info.steps(1)/2,:)+peak_s(:,info.steps(1)/2+1:end,:))/2;
% peak_r=(peak_r(:,1:info.steps(1)/2,:)+peak_r(:,info.steps(1)/2+1:end,:))/2;
% 
% [~,bestOri_s]=nanmax(nanmean(peak_s,1),[],2);%find best ori for all different contrast for each cell
% bestOri_s=reshape(bestOri_s,1,ncell);
% [~,bestCtr_s]=nanmax(nanmean(peak_s,2),[],1);%find best contrast for all different ori for each cell
% bestCtr_s=reshape(bestCtr_s,1,ncell);
% 
% [~,bestOri_r]=nanmax(nanmean(peak_r,1),[],2);%find best ori for all different contrast for each cell
% bestOri_r=reshape(bestOri_r,1,ncell);
% [~,bestCtr_r]=nanmax(nanmean(peak_r,2),[],1);%find best contrast for all different ori for each cell
% bestCtr_r=reshape(bestCtr_r,1,ncell);
% 
% 
% hparse(2)=figure('Name','best Ori and best Contrast','position',[700 400 500 500]);
% sub1=subplot(2,2,1);roseplt2(sub1,bestOri_r,bestOri_s);title('best Orientation');
% sub2=subplot(2,2,2);histplt2(bestCtr_r,bestCtr_s);title('best Contrast');
% 
% subplot(2,2,3);hold on
% scatter(bestOri_s,bestOri_r,'jitter','on', 'jitterAmount', 0.15);
% plot([0 max(bestOri_s)],[0 max(bestOri_s)],'--');
% xlabel('best Orientation still');ylabel('best Orientation running');
% subplot(2,2,4);hold on
% scatter(bestCtr_s,bestCtr_r,'jitter','on', 'jitterAmount', 0.15);
% plot([0 max(bestCtr_s)],[0 max(bestCtr_s)],'--');
% xlabel('best Contrast still');ylabel('best Contrast running');

%% averaged Ori and averaged contrast
% hparse(3)=figure('Name','averaged Ori and averaged contrast','position',[ 200 100 800 200]);
% subplot(1,2,1);imagesc(squeeze(nanmean(peak_r,1)));colorbar;
% subplot(1,2,2);imagesc(squeeze(nanmean(peak_r,2)));colorbar;

% choice = questdlg('Which aspect to analyze?','To extract feature',info.var{1},info.var{2},'both','both');
choice ='both';
if size(sigF,3)> prod(info.steps)
    sigF=sigF(:,:,1:prod(info.steps),:);
    display(sprintf('%d blank conditions removed.',size(sigF,3)-prod(info.steps)))
end
sig_p=reshape(sigF,seg,rep,info.steps(2),info.steps(1),ncell);% peak:'contrast','Orientation',ncell
switch choice
    case info.var{1}   %best orientation or all orientation?
        sigP=sig_p(:,:,:,sub2ind([info.steps(1) ncell],bestOri,1:ncell));%data:seg,rep,Contrast,orientation,ncell
        %sigP=squeeze(mean(sig_P,4));%data:seg,rep,Contrast,orientation,ncell
        variant=0;
    case info.var{2} %define best contrast for each cell
        bestCtr=1;
        sig_p=permute(sig_p,[1 2 4 3 5]);
        sigP=sig_p(:,:,:,sub2ind([info.steps(2) ncell],bestCtr,1:ncell)); %sigF:seg,rep,Ori,ncell
        variant=1;
    case 'both'
        sigP=sigF;
        variant=2;        
end
