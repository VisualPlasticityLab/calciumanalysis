function [sigP,variant,best,hruncmp]=runcmp(sigT,window,matrix)
global info;
seg=size(sigT,1);
rep=size(sigT,2);
ncell=size(sigT,4);
% sigF:seg,rep,Var,ncell

%%  calculate best direction and best Contrast
[~,peakR,errorR,~,peakS,errorS]= sigFcmp(sigT,window,matrix);  % calculate peak from prestim:prestim+stimON-1
% use all the peaks to calculate bestOri and bestContrast as reference
[peak,~,error]=cal_ER(sigT,window);%peak: 1*Var*ncell


peak_s=reshape(peakS,info.steps(2),info.steps(1),ncell);% peak:'contrast','Orientation',ncell
peak_r=reshape(peakR,info.steps(2),info.steps(1),ncell);% peak:'contrast','Orientation',ncell
peak=reshape(peak,info.steps(2),info.steps(1),ncell);


[~,bestDir_s]=nanmax(nanmean(peak_s,1),[],2);%find best ori for all different contrast for each cell
bestDir_s=reshape(bestDir_s,1,ncell);
[~,bestCtr_s]=nanmax(nanmean(peak_s,2),[],1);%find best contrast for all different ori for each cell
bestCtr_s=reshape(bestCtr_s,1,ncell);

[~,bestDir_r]=nanmax(nanmean(peak_r,1),[],2);%find best ori for all different contrast for each cell
bestDir_r=reshape(bestDir_r,1,ncell);
[~,bestCtr_r]=nanmax(nanmean(peak_r,2),[],1);%find best contrast for all different ori for each cell
bestCtr_r=reshape(bestCtr_r,1,ncell);


[~,bestDir]=nanmax(nanmean(peak,1),[],2);%find best ori for all different contrast for each cell
bestDir=reshape(bestDir,1,ncell);
[~,bestCtr]=nanmax(nanmean(peak,2),[],1);%find best contrast for all different ori for each cell
bestCtr=reshape(bestCtr,1,ncell);


hruncmp(1)=figure('Name','best Direction and best Contrast','position',[200 400 500 500]);
sub1=subplot(2,2,1);roseplt2(sub1,bestDir_r-bestDir,bestDir_s-bestDir);title('best Direction');
sub2=subplot(2,2,2);histplt2(bestCtr_r-bestCtr,bestCtr_s-bestCtr);title('best Contrast');

subplot(2,2,3);hold on;
scatter(bestDir_s-bestDir,bestDir_r-bestDir,'jitter','on', 'jitterAmount', 0.15);
plot([-info.steps(1) info.steps(1)],[-info.steps(1) info.steps(1)],'--')
xlabel('best Direction shift when still');ylabel('best Direction shift when running');

subplot(2,2,4); hold on;
scatter(bestCtr_s-bestCtr,bestCtr_r-bestCtr,'jitter','on', 'jitterAmount', 0.15);
plot([-info.steps(2) info.steps(2)],[-info.steps(2) info.steps(2)],'--')
xlabel('best Contrast still');ylabel('best Contrast running');

%% bin the Direction to Orientation, calculate best Orietation and best Contrast
peak_s=(peak_s(:,1:info.steps(1)/2,:)+peak_s(:,info.steps(1)/2+1:end,:))/2;
peak_r=(peak_r(:,1:info.steps(1)/2,:)+peak_r(:,info.steps(1)/2+1:end,:))/2;
peak=(peak(:,1:info.steps(1)/2,:)+peak(:,info.steps(1)/2+1:end,:))/2;

[~,bestOri_s]=nanmax(nanmean(peak_s,1),[],2);%find best ori for all different contrast for each cell
bestOri_s=reshape(bestOri_s,1,ncell);
[~,bestCtr_s]=nanmax(nanmean(peak_s,2),[],1);%find best contrast for all different ori for each cell
bestCtr_s=reshape(bestCtr_s,1,ncell);

[~,bestOri_r]=nanmax(nanmean(peak_r,1),[],2);%find best ori for all different contrast for each cell
bestOri_r=reshape(bestOri_r,1,ncell);
[~,bestCtr_r]=nanmax(nanmean(peak_r,2),[],1);%find best contrast for all different ori for each cell
bestCtr_r=reshape(bestCtr_r,1,ncell);

[~,bestOri]=nanmax(nanmean(peak_r,1),[],2);%find best ori for all different contrast for each cell
bestOri=reshape(bestOri,1,ncell);
[~,bestCtr]=nanmax(nanmean(peak,2),[],1);%find best contrast for all different ori for each cell
bestCtr=reshape(bestCtr,1,ncell);


hruncmp(2)=figure('Name','best Ori and best Contrast','position',[700 400 500 500]);
sub1=subplot(2,2,1);roseplt2(sub1,bestOri_r-bestOri,bestOri_s-bestOri);title('best Orientation');
sub2=subplot(2,2,2);histplt2(bestCtr_r-bestCtr,bestCtr_s-bestCtr);title('best Contrast');

subplot(2,2,3);hold on
scatter(bestOri_s-bestOri,bestOri_r-bestOri,'jitter','on', 'jitterAmount', 0.15);
plot([-info.steps(1)/2 info.steps(1)/2],[-info.steps(1)/2 info.steps(1)/2],'--')
xlabel('best Orientation still');ylabel('best Orientation running');
subplot(2,2,4);hold on
scatter(bestCtr_s-bestCtr,bestCtr_r-bestCtr,'jitter','on', 'jitterAmount', 0.15);
plot([-info.steps(2) info.steps(2)],[-info.steps(2) info.steps(2)],'--')
xlabel('best Contrast still');ylabel('best Contrast running');
%% 
if ~strcmp(info.var{2},'none')
    choice = questdlg('Which aspect to analyze?','To extract feature',info.var{1},info.var{2},'both','both');
    sig_p=reshape(sigT,seg,rep,info.steps(2),info.steps(1),ncell);% peak:'contrast','Orientation',ncell
    switch choice
        case 'Contrast'   %best orientation or all orientation?
            sigP=sig_p(:,:,:,sub2ind([info.steps(1) ncell],bestOri,1:ncell));%data:seg,rep,Contrast,orientation,ncell
            %sigP=squeeze(mean(sig_P,4));%data:seg,rep,Contrast,orientation,ncell
            variant=0;
            best=bestCtr;
        case 'Orientation' %define best contrast for each cell
            sig_p=permute(sig_p,[1 2 4 3 5]);
            sigP=sig_p(:,:,:,sub2ind([info.steps(2) ncell],bestCtr,1:ncell)); %sigF:seg,rep,Ori,ncell
            variant=1;
            best=bestOri;
        case 'both'
            sigP=sigT;
            variant=2;
            best=[bestCtr bestOri];
    end
else
    sigP=sigT;
    switch info.var{1}
        case 'Contrast'
            variant=0;
            best=bestCtr;
        case 'Orientation' %define best contrast for each cell
            variant=1;
            best=bestOri;
    end
end