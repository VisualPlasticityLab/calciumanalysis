function [peakA,sigA,STE]=cal_spER(sig,win_sig)
%calculate responses of each stimulus after averaging
%sigF:seg,rep,Var,ncell
%peak 1*Var*ncell
%error:1,Var*ncell

debug=0;
seg=size(sig,1);
rep=size(sig,2);
Var=size(sig,3); % Var=1
ncell=size(sig,4);
sigF = sig;
%% For spiking data, only use 1s data 
win_sig(end) =  ceil((win_sig(end)-win_sig(1)+1)/3)+win_sig(1)-1;
%% Normalize with baseline
% baseline= mean(sigF(1:win_sig(1)-1,:,:,:));
% sigF=sigF-repmat(baseline,[seg 1 1 1]);

%% Calculate standard error based on each trial
peakF = nanmean(sigF(win_sig(1):win_sig(end),:,:,:),1);
STE=nanstd(peakF,0,2)/sqrt(rep); 
STE=reshape(STE,1,Var,ncell);
%% Total response OR maxium response    
sigA = squeeze(nanmean(sigF,2));
peakA = nansum(sigA(win_sig(1):win_sig(end),:,:,:),1); %total three sec response

% %%
% sigA(sigA<0)=0;
% peakA(peakA<0)=0;

% onesecbin = round((win_sig(end)-win_sig(end))/3);
% peakA1=nanmean(sigA(win_sig(1):win_sig(1)+onesecbin,:,:,:),1);
% peakA2=nanmean(sigA(win_sig(1)+onesecbin+1:win_sig(1)+onesecbin*2,:,:,:),1);
% peakA3=nanmean(sigA(win_sig(1)+onesecbin*2+1:win_sig(end),:,:,:),1);
% peakA = nanmax(cat(1,peakA1,peakA2,peakA3));

peakA=reshape(peakA,1,Var,ncell); %%just changed to match the OPTION.A dimension. 10/7
%% Plot
if debug
    figure;hold on;
    i=1;
    plot(squeeze(sigF(:,:,i,14)),'Linewidth',2);
    plot(squeeze(sigA(:,:,i,14)),'k-','Linewidth',3);
end
