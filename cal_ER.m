function [baseA,baseSD,peakA,STD,cal_ER_option]=cal_ER(sig,win_sig,cal_option)
%calculate responses of each stimulus after averaging
%sigF:seg,rep,Var,ncell
%peak 1*Var*ncell
%error:1,Var*ncell
% cal_ER_option = 'ca_w_baseadj';
% if nvargin<3
%     cal_option ='';
% end
% cal_ER_option = [cal_option '_wo_baseadj'];
cal_ER_option = ['sp_wo_baseadj'];
debug=0;
%% clean data
% sigF = cleandata(sigF,win_sig);
%%
seg=size(sig,1);
rep=size(sig,2);
Var=size(sig,3); % Var=1
ncell=size(sig,4);
sigF = sig;
%% Each response and its deviation
time_3s = win_sig(end)-win_sig(1)+1;

%all 3s after stim onset
peakFull = nanmean(sigF(win_sig(1):win_sig(end),:,:,:),1);

%1s before stim onset
baseline= nanmean(sigF(1:win_sig(1)-1,:,:,:),1);
baseSTD = nanstd(sigF(1:win_sig(1)-1,:,:,:));

%first 3s after stim onset - BR 120320
peakF = nanmean(sigF(win_sig(1):win_sig(1)+time_3s,:,:,:),1);

%standard deviation among different reptition
STD=nanstd(peakF,0,2); 
STD=reshape(STD,1,Var,ncell);

%standard error among different reptition
STE=nanstd(peakF,0,2)/sqrt(rep); 

% mean response of all reptitions + rectify
baseA = nanmean(baseline,2);
baseSD= nanmean(baseSTD,2);
peakA = nanmean(peakF,2);
peakAll = nanmean(peakFull,2);

STE=reshape(STE,1,Var,ncell);
baseSD= reshape(baseSD,1,Var,ncell);
baseA = reshape(baseA,1,Var,ncell);
peakA = reshape(peakA,1,Var,ncell); %%just changed to match the OPTION.A dimension. 10/7
peakAll = reshape(peakAll,1,Var,ncell);

%% Plot
if debug
    figure;hold on;
    i=1;
    plot(squeeze(sigF(:,:,i,14)),'Linewidth',2);
    plot(squeeze(sigA(:,:,i,14)),'k-','Linewidth',3);
end
