function [peak,peakF,deviation]=cal_ER(sig,win_sig)
%calculate responses of each stimulus after averaging
debug=0;
%sigF:seg,rep,Var,ncell
% peak 1*Var*ncell
%error:1,Var*ncell
%idx

%% filter trace : smoothing and subtract baseline
%sigF=interp1(1:seg,sigF,1:1/mag:seg);
seg=size(sig,1);
rep=size(sig,2);
Var=size(sig,3);
ncell=size(sig,4);
sigF = sig;

[b,a] = butter(12,.2,'low');           % IIR filter design
sigF = filtfilt(b,a,sigF);
baselineF=nanmin(sigF(1:win_sig(1),:,:,:));
sigF=sigF-repmat(baselineF,[seg 1 1 1 ]);
% for i = 1:ncell
%     curve=func_preprocess_data(reshape(sigF(:,:,:,i),seg*rep*Var,1),win_sig(1)+(0:rep*Var-1)*seg,win_sig(end)+(0:rep*Var-1)*seg);
%     sigF(:,:,:,i) = reshape(curve,seg,rep,Var,1);
% end

% dsigF = reshape(detrend(reshape(sigF,seg,rep*Var*ncell)),seg,rep,Var,ncell);
% dthreshold=nanmean(dsigF(1:win_sig(1),:,:,:))+2*nanstd(dsigF(win_sig(1)-3:win_sig(1),:,:,:));% threshold size: 1,1,Var,ncell
% dsigF(repmat(nanmax(dsigF(win_sig(1):win_sig(end),:,:,:))<dthreshold,[seg 1 1 1 ]))=0; % set the  trace<threshold as NaN;seg,1,Var,ncell
% dsigF(dsigF<0)=0;
% dpeakF=nanmean(dsigF(win_sig(1):win_sig(end),:,:,:),1);

%define response: prestim mean+2*std
threshold=nanmean(sigF(1:win_sig(1),:,:,:))+2*nanstd(sigF(win_sig(1)-3:win_sig(1),:,:,:));% threshold size: 1,1,Var,ncell
sigF(repmat(nanmax(sigF(win_sig(1):win_sig(end),:,:,:))<threshold,[seg 1 1 1 ]))=0; % set the  trace<threshold as NaN;seg,1,Var,ncell
sigF(sigF<0)=0;
peakF=nanmean(sigF(win_sig(1):win_sig(end),:,:,:),1);
%peak=reshape(peakF,1,Var,ncell); %%just changed to match the OPTION.A dimension. 10/7

%% exclude non-significant response 

%averaging over #of repetitions
sigA= nanmean(sig,2);  %sigA:seg,1,Var,ncell
[b,a] = butter(12,.2,'low');           % IIR filter design
sigA = filtfilt(b,a,sigA);

% baselineA=min(nanmean(sigA(1:win_sig(1),:,:,:)), nanmin(sigA));
baselineA=nanmin(sigA(1:win_sig(1),:,:,:));
sigA=sigA-repmat(baselineA,[seg 1 1 1 ]);
%define response: prestim mean+2*std
% threshold=nanmean(sigA(1:win_sig(1),:,:,:))+2*nanstd(sigA(win_sig(1)-3:win_sig(1),:,:,:));% threshold size: 1,1,Var,ncell
threshold=nanmin(sigA(1:win_sig(1),:,:,:))+3*nanstd(sigA(win_sig(1)-3:win_sig(1),:,:,:));% threshold size: 1,1,Var,ncell
sigA(repmat(nanmax(sigA(win_sig(1):win_sig(end),:,:,:))<threshold,[seg 1 1 1 ]))=0; % set the  trace<threshold as NaN;seg,1,Var,ncell
sigA(sigA<0)=0;

if debug
    testcell=debug;
temp=squeeze(sig(:,:,:,testcell));
ymax = max(reshape(temp,[],1));
tempF=squeeze(sigF(:,:,:,testcell));

figure;hold on;
plot([win_sig(1) win_sig(1)] ,[0 ymax] ,'b-','linewidth',1);% stimulus on time
plot([win_sig(end) win_sig(end)], [0 ymax],'b--','linewidth',1);% stimulus off time
% plot(temp,'--');
plot(tempF,'-');
plot(nanmean(temp,2),'r--','Linewidth',1);
plot(nanmean(tempF,2),'r--','Linewidth',2);
plot(squeeze(sigA(:,:,:,testcell)),'r-','Linewidth',2);

end
%% OPTION A. maxium
% [peak,idx]=max(sigA(window,:,:,:),[],1);  % find peak from local maxium only 1*1*Var*ncell
% idx=idx+window(1)-1;
% peak=reshape(peak,1,Var,ncell);% peak 1*Var*ncell
% idx=reshape(idx,1,Var,ncell);  %idx 1*Var*ncell

%% OPTION B. mean of the response window
peak=nanmean(sigA(win_sig(1):win_sig(end),:,:,:),1);
peak=reshape(peak,1,Var,ncell); %%just changed to match the OPTION.A dimension. 10/7
idx=[];
%figure;hold on;plot(reshape(sigA,[],Var*ncell));plot(idx(:),peak(:),'*');
%% error calculation need to be updated, because of the zeros

%error_all=std(sigF(repmat(~badmatrix,[seg 1 1 1]),2))  %error:seg,1,Var,ncell
%temp=reshape(error_all,seg,Var,ncell);%temp:seg,Var,ncell
%error=temp(:,idx);
% if rep =1 , or using mean response would cause error
try
    temp=permute(sigF,[2 1 3 4 ]);  %temp:rep,seg,Var,ncell
    pos=sub2ind([seg,Var,ncell],idx,reshape(1:prod(size(idx)),1,Var,ncell));
    deviation=nanstd(temp(:,pos)); % deviation: 1,1,Var,ncell
    deviation=reshape(deviation,1,Var,ncell);
catch
%     deviation=zeros(size(peak)); % rep=1
    deviation=nanstd(sigA,0,1); 
    deviation=reshape(deviation,1,Var,ncell);

end
 