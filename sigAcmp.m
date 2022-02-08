function [sigR,peakR,errorR,sigS,peakS,errorS]= sigAcmp(sigF,win_sig,matrix)
%sigF  seg,rep,Var,ncell;
%sigR,sigS seg*Var*ncell;
%maxtrix rep*Var
seg = size(sigF,1);
Var=size(sigF,3);
ncell = size(sigF,4);
sigR = nan(seg,Var,ncell);
sigS = nan(seg,Var,ncell);
peakR = nan(1,Var,ncell);
peakS = nan(1,Var,ncell);

for k=1:Var
    %%    
    tempR=sigF(:,matrix(:,k),k,:);        %tempR  seg,some rep,1,ncell
    sigR(:,k,:)=squeeze(mean(tempR,2));
    [peakR(1,k,:),~,errorR(1,k,:)]=cal_ER(tempR,win_sig); %idxR(1,k,:)  
    %%
    tempS=sigF(:,~matrix(:,k),k,:);
%     if ~isempty(tempS)
%         sigS(:,k,:)=squeeze(mean(tempS,2));
%         [peakS(1,k,:),~,errorS(1,k,:)]=cal_ER(tempS,win_sig); %idxS(1,k,:)
%     else
        peakS(1,k,:)=NaN;
        errorS(1,k,:)=NaN;
        sigS(:,k,:)=NaN;
%     end
end