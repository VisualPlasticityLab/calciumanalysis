function [peakR,xR,yR]= sigRcmp(sigF,win_sig,matrix)
%sigF  seg,rep,Var,ncell;
%sigR,sigS seg*Var*ncell;
%maxtrix rep*Var
Var=size(sigF,3);
ncell = size(sigF,4);
xR=[];
yR=[];

for k=1:Var
    tempR=sigF(:,matrix(:,k),k,:);        %tempR  seg,some rep,1,ncell
    rep = size(tempR,2);
    [peakR(1,k,:),tempPeak,~]=cal_ER(tempR,win_sig); %idxR(1,k,:)
    xR = [xR , ones(1,rep)*k];
    yR = cat(1,yR,reshape(tempPeak,rep,ncell));
end