function [peakR,sigR,errorR,peakS,sigS,errorS]= sigspFcmp(sigspF,sigwin,matrix)
%sigF  seg,rep,Var,ncell;
%sigR,sigS seg*Var*ncell;
%maxtrix rep*Var
seg=size(sigspF,1);
rep=size(sigspF,2);
Var=size(sigspF,3);
Ncell= size(sigspF,4);

for k=1:Var
    tempR=sigspF(:,matrix(:,k),k,:);        %tempR  seg,some rep,1,ncell
    sigR(:,k,:)=squeeze(mean(tempR,2));
    if ~isempty(tempR)
        [peakR(1,k,:),~,errorR(1,k,:)]=cal_spER(tempR,sigwin);
    else
         peakR(1,k,1:Ncell)= NaN;
         errorR(1,k,1:Ncell)= NaN;  
    end
    
    tempS=sigspF(:,~matrix(:,k),k,:);
    sigS(:,k,:)=squeeze(mean(tempS,2));
    if ~isempty(tempS)
        [peakS(1,k,:),~,errorS(1,k,:)]=cal_spER(tempS,sigwin);
    else
         peakS(1,k,1:Ncell)= NaN;
         errorS(1,k,1:Ncell)= NaN;
    end
    
end