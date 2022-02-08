function [baseR, baseSdR,peakR,errorR,peakAllR,...
          baseS, baseSdS,peakS,errorS,peakAllS,cal_ER_option]= sigFcmp_onestimloop(sigF,sigwin,matrix)
%sigF  seg,rep,Var,ncell;
%sigR,sigS seg*Var*ncell;
%maxtrix rep*Var
seg=size(sigF,1);
rep=size(sigF,2);
Var=size(sigF,3);
Ncell= size(sigF,4);
% if no ball info means all still
if isempty(matrix)
    matrix = logical(zeros(rep,Var));
end

for k=1:Var
    tempR=sigF(:,matrix(:,k),k,:);        %tempR  seg,some rep,1,ncell
    sigR(:,k,:)=squeeze(mean(tempR,2));
    if ~isempty(tempR)                    %if there is more than 1 running trials
        [baseR(1,k,:),baseSdR(1,k,:),peakR(1,k,:),errorR(1,k,:),peakAllR(1,k,:),cal_ER_option]=cal_ER_onestimloop(tempR,sigwin);
    else
         baseR(1,k,1:Ncell)= NaN;
         baseSdR(1,k,1:Ncell)= NaN;
         peakR(1,k,1:Ncell)= NaN;
         errorR(1,k,1:Ncell)= NaN;  
         peakAllR(1,k,1:Ncell)= NaN;
    end
    
    tempS=sigF(:,~matrix(:,k),k,:);
    sigS(:,k,:)=squeeze(mean(tempS,2));
    if ~isempty(tempS)
        [baseS(1,k,:),baseSdS(1,k,:),peakS(1,k,:),errorS(1,k,:),peakAllS(1,k,:),cal_ER_option]=cal_ER_onestimloop(tempS,sigwin);
    else
         baseS(1,k,1:Ncell)= NaN;
         baseSdS(1,k,1:Ncell)= NaN;
         peakS(1,k,1:Ncell)= NaN;
         errorS(1,k,1:Ncell)= NaN;  
         peakAllS(1,k,1:Ncell)= NaN;
    end
    
end