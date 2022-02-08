function [gdid,X]= AutoCellSelection(folder)
%%
if nargin<1
    [~,folder]=uigetfile('peakSI.mat');
end
load(fullfile(folder,'peakSI.mat'))
if isempty(run.matrix)
    run.matrix = logical(zeros(size(sigF,2),size(sigF,3)));
end
%% calculate correlation and pvalue
for nth = 1:size(sigF,4) %the last is noise
    for kth=1:size(sigF,3)
        if sum(~run.matrix(:,kth))>=2
            tempS=sigF(:,~run.matrix(:,kth),kth,nth);
            [rhoS,pS] = corr(tempS);
            corS(kth) = mean(rhoS(~triu(ones(size(rhoS,1)))));
            pvalS(kth) = mean(pS(~triu(ones(size(rhoS,1)))));
            peakS(kth) = max(sum(tempS));
        else
            corS(kth) = NaN;
            pvalS(kth) = NaN;
            peakS(kth) = NaN;
        end
        if sum(run.matrix(:,kth))>=2
            tempR=sigF(:,run.matrix(:,kth),kth,nth);
            [rhoR,pR]  = corr(tempR);
            corR(kth) = mean(rhoR(~triu(ones(size(rhoR,1)))));
            pvalR(kth) = mean(pR(~triu(ones(size(rhoR,1)))));
            peakR(kth) = max(sum(tempR));
        else
            corR(kth)= NaN;
            pvalR(kth)=NaN;
            peakR(kth) = NaN;
        end
    end
    [corV(nth),idx] = nanmax(cat(2,corS,corR));
    temp = cat(2,pvalS,pvalR);
    temp2 = cat(2,peakS,peakR);
    pV(nth) = temp(idx);
    peak(nth) = temp2(idx);
end

%% k-means
X=[corV;pV;peak]';
opts = statset('Display','final');
[cidx, ctrs] = kmeans(X(:,1:2), 2, 'Distance','city', ...
    'Replicates',10, 'Options',opts);
ctrs(:,3) =[ mean(X(cidx==1,3));mean(X(cidx==2,3))];
figure;
plot3(X(cidx==1,1),X(cidx==1,2),X(cidx==1,3),'r.', ...
     X(cidx==2,1),X(cidx==2,2),X(cidx==2,3),'b.', ctrs(:,1),ctrs(:,2),ctrs(:,3),'kx','MarkerSize',10);%ctrs(:,3)
xlabel('Correlation')
ylabel('P-value')
 zlabel('Amplitude')
legend(['group1=red,N=' num2str(sum(cidx==1))] ,['group2=blue,N=' num2str(sum(cidx==2))]);
% if ctrs(1,3)>ctrs(2,3)
if ctrs(1,2)>ctrs(2,2)
    gdid = cidx==1;
else
    gdid = cidx==2;
end
