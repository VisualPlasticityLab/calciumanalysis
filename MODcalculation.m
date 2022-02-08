function MOD=MODcalculation(ori,fit)
% Given fitted values to caluclate orientation to calculate ODI for A and B


nSteps = numel(peak);
nSteps1 = nSteps-mod(nSteps,2);
%oppo_ori = mod(pref_ori+nSteps1/2-1,nSteps1)+1;

A = peakR(pref_ori);
A180=peakR(pref_ori+nSteps1/2);

B = peakS(pref_ori);
B180=peakS(pref_ori+nSteps1/2);


MOD.pref_ori = (A+A180-B-B180)/(A+A180+B+B180);