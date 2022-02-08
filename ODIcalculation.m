function [ODI,ODI_F0]=ODIcalculation(pref_ori,peakA,peakB)
% Given preferred orientation to calculate ODI for A and B
% contraeye A, ipsi B, 

nSteps = size(peakA,2);
nSteps = nSteps-mod(nSteps,2);

if pref_ori ==0
    ODI =NaN;
    ODI_F0 = NaN;
    return;
end

peakA(isnan(peakA))=0;
peakB(isnan(peakB))=0;

oppo_ori = mod(pref_ori+nSteps/2-1,nSteps)+1;

A = peakA(pref_ori)+ peakA(oppo_ori);
%A = sqrt(peakA(pref_ori));
A_F0=min(peakA) ;

B = peakB(pref_ori)+peakB(oppo_ori);
%B = sqrt(peakB(pref_ori));
B_F0=min(peakB);

ODI = (A-B)/(A+B); %just use the pref_direction
ODI_F0= (A-A_F0-B+B_F0)/(A-A_F0-B+B_F0);  %USE F0-F1 modulation to 