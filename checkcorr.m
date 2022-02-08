function [corr_v,corr_s] = sigcorr(sig0,vel)

%corr_v: correlation value
%corr_s : correlation significance

nn = 63;
 ncell = size(sig0,2);

for nn=1:ncell
 
 corr_v(nn) = corr(sig0(:,nn),vel);
 clear tempv
 for rep = 1:1000
     tempv(rep)= corr(shuffle(sig0(:,nn)),vel);
 end
 
 threshold = prctile(tempv,[5 95]);
 
 if corr_v(nn)<threshold(1) || corr_v(nn)>threshold(1)
     corr_s(nn) =1;
 else
     corr_s(nn) = 0;
 end
end