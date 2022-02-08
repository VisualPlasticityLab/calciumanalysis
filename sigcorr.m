function [corr_v,corr_p] = sigcorr(sig0,vel,fast)

%corr_v: correlation value
%corr_p : correlation probability ->significance
if nargin<3
    fast = 0;
end

if size(sig0,1)~=size(vel)
    vel = vel';
end

ncell = size(sig0,2);

for nn=1:ncell
     corr_v(nn) = corr(sig0(:,nn),vel);
end

if ~fast 
     for nn=1:ncell
        for rep = 1:1000
                 tempv(rep)= corr(shuffle(sig0(:,nn)),vel);
             end
             threshold = prctile(tempv,[1 99]);
             corr_p(nn) = corr_v(nn)<threshold(1) || corr_v(nn)>threshold(2);
     end
else
    corr_p = zeros(size(corr_v));
end

figure;
if ~fast
    histplt2(corr_v(corr_p),corr_v(~corr_p))
else
    hist(corr_v,-max(abs(corr_v)):.05:max(abs(corr_v)))
end
 