function extractWF(fn)
sbxread(fn,1,1);
global info;
if info.volscan ==1
    nplanes=info.otparam(3);
else
    nplanes=1;
end
for k=1:nplanes
    ct=0;
    for i=k:nplanes:info.max_idx
        ct=ct+1;
        m=double(sbxread(fn,i-1,1));
        WF(:,ct,k) = sum(sum(m,3),2);
    end
end
info.WF=WF;
save(fn,'info','-append');

    