function [decay_each,decay]=fitexp1(sig_chunk,num,debug)
% sig_chunk input y
% num: percentage of drop?

if nargin<3
debug=0;
end

[amp1,pos1] = findpeaks(sig_chunk);
% hold on;plot(pos1,amp1,'*')
[amp2,pos2]= findpeaks(-1*sig_chunk);
amp2= amp2*-1;
while pos2(1) <  pos1(1)
    pos2 = pos2(2:end);
    amp2 = amp2(2:end);
end
while pos1(end) >  pos2(end)
    pos1 = pos1(1:end-1);
    amp1 =  amp1(1:end-1);
end
drop = amp1-amp2;
dur = pos2 - pos1;

dur_thr=round(prctile(dur,(1-num/numel(drop))*100));
peak_thr = prctile(drop,50);
gd = find(dur>dur_thr & drop > peak_thr);

temp= sig_chunk(ones(numel(gd),1)*[0:dur_thr-1]+repmat(pos1(gd),1,dur_thr));
totaltrace =0;
for i=1:size(gd,2)
    totaltrace = totaltrace+ temp(i,:)'/max(temp(i,:));
    [f,~,o] = fit([0:dur_thr-1]',temp(i,:)','exp1');
    a(i)=f.a;
    b(i)=f.b;
end
meantrace = totaltrace/numel(gd);
meanf=fit([0:dur_thr-1]',meantrace,'exp1');
decay=meanf.b;
decay_each=b(i);
%% Do some plotting if debug
if debug 
figure;hold on;
plot(meanf([0:dur_thr-1]')*a,'--b')
legend(['mean&fit:' num2str(decay) 'meanoffit:' num2str(mean(decay_each)) ])
plot(temp','r')
end