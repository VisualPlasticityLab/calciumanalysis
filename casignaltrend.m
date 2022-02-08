function casignaltrend(sigfn)
if nargin==0
    sigfn = uigetfile('*.signals','processed combined cell signals');
end

pos_ = strfind(sigfn,'_');
fn = sigfn(1:pos_(end-3));
mp = matfile([fn 'memmap.mat']);

h1 = figure;
try
    imagesc(mp.V) ;
catch
    tmpf = uigetfile('*.align','find aligned data for m');
    tmp = load(tmpf,'m','-mat');
    imagesc(tmp.m);
end
colormap summer;
maskpic = roipoly;
try
    imshowpair(mp.V,maskpic*1000,'Scaling','joint');
catch
    imshowpair(tmp.m,maskpic*1000,'Scaling','joint');
end
title(strrep(fn(1:end-1),'_','-'))


%%
tic;
Yr = mp.Yr;
Y_noise = mean(Yr(maskpic,:));
toc;
load(sigfn,'-mat','sig','velocity');
if ~exist('velocity')
    cal_SNR;
    load(sigfn,'-mat','sig','velocity');
end
%%
S = mean(sig,2);
V = downsample(velocity,round(numel(velocity)/numel(S)));
V = V(1:numel(S));
clear velocity
S_norm = S/prctile(S,99)*prctile(V,99);
Y_norm = Y_noise/prctile(Y_noise,99)*prctile(V,99);

thr = [1,99];
[S_up,S_up_loc,S_up_f] = envelope2(S_norm,thr(1));
[Y_up,Y_up_loc,Y_up_f] = envelope2(Y_norm,thr(1));

[S_dn,S_dn_loc,S_dn_f] = envelope2(S_norm,thr(2));
[Y_dn,Y_dn_loc,Y_dn_f] = envelope2(Y_norm,thr(2));

sig1 = sig./ S_dn_f'*mean(S_dn_f); 
S1 = mean(sig1,2);
S1_norm = S1/prctile(S1,99)*prctile(V,99);
[S1_up,S1_up_loc,S1_up_f] = envelope2(S1_norm,thr(1));
[S1_dn,S1_dn_loc,S1_dn_f] = envelope2(S1_norm,thr(2));

%%
h = figure('Name',strrep(fn(1:end-1),'_','-'),'Position',[558 415 787 282]);
hold on;
plot(V,'Color',[0.5 .5 .5 .5]);
plot(Y_norm,'Color',[0 1 0 .3])
plot(S_norm,'Color',[1 0 0 .3])
plot(S1_norm,'Color',[0 0 1 .3])

plot(S_up_loc,S_up,'ro');
plot(S_up_f,'r--','Linewidth',2)
plot(S_dn_loc,S_dn,'ro');
plot(S_dn_f,'r--','Linewidth',2)

plot(Y_up_loc,Y_up,'go');
plot(Y_up_f,'g--','Linewidth',2)
plot(Y_dn_loc,Y_dn,'go');
plot(Y_dn_f,'g--','Linewidth',2)

plot(S1_up_loc,S1_up,'bo');
plot(S1_up_f,'b--','Linewidth',2)
plot(S1_dn_loc,S1_dn,'bo');
plot(S1_dn_f,'b--','Linewidth',2)


ylim([0 max(V)])
xlabel('Time (s)')
ylabel('Normalized response')
legend('Velocity','Noise\_norm','Signal\_norm','Signal1\_norm')
legend('Boxoff')
title(sprintf('Corr:sig-v=%.2f,noise-v=%.2f,sig1-v= %.2f ',...
    corr(S,V'),corr(Y_noise',V'),corr(S1,V')))
%


%%
saveas(h1,[fn 'noiseselection.fig'])
saveas(h,sprintf('%scasignal.fig',fn))
saveas(h,sprintf('%scasignal.png',fn))

save(sigfn,'Y_dn_f','sig1','-append');
if contains(sigfn,'&') % derive spike signals for multiple files
    fnames = splitfile(sigfn);
    try
        T=0;
        eachsize = size(sig,1)/numel(fnames);
        sig1ALL=sig1;
        for n=1:numel(fnames)
            sig1=sig1ALL(T+1:T+eachsize,:);
            save(fnames{n},'sig1','-append');
            T=T+eachsize;
        end
    catch
        disp('check maybe the file size is not evenly distributed')
    end
end