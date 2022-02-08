function sigF=cleandata(sigF,win_sig)
%win_sig [prestim+1 prestim+stimON seg]

baseline=nanmin(sigF(win_sig(1)-5:win_sig(1)+5,:,:,:));

seg=size(sigF,1);
sigF=sigF-repmat(baseline,[seg 1 1 1 ]);
noisewin=win_sig(1)-5:win_sig(1)-1;
threshold=2*nanstd(sigF(noisewin,:,:,:));% threshold size: 1,1,Var,ncell
if numel(win_sig)<3
    sigF(repmat(mean(sigF(win_sig(1):win_sig(2),:,:,:))<threshold,[seg 1 1 1 ]))=0; % set
else
    sigF(repmat(mean(sigF(win_sig(1):win_sig(end),:,:,:))<threshold,[seg 1 1 1 ]))=0; % set
end