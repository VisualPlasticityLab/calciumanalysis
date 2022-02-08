function newfile=func_concat_signals()

f=uigetfile('*signals.mat','Pick a synced signals');
load(f);

j = menu('more file?','Yes','No');
while j==1
    fnew=uigetfile('*signals.mat','Pick another synced signals');
    dnew=load(fnew);
    % sigF:seg,rep,Var,ncell
    % matrix=rep*stim(Var)
    sigF=cat(2,sigF,dnew.sigF);
    matrix=cat(2,matrix,dnew.sigF);
    assert(Cor==dnew.Cor)
    j = menu('one more file?','Yes','No');
end

[~,SI.peakR,SI.errorR,~,SI.peakS,SI.errorS]= sigFcmp(sigF,sigwin,matrix);  % sigR:seg,1,Var,ncell

SI.OSI_R=calOSI(SI.peakR);
SI.OSI_S=calOSI(SI.peakS);
SI.gOSI_R=calgOSI(SI.peakR);
SI.gOSI_S=calgOSI(SI.peakS);
SI.DSI_R=calDSI(SI.peakR);
SI.DSI_S=calDSI(SI.peakS);

hpk(1)=figure('Name','OSIcomparison');hold on;
s1=scatter(SI.gOSI_S(:),SI.gOSI_R(:),'jitter','on', 'jitterAmount', 0.15);
s2=scatter(SI.OSI_S(:),SI.OSI_R(:),'g','jitter','on', 'jitterAmount', 0.15);
l=legend([s1 s2],'gOSI','OSI');
plot([0 max(SI.gOSI_R(:))],[0 max(SI.gOSI_R(:))],'--');
xlabel('gOSI still');ylabel('gOSI running');

hpk(2)=figure('Name','DSIcomparison');
scatter(SI.DSI_S(:),SI.DSI_R(:),'jitter','on', 'jitterAmount', 0.15);
hold on;
plot([0 max(SI.DSI_R(:))],[0 max(SI.DSI_R(:))],'--');
xlabel('DSI still');ylabel('DSI running');
%hpk(3)=polarplt(cat(1,SI.peakR,SI.peakS),Cor);

%% save figures and parameters
newfile=uiputfile('','Pick a name to save file','signals.mat');
folder='';
save(newfile,'sigF','matrix','sigwin','Cor','variant','SI');
for i=1:numel(hpk)
    savefig(hpk(i),fullfile(folder,hpk(i).Name));
    saveas(hpk(i),fullfile(folder,hpk(i).Name),'png');
end
