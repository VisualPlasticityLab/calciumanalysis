function func_concat_signals()

% sigF:seg,rep,Var,ncell
f=uigetfile('*signals.mat','Pick a synced signals');

j = menu('more file?','Yes','No');
while j==1
  temp=func_sync_signals();
  allextractedtraces=cat(2,allextractedtraces,temp);
  j = menu('one more file?','Yes','No');
end

newfile=uiputfile('','Pick a name to save file',[eyetype '_traces.mat']);
save(newfile,'SI','Gd');


drawingoption=(~isempty(matrix))*10+variant;

switch drawingoption
      
    case 11  %'running with orientation'
        [hsigF,Gd]=sigFplt(sigF,matrix,window,Cor);  
        [~,SI.peakR,SI.errorR,~,SI.peakS,SI.errorS]= sigFcmp(sigF,window,matrix);  % sigR:seg,1,Var,ncell
        hpk(1)=polarplt(cat(1,SI.peakR,SI.peakS),Cor);
        
        
        SI.OSI_R=calOSI(SI.peakR);        
        SI.OSI_S=calOSI(SI.peakS);        
        SI.gOSI_R=calgOSI(SI.peakR);
        SI.gOSI_S=calgOSI(SI.peakS);
        SI.DSI_R=calDSI(SI.peakR);
        SI.DSI_S=calDSI(SI.peakS);  
        
        hpk(2)=figure('Name','OSIcomparison');hold on;
        scatter(SI.gOSI_S(:),SI.gOSI_R(:),'jitter','on', 'jitterAmount', 0.15);
        scatter(SI.OSI_S(:),SI.OSI_R(:),'g','jitter','on', 'jitterAmount', 0.15);
        plot([0 max(SI.gOSI_R(:))],[0 max(SI.gOSI_R(:))],'--');
        xlabel('gOSI still');ylabel('gOSI running');
        
        hpk(3)=figure('Name','DSIcomparison');
        scatter(SI.DSI_S(:),SI.DSI_R(:),'jitter','on', 'jitterAmount', 0.15);
        hold on;
        plot([0 max(SI.DSI_R(:))],[0 max(SI.DSI_R(:))],'--');
        xlabel('DSI still');ylabel('DSI running');
        
        
        
    case 12%'running with both contrast and orientation'
        hsigF=sigVarplt(sigF,window,matrix,Cor);  % sigF:seg,rep,Var,ncell
        [~,SI.peakR,SI.errorR,~,SI.peakS,SI.errorS]= sigFcmp(sigF,window,matrix);  % calculate peak from prestim:prestim+stimON-1
        SI.peakR=reshape(SI.peakR,info.steps(2),info.steps(1),ncell);
        SI.peakR=fixpeak(SI.peakR);
        SI.errorR=reshape(errorR,info.steps(2),info.steps(1),ncell);
        
        
        SI.peakS=reshape(SI.peakS,info.steps(2),info.steps(1),ncell);
        SI.peakS=fixpeak(SI.peakS);
        SI.errorS=reshape(SI.errorS,info.steps(2),info.steps(1),ncell);
        
        hpk(1)=fitplt(cat(1,SI.peakR,SI.peakS),cat(1,SI.errorR,SI.errorS),Cor);
        hpk(1).Name='Orientation tuning under diff contrast-running';
        hpk(2)=polarplt(cat(1,SI.peakR,SI.peakS),Cor);
        % SI.peakR contrast*ori*ncell
        
        SI.OSI_R=calOSI(SI.peakR);        
        SI.OSI_S=calOSI(SI.peakS);
        
        SI.gOSI_R=calgOSI(SI.peakR);
        SI.gOSI_S=calgOSI(SI.peakS);
        
        hpk(3)=fitplt(cat(1,SI.gOSI_R,SI.gOSI_S),[],Cor);
        hpk(3).Name='gOSI under diff contrast and running state';
        
end

%% save figures and parameters
% prompt = 'Do you want more? Y/N [Y]: ';
% str = input(prompt,'s');
% if isempty(str)
%     str = 'Y';
% end


for i=1:numel(hsigF)
    saveas(hsigF(i),fullfile(folder,hsigF(i).Name),'png');
end
for i=1:numel(hpk)
    savefig(hpk(i),fullfile(folder,hpk(i).Name));
    saveas(hpk(i),fullfile(folder,hpk(i).Name),'png');
end
if exist('hrun')
for i=1:numel(hrun)
    saveas(hrun(i),fullfile(folder,hrun(i).Name),'png');
end
end

if exist('hparse')c
    for i=1:numel(hparse)
    saveas(hparse(i),fullfile(folder,hparse(i).Name),'png');
    end
end
if exist('hpick')
    saveas(hpick,fullfile(folder,hpick.Name),'png');
end