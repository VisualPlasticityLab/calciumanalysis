
%function ca_analysis;
% load calcium signal and sync with the ball motion/ pupil dilation data, used for close-loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. load calcium signal
% 2. syn with external stimulus
% 3. syn with internal state: ball motion/ eye motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clear all; global info;

%% load calcium signal and open recorded image
[f,p]=uigetfile('.signals','load calcium signal data .signals');
[fname,path]=uigetfile('*.sbx','Pick the corresponding sbx file');
multiplane=3-menu('is it part of the multiplane data?','Yes-Part1','Yes-Part2','No');

load(fullfile(p,f),'-mat');
fn=fullfile(path,strtok(fname,'.'));
if strfind(f,'pick')
    load([f(1:strfind(f,'pick')-1) '.signals'], '-mat');
end
if ~exist('sig')
    sig=sig_chunk;
end
if size(sig,1)<size(sig,2)
    sig=sig';
end
if ~exist('Cor')
    Cor=1:size(sig,2)-1;
else
    bad=setdiff(1:total,Cor);
end
sig=sig(:,Cor);
if multiplane
    part=3-multiplane;
    x=part:2:size(sig,1)*2;
    xx=1:size(sig,1)*2;
    sig=interp1(x,sig,xx,'spline');
end
%% sync with stimulus data and movement
movement=preprocess(fn);
CAframeHz =info.resfreq/info.recordsPerBuffer;
%just for plotting purpose,assign stimtype to all the ON frames
stim=zeros(1,info.max_idx+1); 
imagingsig=info.frame(1):info.frame(end);
stimsig= (imagingsig-info.frame(1))*(numel(info.sttype)-1)/(info.frame(end)-info.frame(1))+1;
stim(imagingsig)= info.sttype(round(stimsig));
% calculate types of stimulus, here orientations,[1 2 3 4 5]=[45 135 225 315 gray] 
% Var=numel(unique(info.stimtype));
Var=5;
% change to [4 3 2 1 0]=[45 135 225 315 gray];
stim = Var-stim;
grey = 0;
Var0= Var-1;

if size(stim,2)> size(sig,1)
    stim(size(sig,1)+1:end)=[];
end

if movement
    load([fn '_ball.mat'],'ball','time');
    speed=abs(ball(1:end-1))/192*2./diff(time)'; %size of the imaging area 192pixel, 2cm
    % load([fname '.eye'],'eye','-mat');
    sp=conv(speed,ones(1,ceil(CAframeHz*1))/ceil(CAframeHz*1),'same'); %conv the speed into with a 1s filter
    if numel(sp) >= info.max_idx*2
        velocity=interp1(0:1/2:info.max_idx, sp(1:info.max_idx*2+1),0:info.max_idx); % downsampling the speed into the frame num
 %    elseif multiplane
 %       velocity=interp1(1:(info.max_idx+1), sp,1:2:(info.max_idx+1));
    else
        velocity=sp;
    end
else
    velocity=[];
end
%% do some plotting
% figure;plot(stim);hold on;plot(info.frame,stim(info.frame),'*')
example=unique(randi(Cor(end),1,10));
hsig=sigplt(sig(:,example),[stim/max(stim);velocity/max(velocity)],Cor(example));
%% reorganize stimulus according to stimulus type
ncell=size(sig,2);
rep0=0;
seg0=0;
Direction=menu('Analysis orientation or direction?', 'Ori','Direction');
stim_ori=stim;
if Direction==1
stim(stim==4)=2;
stim(stim==3)=1;
Var0=2;
end
figure;hold on;plot(velocity);plot(stim);plot(stim_ori);

for j=1:Var0
    ori_start{j}=find(stim==j&circshift(stim,[0,1])==grey);
    rep=numel(ori_start{j});
    if rep>=1
        oriends=find(stim==j&[diff(stim) 0]~=0);
        for i=1:rep
            temp= oriends(find(oriends>ori_start{j}(i),1));
            if ~isempty(temp)
                ori_end{j}(i)=temp;
            else
                ori_start{j}(i)=[];
            end
        end
        disp(sprintf('stimulus type %d has %d reptitions',j,rep));
        if rep>rep0;rep0=rep;end
        if max(ori_end{j}-ori_start{j})>seg0; seg0=max(ori_end{j}-ori_start{j});end
        plot(ori_start{j},j*ones(size(ori_start{j})),'*');
        plot(ori_end{j},j*ones(size(ori_end{j})),'*');
    end

end
prestim = ceil(CAframeHz*1);
sigT = nan(prestim+seg0,rep0,Var0,ncell); % build sigT based on signals from sig and organized according to segment
duration= nan(rep0,Var0);
for j = 1:Var0
    for i = 1:numel(ori_start{j})
        duration(i,j)= ori_end{j}(i)-ori_start{j}(i)+1; 
        sigT(1:duration(i,j)+prestim,i,j,:) = sig(ori_start{j}(i)-prestim:ori_end{j}(i),:);
    end
end

%% select cells based on their peak responses

if ~isempty(strfind(f,'memmap')) & isempty(strfind(f,'pick2'))
    hcell=openfig([f(1:findstr(f,'cell')+3) '.fig']);
    if ~exist('bad','var')
        bad=pick(f);
    end
    [bad2,hpick]=pick2(sigT,window,Cor);
    bad=unique([bad bad2]);
    if ~isempty(bad)
        Cor = setdiff(1:ncell,bad);
        ncell= numel(Cor);
        newf=[f(1:end-8) 'pick2ed' num2str(ncell)];
        
        axesObjs = get(hcell, 'Children');  %axes handles
        dataObjs = get(axesObjs, 'Children');
        temp=dataObjs(total+1-bad);
        delete(temp);
        savefig(hcell,[newf '.fig']);
        
        sig=sig(:,Cor);
        sigT=sigT(:,:,:,Cor);
        save([newf '.signals'],'Cor');
    end
end

%% syn with ball/eye motion
% if movement
%     run=[];%frame_on=zeros(Var,rep)
%     for k=1:seg0
%         run=cat(3,run,velocity(frame_on'+k-1));% run:rep,Var,stimON
%     end
%     [matrix,hrun]=runplt(run);
%     %[matrix,hrun]=runplt2(run);
% else
    matrix=[];
% end
%% parse sigF according to different types of info.var

if info.steps(2)>1
    [ sigF,variant,hparse]=parseF(sigT,window,matrix);
else
    sigF=sigT; % sigF:seg,rep,Var,ncell
    variant=strcmp(info.var(1),'Orientation');
end
%% plot figures and calculate Orientation selection or Contrast invariant
%save([fname '_' num2str(numel(Cor)) 'signals.mat'],'sigF','matrix','window','Cor');

drawingoption=(~isempty(matrix))*10+variant;
win_sig = prestim:prestim+nanmedian(duration(:));
switch drawingoption
    case 0%'no running with contrast map'
        hsigF=sigFplt(sigF,[],win_sig,Cor);  % sigF:seg,rep,Var,ncell
        [SI.peak,~,SI.error]=cal_ER(sigF,win_sig);%peak: 1*Var*ncell
        
        hpk(1)=fitplt(SI.peak,SI.error,Cor);
        hpk(1).Name='Contrast tuning under best orientation';
        
    case 1 %'no running with orientation map'
        [SI.peak,~,SI.error]=cal_ER(sigF,win_sig);  %peak:1*
        hsigF=sigFplt(sigF,[],win_sig,Cor);  % sigF:seg,rep,Var,ncell
        
        hpk(1)=polarplt(SI.peak,Cor);
        hpk(1).Name='ori tuning under best contrast';
        
        SI.OSI=calOSI(SI.peak);
        SI.DSI=calDSI(SI.peak);
        SI.gOSI=calgOSI(SI.peak);
        
        hpk(2)=figure('Name','orientation tuning under best contrast');hold on;
        fitplt(SI.peak,[],Cor);
        
    case 2 %'no running with both contrast&orientation map'
        [hsigF,Gd]=sigVarplt(sigF,win_sig,[],Cor);  % sigF:seg,rep,Var,ncell
        [SI.peak,~,SI.error]=cal_ER(sigF,win_sig);%peak:1*Var*ncell
        SI.peak=reshape(SI.peak,info.steps(2),info.steps(1),ncell);
        SI.error=reshape(SI.error,info.steps(2),info.steps(1),ncell);
        
        hpk(1)=fitplt(SI.peak,SI.error,Cor);           %peak: 1*Var*ncell
        hpk(1).Name='Orientation tuning under diff contrast-running';
        
        hpk(2)=polarplt(SI.peak,Cor);
        
        SI.OSI=calOSI(SI.peak);
        SI.DSI=calDSI(SI.peak);
        SI.gOSI=calgOSI(SI.peak);
        
        hpk(3)=fitplt(SI.gOSI,[],Cor);
        hpk(3).Name='gOSI under diff contrast';
        
    case 10  %'running with contrast'
        [~,SI.peakR,SI.errorR,~,SI.peakS,SI.errorS]= sigFcmp(sigF,win_sig,matrix);  % sigR:seg,1,Var,ncell
        hsigF=sigFplt(sigF,matrix,win_sig,Cor);  % sigF:seg,rep,Var,ncell
        SI.peakR=fixpeak(SI.peakR);
        SI.peakS=fixpeak(SI.peakS);
        
        
        hpk(1)=figure('Name','Contrast tuning under best orientation-running');hold on;
        hpk(1)=fitplt(cat(1,SI.peakR,SI.peakS),cat(1,SI.errorR,SI.errorS),Cor);
        
        
    case 11  %'running with orientation'
        hsigF=sigFplt(sigF,matrix,win_sig,Cor);  % sigF:seg,rep,Var,ncell
        [~,SI.peakR,SI.errorR,~,SI.peakS,SI.errorS]= sigFcmp(sigF,win_sig,matrix);  % sigR:seg,1,Var,ncell
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
        hsigF=sigVarplt(sigF,win_sig,matrix,Cor);  % sigF:seg,rep,Var,ncell
        [~,SI.peakR,SI.errorR,~,SI.peakS,SI.errorS]= sigFcmp(sigF,win_sig,matrix);  % calculate peak from prestim:prestim+stimON-1
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

mark=findstr(fn,'_');
date=fn(mark(1)+1:mark(2)-1);
serial=fn(mark(2)+1:end);
if multiplane
    folder=fullfile([ date '-' serial '-' num2str(3-multiplane) ],['picked' num2str(numel(Cor))],['drawingoption' num2str(drawingoption)]);
else
    folder=fullfile([ date '-' serial ],['picked' num2str(numel(Cor))],['drawingoption' num2str(drawingoption)]);
end
if Direction==1
    folder = [ folder '_dir'];
else
    folder = [ folder '_ori'];
end

if exist(folder,'dir')
    if menu('Existing','Override','New')==2
           folder = [folder '_2'];
    end
end
   mkdir(folder);


SI.sigF = sigF;
if ~exist('Gd')
    Gd=Cor;
end
save(fullfile(folder,'peakSI.mat'),'SI','Gd','matrix','win_sig','Cor');

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
close all;

