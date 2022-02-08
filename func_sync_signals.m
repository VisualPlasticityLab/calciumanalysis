function newdatafile=func_sync_signals()
% load calcium signal and sync with the ball motion/ pupil dilation data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. load calcium signal
% 2. syn with external stimulus
% 3. syn with internal state: ball motion/ eye motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global info;
%% load calcium signal and open recorded image
[f,p]=uigetfile('.signals','load calcium signal data .signals');
cd(p);
load(f,'-mat');
if strfind(f,'pick')
    load([f(1:strfind(f,'pick')-1) '.signals'], '-mat');
end
if strfind(f,'_memmap')
    fname=f(1:strfind(f,'_memmap')-1)
end
if strfind(f,'x')
    temp=uigetfile('*.sbx','Pick the correspondent recording file');
    fname=strtok(temp,'.');
end

if ~exist('sig')
    sig=sig_chunk;
end
if size(sig,1)<size(sig,2)
    sig=sig';
end
total=size(sig,2);

if ~exist('Cor')
    Cor=1:total;
else
    bad=setdiff(1:total,Cor);
end

sig=sig(:,Cor);
%% sync with stimulus data and plot

movement=preprocess(fname);
CAframeHz =info.resfreq/info.recordsPerBuffer;
CAlineHz = info.config.lines;
stim=zeros(1,info.max_idx); %just for plotting purpose,assign stimtype to all the ON frames
for i=1:(numel(info.frame)/2)
    stim(info.frame(i*2-1):info.frame(i*2))=info.stimtype(i);
end

if movement
    try
        load([fname '.ball'],'-mat');
    catch
        load([fname '_ball'],'-mat');
        speed=abs(ball(1:end-1))/192*2./diff(time)'; %size of the imaging area 192pixel, 2cm
    end
    % load([fname '.eye'],'-mat');
    sp=conv(speed,ones(1,ceil(CAframeHz*1))/ceil(CAframeHz*1),'same'); %conv the speed into with a 1s filter
    velocity=interp1(0:1/2:info.max_idx, sp(1:info.max_idx*2+1),1:info.max_idx); % downsampling the speed into the frame num
else
    velocity=[];
end
hsig=sigplt(sig,[stim/max(stim);velocity/max(velocity)],Cor); %
% baseline=prctile(sig,20,1);
% sig=sig./repmat(baseline,size(sig,1),1)-1;

%% set paramenter for stimulus
ncell=size(sig,2);
prestim=ceil(CAframeHz*1); %use 1s baseline
stimON=info.frame(2)-info.frame(1); % StimON duration
seg=prestim+min(info.frame(3:2:end)-info.frame(1:2:(end-2))); % segment=prestim+info.frame(TTL ON+TTL off )

Var=numel(unique(info.stimtype));% calculate types of stimulus, orientations, contrast, etc
rep=floor(min(2*numel(info.stimtype),numel(info.frame))/2/Var); %calculate repetitions

sigT=zeros(seg*rep,Var,ncell); % build sigT based on signals from sig and organized according to segment
frame_on=zeros(Var,rep); %tag frame_on of each seg
sigwin=prestim:prestim+stimON;
%% reorganize stimulus according to stimulus type
last=info.frame(end-1)+seg-prestim;
if last>size(sig,1)
    sig(end+1:last,:)=0;
    disp(sprintf('padding zeros to the last %d frames', last-size(sig,1)));
end

for j=1:Var
    ori=find(info.stimtype==j);
    on=ori*2-1;
    frame_on(j,:)=info.frame(on(1:rep)); % the frame marking the prestim timepoint size 1*rep
    temp=ones(seg,1)*frame_on(j,:)+ (0:seg-1)'*ones(1,rep)-prestim; % array seg*rep6
    sigT(:,j,:)=sig(temp,:);   %(seg*rep)*Var*ncell
end

sigT=reshape(sigT,seg,rep,Var,ncell);
%% select cells based on their peak responses

if ~isempty(strfind(f,'memmap')) & isempty(strfind(f,'pick2'))
    hcell=openfig([f(1:findstr(f,'cell')+3) '.fig']);
    if ~exist('bad','var')
        bad=pick(f);
    end
    [bad2,hpick]=pick2(sigT,sigwin,Cor);
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
if movement
    run=[];%frame_on=zeros(Var,rep)
    for k=1:stimON
        run=cat(3,run,velocity(frame_on'+k-1));% run:rep,Var,stimON
    end
    [matrix,hrun]=runplt(run);
    %[matrix,hrun]=runplt2(run);
else
    matrix=[];
end
%% parse sigF according to different types of info.var

if info.steps(2)>1
    [ sigF,variant,hparse]=parseF(sigT,sigwin,matrix);

else
    sigF=sigT; % sigF:seg,rep,Var,ncell
    variant=strcmp(info.var(1),'Orientation');
end


%% plot figures and calculate Orientation selection or Contrast invariant
newdatafile=[fname '_' num2str(numel(Cor)) 'signals.mat'];
save(newdatafile,'sigF','matrix','sigwin','Cor','variant');

