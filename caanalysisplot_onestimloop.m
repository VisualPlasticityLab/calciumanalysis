 function caanalysisplot_onestimloop(varargin)
 % caanalysisplot(foldername, Cor, figskip)
 % input: foldername
 % input: Cor, selected good cells
 % input: figskip, not plotting each cell TC
 % example caanalysisplot('.',f1.pair0&f1.pair1&pair2,0)
if nargin<1
    [~,folder]=uigetfile('peakSI.mat');
    rewrite=1;
else
    folder = varargin{1};
    rewrite =0;
end
load(fullfile(folder,'peakSI.mat'));

if nargin<2 
    Cor =1:size(sigF,4); % keep the last one even most times it'ss noise background
 % sigF: seg,rep,Var,ncell
else 
    Cor = varargin{2};
    if isempty(Cor) %
    Cor =1:size(sigF,4);
    end        
end  

if nargin<3
    figskip = 0;  % default to print all tuning Br changed on 120320 from 1
else
    figskip = varargin{3};
end

% sigF = sigF(:,:,:,Cor);
%% 
%save([fname '_' num2str(numel(Cor)) 'signals.mat'],'sigF','matrix','win_sig','Cor');
if ~exist('matrix','var') 
    matrix = run.matrix;
end

if ~exist('sigF','var')
    variant=1;
    try
        sigF = SI.sigF;
        SI = rmfield(SI,'sigF');
    catch
        [file,path]=uigetfile('*signals.mat','select sigF file');
        signals=load(fullfile(path,file))
        matrix = signals.matrix;
        win_sig = signals.window;
        sigF= signals.sigF;
    end
end
 
sigF(isnan(sigF))=0;

if win_sig(end) == size(sigF,1)
    temp=inputdlg(sprintf('segment=%d %d Manuel:',win_sig));
    win_sig = str2num(temp{1});
end
%% plot figures and calculate Orientation selection or Contrast invariant

wrunning = (~isempty(matrix)) & sum(matrix(:));
drawingoption=double(wrunning)*10+variant;
% variant=11;
% drawingoption=11;
% drawingoption=(~isempty(matrix)& sum(matrix(:)))*10+variant;

switch drawingoption
    case 0%'no running with contrast map'
        hsigF=sigFplt(sigF,[],win_sig,Cor,folder);  % sigF:seg,rep,Var,ncell
        [SI.peak,~,SI.error,cal_ER_option]=cal_ER_onestimloop(sigF,win_sig);%peak: 1*Var*ncell
        
        hpk(1)=fitplt(SI.peak,SI.error,Cor);
        hpk(1).Name='Contrast tuning under best orientation';
        
    case 1 %'no running with orientation map'
         hsigF=sigFplt(sigF,[],win_sig,Cor,folder);  % sigF:seg,rep,Var,ncell
       [SI.peak,~,SI.error,cal_ER_option]=cal_ER_onestimloop(sigF,win_sig);  %peak:1*Var*ncell
        [SI.OSI,SI.pref_dir]=calOSI(SI.peak);
        SI.DSI=calDSI(SI.peak);
        SI.gOSI=calgOSI(SI.peak);
        
%         hpk(1)=polarplt(SI.peak,Cor);
%         hpk(1).Name='ori tuning under best contrast';
% %         hpk(2)=fitplt(SI.peak,[],Cor);
%         hpk(2).Name='ori tuning under best contrast';
        
        orientation = floor(size(sigF,3)/2);
        SI.pref_ori=SI.pref_dir;
        SI.pref_ori(SI.pref_ori>orientation)=SI.pref_ori(SI.pref_ori>orientation)-orientation;
        hpk(3)=figure('Name','Preferred Orientation')
        h1 = histogram(SI.pref_ori);
        h1.FaceColor=[ 0.5 0.5 0.5];
        legend('anethtized','Location','north','Orientation','horizontal');
        xlabel('Orientation');ylabel('cell count');

        hpk(3)=figure('Name','Preferred Direction')
        h1 = histogram(SI.pref_dir);
        h1.FaceColor=[ 0.5 0.5 0.5];
        legend('anethtized','Location','north','Orientation','horizontal');
        xlabel('Orientation');ylabel('cell count');
        
    case 2 %'no running with both contrast&orientation map'
        %hsigF=sigVarplt(sigF,win_sig,[],Cor,folder);  % sigF:seg,rep,Var,ncell
if ~figskip
        hsigF=sigVarplt(sigF,info.steps,win_sig,matrix,Cor,folder);  % sigF:seg,rep,Var,ncell
end
        [SI.peak,~,SI.error,cal_ER_option]=cal_ER_onestimloop(sigF,win_sig);%peak:1*Var*ncell
%         SI.peak=reshape(SI.peak,info.steps(2),info.steps(1),[]);
%         SI.error=reshape(SI.error,info.steps(2),info.steps(1),[]);
%         
        hpk(1)=fitplt(SI.peak,SI.error,Cor);           %peak: 1*Var*ncell
        hpk(1).Name='Orientation tuning under diff contrast-running';
        
%         hpk(2)=polarplt(SI.peak,Cor);
        
        [SI.OSI,SI.pref_dir,hpk(2)]=calOSI(SI.peak);
        hpk(2).Name = 'Still only';
        SI.DSI=calDSI(SI.peak);
        SI.gOSI=calgOSI(SI.peak);
        hpk(3)=fitplt(SI.gOSI,[],Cor);
        hpk(3).Name='gOSI under diff contrast';
        
        hpk(4)=figure('Name','Preferred Orientation');
        h1 = histogram(SI.pref_dir);
        h1.FaceColor=[ 0.5 0.5 0.5];
        legend('anethtized','Location','north','Orientation','horizontal');
        xlabel('Orientation');ylabel('cell count');
    case 10  %'running with contrast'
        [~,SI.peakR,SI.errorR,~,SI.peakS,SI.errorS]= sigFcmp_onestimloop(sigF,win_sig,matrix);  % sigR:seg,1,Var,ncell
if ~figskip
        hsigF=sigFplt(sigF,matrix,win_sig,Cor,folder);  % sigF:seg,rep,Var,ncell
end
        SI.peakR=fixpeak(SI.peakR);
        SI.peakS=fixpeak(SI.peakS);
        hpk(1)=figure('Name','Contrast tuning under best orientation-running');hold on;
        hpk(1)=fitplt(cat(1,SI.peakR,SI.peakS),cat(1,SI.errorR,SI.errorS),Cor);
        
    case 11  %'running with orientation'    
if ~figskip
          hsigF=sigFplt(sigF,matrix,win_sig,Cor,folder);  % sigF:seg,rep,Var,ncell
end
    [SI.peakR,~,SI.errorR,SI.peakS,~,SI.errorS]= sigFcmp_onestimloop(sigF,win_sig,matrix);  % sigR:seg,1,Var,ncell
%         SI.peakR=fixpeak(SI.peakR);
%         SI.peakS=fixpeak(SI.peakS);
%         [SI.OSI_R,SI.pref_dir_R]=calOSI(SI.peakR);
%         [SI.OSI_S,SI.pref_dir_S]=calOSI(SI.peakS);
%         
%         SI.gOSI_R=calgOSI(SI.peakR);
%         SI.gOSI_S=calgOSI(SI.peakS);
%         SI.DSI_R=calDSI(SI.peakR);
%         SI.DSI_S=calDSI(SI.peakS);
% %matrix(:)=0;

%         hpk(1)=fitplt(cat(1,SI.peakR,SI.peakS),cat(1,SI.errorR,SI.errorS),Cor);
%         hpk(2)=figure('Name','OSIcomparison');
%         subplot(1,2,1);hold on;
%         scatter(SI.gOSI_S(:),SI.gOSI_R(:),'b','jitter','on', 'jitterAmount', 0.15);
%         scatter(SI.OSI_S(:),SI.OSI_R(:),'g','jitter','on', 'jitterAmount', 0.15);
%         plot([0 max(SI.gOSI_R(:))],[0 max(SI.gOSI_R(:))],'--');
%         xlabel('OSI still');ylabel('OSI running');
%         legend('gOSI','OSI')
%         subplot(1,2,2);hold on;
%         scatter(SI.DSI_S(:),SI.DSI_R(:),'jitter','on', 'jitterAmount', 0.15);
%         plot([0 max(SI.DSI_R(:))],[0 max(SI.DSI_R(:))],'--');
%         xlabel('DSI still');ylabel('DSI running');
        
%         orientation = floor(size(sigF,3));
%         SI.pref_ori_R=SI.pref_dir_R;
%         SI.pref_ori_R(SI.pref_ori_R>orientation)=SI.pref_ori_R(SI.pref_ori_R>orientation)-orientation;
%         SI.pref_ori_S=SI.pref_dir_S;
%         SI.pref_ori_S(SI.pref_ori_S>orientation)=SI.pref_ori_S(SI.pref_ori_S>orientation)-orientation;
% 
%         hpk(3)=figure('Name','Preferred Orientation & Direction');
%         subplot(1,2,1);histplt2(SI.pref_ori_S,SI.pref_ori_R,.5:1:orientation+1.5);
%         legend({'still','run'},'orientation','horizontal')
%         xlabel('Orientation');ylabel('cell count');
%         subplot(1,2,2);histplt2(SI.pref_dir_S,SI.pref_dir_R,.5:1:orientation*2+1.5);
%         legend({'still','run'},'orientation','horizontal')
%         xlabel('Direction');ylabel('cell count');
%         
        
        [SI.peak,~,SI.error,cal_ER_option]=cal_ER_onestimloop(sigF,win_sig);  %peak:1*Var*ncell
        %[SI.OSI,SI.pref_dir]=calOSI(SI.peak);
        %SI.gOSI=calgOSI(SI.peak);
        %SI.DSI=calDSI(SI.peak);
        orientation = floor(size(sigF,3)/2);
        %SI.pref_ori=SI.pref_dir;
        %SI.pref_ori(SI.pref_ori>orientation)=SI.pref_ori(SI.pref_ori>orientation)-orientation;
%         hpk(4)=figure('Name','Preferred Orientation&Direction')
%         subplot(1,2,1);h1 = histogram(SI.pref_ori);
%         h1.FaceColor=[ 0.5 0.5 0.5];
%         legend('anethtized','Location','north','Orientation','horizontal');
%         xlabel('Orientation');ylabel('cell count');
%         subplot(1,2,2);h1 = histogram(SI.pref_dir);
%         h1.FaceColor=[ 0.5 0.5 0.5];
%         legend('anethtized','Location','north','Orientation','horizontal');
%         xlabel('Orientation');ylabel('cell count');
    case 12 %'running with both orientation and contrast/eye/etc'
if ~figskip
        hsigF=sigVarplt(sigF,info.steps,win_sig,matrix,Cor,folder);  % sigF:seg,rep,Var,ncell
end          
        [SI.baseR,SI.SdR,SI.peakR,SI.errorR,~,...
            SI.baseS,SI.SdS,SI.peakS, SI.errorS,~, cal_ER_option]= sigFcmp_onestimloop(sigF,win_sig,matrix);  % calculate peak from prestim:prestim+stimON-1
%         SI.peakR=reshape(SI.peakR,info.steps(2),info.steps(1),[]);
%         SI.errorR=reshape(SI.errorR,info.steps(2),info.steps(1),[]);
%         
%         SI.peakS=reshape(SI.peakS,info.steps(2),info.steps(1),[]);
%         SI.errorS=reshape(SI.errorS,info.steps(2),info.steps(1),[]);
%         
%         hpk(1)=fitplt(cat(1,SI.peakR,SI.peakS),cat(1,SI.errorR,SI.errorS),Cor);
%         hpk(1).Name='Orientation tuning under various condition-running';
%        
%        [SI.OSI_R,SI.gOSI_R,SI.DSI_R,SI.gDSI_R,SI.pref_dir_R]=calOSI(SI.peakR);
        %hpk(2).Name = 'Running condition';
%        [SI.OSI_S,SI.gOSI_S,SI.DSI_S,SI.gDSI_S,SI.pref_dir_S]=calOSI(SI.peakS);
        %hpk(3).Name = 'Still condition';
        
%        SI.gOSI_R=calgOSI(SI.peakR);
%        SI.gOSI_S=calgOSI(SI.peakS);

%         hpk(4)=fitplt(cat(1,SI.gOSI_R,SI.gOSI_S),[],Cor);
%         hpk(4).Name='gOSI under diff contrast and running state';
end
%% save/rewrite some figure & parameters
%folder = fullfile(folder,cal_ER_option);
% if ~exist(folder)
%     mkdir(folder)
% end

% if rewrite
%     baselinewin=win_sig(1)-5:win_sig(1);
%     baseline=nanmean(sigF(baselinewin,:,:,:));
%     seg=size(sigF,1);
%     sigF=sigF-repmat(baseline,[seg 1 1 1 ]);
%     threshold=2*nanstd(sigF(baselinewin,:,:,:));% threshold size: 1,1,Var,ncell
%     sigF(repmat(mean(sigF(win_sig(1):win_sig(end),:,:,:))<threshold,[seg 1 1 1 ]))=NaN; % set
%     [b,a] = butter(12,.2,'low');           % IIR filter design
%     sigF = filtfilt(b,a,sigF);
%     save(fullfile(folder,'peakSI.mat'),'SI','sigF','variant','Cor','matrix','win_sig');
    %     picked= dir('*pick2ed*.fig');
    %     pickedimg=picked.name(1:end-4);
    %     cd ../..;
    %     img=openfig([strtok(pickedimg,'pi') '.fig']);
    %     axesObjs = get(img, 'Children');  %axes handles
    %     dataObjs = get(axesObjs, 'Children');
    %     ncell =( numel( dataObjs)-1)/2;
    %     bad = setdiff(1:ncell,Cor);
    %     temp=[dataObjs(ncell+1-bad),dataObjs(2*ncell+1-bad)];
    %     delete(temp);
    %     title(strrep(pickedimg,'_','\_'),'FontSize',16);
    %     saveas(img, fullfile(folder,[pickedimg '.fig']));
    %     saveas(img, fullfile(folder,[pickedimg '.png']));
    %     caanalysisplot(folder);
% else
    if ~exist(fullfile(folder,'peakSI.mat'))
        save(fullfile(folder,'peakSI.mat'),'SI','sigF','variant','Cor','matrix','win_sig','-v7.3');
    else
        save(fullfile(folder,'peakSI.mat'),'SI','-append'); %need to update SI calculation JSUN
    end
    if exist('hpk','var')
        for i=1:numel(hpk)
            try
        %         savefig(hpk(i),fullfile(folder,hpk(i).Name));
                saveas(hpk(i),fullfile(folder,hpk(i).Name),'png');
                end
            end
    end
%  end
%close all


    
