 function caanalysisplot_s(varargin)
if nargin<1
    [~,folder]=uigetfile('peakSI.mat');
    rewrite=1;
else
    folder = varargin{1};
    rewrite =0;
end
load(fullfile(folder,'peakSI.mat'));
folder = fullfile(folder,'s')
if ~exist(folder)
    mkdir(folder)
end
%% plot figures and calculate Orientation selection or Contrast invariant
%save([fname '_' num2str(numel(Cor)) 'signals.mat'],'sigsF','matrix','win_sig','Cor');
if ~exist('matrix','var') 
    matrix = run.matrix;
end

if ~exist('sigsF','var')
    variant=1;
    try
        sigsF = SI.sigsF;
        SI = rmfield(SI,'sigsF');
    catch
        [file,path]=uigetfile('*signals.mat','select sigsF file');
        signals=load(fullfile(path,file))
        matrix = signals.matrix;
        win_sig = signals.window;
        sigsF= signals.sigsF;
    end
end

sigF = cleandata(sigF,win_sig);
 
sigsF(isnan(sigsF))=0;
% Gd = [2 9 3 10 4 11 18 65 16 19 17 23 24 31 32 30 39 59 50 51 53 54 62 40 60 47 41 42 43];
% if exist('Gd','var')
%     Cor =Gd;
% end
Cor =1:size(sigsF,4);

% variant=11;
% drawingoption=11;
drawingoption=(~isempty(matrix)& sum(matrix(:)))*10+variant;

switch drawingoption
    case 0%'no running with contrast map'
        hsigsF=sigFplt(sigsF,[],win_sig,Cor);  % sigsF:seg,rep,Var,ncell
        [SI.peak,~,SI.error]=cal_ER(sigsF,win_sig);%peak: 1*Var*ncell
        
        hpk(1)=fitplt(SI.peak,SI.error,Cor);
        hpk(1).Name='Contrast tuning under best orientation';
        
    case 1 %'no running with orientation map'
        [SI.peak,~,SI.error]=cal_ER(sigsF,win_sig);  %peak:1*Var*ncell
        [SI.OSI,SI.pref_dir]=calOSI(SI.peak);
        SI.DSI=calDSI(SI.peak);
        SI.gOSI=calgOSI(SI.peak);
        
         hsigsF=sigFplt(sigsF,[],win_sig,Cor);  % sigsF:seg,rep,Var,ncell
        hpk(1)=polarplt(SI.peak,Cor);
        hpk(1).Name='ori tuning under best contrast';
        hpk(2)=fitplt(SI.peak,[],Cor);
        hpk(2).Name='ori tuning under best contrast';
        
        orientation = floor(max(SI.pref_dir)/2);
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
        %[hsigsF,Gd]=sigVarplt(sigsF,win_sig,[],Cor);  % sigsF:seg,rep,Var,ncell
        [SI.peak,~,SI.error]=cal_ER(SI.sigsF,win_sig);%peak:1*Var*ncell
        SI.peak=reshape(SI.peak,info.steps(2),info.steps(1),ncell);
        SI.error=reshape(SI.error,info.steps(2),info.steps(1),ncell);
        
        hpk(1)=fitplt(SI.peak,SI.error,Cor);           %peak: 1*Var*ncell
        hpk(1).Name='Orientation tuning under diff contrast-running';
        
        hpk(2)=polarplt(SI.peak,Cor);
        [SI.OSI,SI.pref_dir]=calOSI(SI.peak);
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
        [~,SI.peakR,SI.errorR,~,SI.peakS,SI.errorS]= sigFcmp(sigsF,win_sig,matrix);  % sigR:seg,1,Var,ncell
        hsigsF=sigFplt(sigsF,matrix,win_sig,Cor);  % sigsF:seg,rep,Var,ncell
        SI.peakR=fixpeak(SI.peakR);
        SI.peakS=fixpeak(SI.peakS);
        hpk(1)=figure('Name','Contrast tuning under best orientation-running');hold on;
        hpk(1)=fitplt(cat(1,SI.peakR,SI.peakS),cat(1,SI.errorR,SI.errorS),Cor);
        
    case 11  %'running with orientation'
        [SI.peakR,~,SI.errorR,SI.peakS,~,SI.errorS]= sigFcmp(sigsF,win_sig,matrix);  % sigR:seg,1,Var,ncell
%         SI.peakR=fixpeak(SI.peakR);
%         SI.peakS=fixpeak(SI.peakS);
        [SI.OSI_R,SI.pref_dir_R]=calOSI(SI.peakR);
        [SI.OSI_S,SI.pref_dir_S]=calOSI(SI.peakS);
        
        SI.gOSI_R=calgOSI(SI.peakR);
        SI.gOSI_S=calgOSI(SI.peakS);
        SI.DSI_R=calDSI(SI.peakR);
        SI.DSI_S=calDSI(SI.peakS);

        hsigsF=sigFplt(sigsF,matrix,win_sig,Cor);  % sigsF:seg,rep,Var,ncell
        hpk(1)=fitplt(cat(1,SI.peakR,SI.peakS),cat(1,SI.errorR,SI.errorS),Cor);
               
        hpk(2)=figure('Name','OSIcomparison');hold on;
        scatter(SI.gOSI_S(:),SI.gOSI_R(:),'b','jitter','on', 'jitterAmount', 0.15);
        scatter(SI.OSI_S(:),SI.OSI_R(:),'g','jitter','on', 'jitterAmount', 0.15);
        plot([0 max(SI.gOSI_R(:))],[0 max(SI.gOSI_R(:))],'--');
        xlabel('OSI still');ylabel('OSI running');
        legend('gOSI','OSI')
        
        hpk(3)=figure('Name','DSIcomparison');
        scatter(SI.DSI_S(:),SI.DSI_R(:),'jitter','on', 'jitterAmount', 0.15);
        hold on;
        plot([0 max(SI.DSI_R(:))],[0 max(SI.DSI_R(:))],'--');
        xlabel('DSI still');ylabel('DSI running');
        
        orientation = floor(max(SI.pref_dir_R)/2);
        SI.pref_ori_R=SI.pref_dir_R;
        SI.pref_ori_R(SI.pref_ori_R>orientation)=SI.pref_ori_R(SI.pref_ori_R>orientation)-orientation;
        SI.pref_ori_S=SI.pref_dir_S;
        SI.pref_ori_S(SI.pref_ori_S>orientation)=SI.pref_ori_S(SI.pref_ori_S>orientation)-orientation;

        hpk(4)=figure('Name','Preferred Orientation');
        histplt2(SI.pref_ori_R,SI.pref_ori_S);
        xlabel('Orientation');ylabel('cell count');
        
        hpk(5)=figure('Name','Preferred Direction');
        histplt2(SI.pref_dir_R,SI.pref_dir_S);
        xlabel('Direction');ylabel('cell count');
    case 12%'running with both contrast and orientation'
        hsigsF=sigVarplt(sigsF,win_sig,matrix,Cor);  % sigsF:seg,rep,Var,ncell
        [~,SI.peakR,SI.errorR,~,SI.peakS,SI.errorS]= sigFcmp(sigsF,win_sig,matrix);  % calculate peak from prestim:prestim+stimON-1
        SI.peakR=reshape(SI.peakR,info.steps(2),info.steps(1),ncell);
        SI.peakR=fixpeak(SI.peakR);
        SI.errorR=reshape(SI.errorR,info.steps(2),info.steps(1),ncell);
        
        
        SI.peakS=reshape(SI.peakS,info.steps(2),info.steps(1),ncell);
        SI.peakS=fixpeak(SI.peakS);
        SI.errorS=reshape(SI.errorS,info.steps(2),info.steps(1),ncell);
        
        hpk(1)=fitplt(cat(1,SI.peakR,SI.peakS),cat(1,SI.errorR,SI.errorS),Cor);
        hpk(1).Name='Orientation tuning under diff contrast-running';
        hpk(2)=polarplt(cat(1,SI.peakR,SI.peakS),Cor);
        % SI.peakR contrast*ori*ncell
        
        [SI.OSI_R,SI.pref_dir_R]=calOSI(SI.peakR);
        [SI.OSI_S,SI.pref_dir_S]=calOSI(SI.peakS);
        
        SI.gOSI_R=calgOSI(SI.peakR);
        SI.gOSI_S=calgOSI(SI.peakS);
        
        hpk(3)=fitplt(cat(1,SI.gOSI_R,SI.gOSI_S),[],Cor);
        hpk(3).Name='gOSI under diff contrast and running state';
        
end
%     answer = inputdlg('Enter gd cells','Pick gd cells',[1 80]);
% if ~isempty(answer)
%         Gd=str2num(answer{1});  rewrite =1;
% elseif menu('rewrite?','Yes','No')==1
%     rewrite=1;
% end

%% rewrite some parameters
if rewrite
%     baselinewin=win_sig(1)-5:win_sig(1);
%     baseline=nanmean(sigsF(baselinewin,:,:,:));
%     seg=size(sigsF,1);
%     sigsF=sigsF-repmat(baseline,[seg 1 1 1 ]);
%     threshold=2*nanstd(sigsF(baselinewin,:,:,:));% threshold size: 1,1,Var,ncell
%     sigsF(repmat(mean(sigsF(win_sig(1):win_sig(end),:,:,:))<threshold,[seg 1 1 1 ]))=NaN; % set
%     [b,a] = butter(12,.2,'low');           % IIR filter design
%     sigsF = filtfilt(b,a,sigsF);
    save(fullfile(folder,'peakSI.mat'),'SI','sigsF','variant','Cor','matrix','win_sig');
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
else
    save(fullfile(folder,'peakSI.mat'),'SI','sigsF','variant','Cor','matrix','win_sig');
    try
        for i=1:numel(hpk)
        savefig(hpk(i),fullfile(folder,hpk(i).Name));
        saveas(hpk(i),fullfile(folder,hpk(i).Name),'png');
        end
        for i=1:numel(hsigsF)
        saveas(hsigsF(i),fullfile(folder,hsigsF(i).Name),'png');
        end
    end
 end
close all


    
