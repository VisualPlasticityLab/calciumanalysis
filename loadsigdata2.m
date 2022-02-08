function [eye1,eye2,sigtype]=loadsigdata2
d=dir('peakSI.mat');
if ~isempty(d)
    botheyepath=pwd;
else
    [~,botheyepath]=uigetfile('.mat','select PeakSI matfile for both eye');
end
variableInfo = who('-file',fullfile(botheyepath,'peakSI.mat'));
pos1 = strfind(botheyepath,filesep);

sigtype = 'sp';
load(fullfile(botheyepath,'peakSI.mat'),'matrix','win_sig','sigF');
eye1.sigF = sigF(:,:,2:2:end,:); % from stimfile ipsieye, starting=2
eye2.sigF = sigF(:,:,1:2:end-1,:);% from stimfile contraeye, starting=1

eye1.win_sig = win_sig;
eye2.win_sig = win_sig;

if isempty(matrix)
    matrix = logical(zeros(size(sigF,2),size(sigF,3)));
end

eye1.matrix = matrix(:,2:2:end);
eye2.matrix = matrix(:,1:2:end-1);

[eye1.baseR,eye1.baseSdR,eye1.peakR,eye1.errorR,eye1.peakAllR,...
    eye1.baseS,eye1.baseSdS,eye1.peakS,eye1.errorS,eye1.peakAllS]= sigFcmp(eye1.sigF,eye1.win_sig,eye1.matrix);
[eye2.baseR,eye2.baseSdR,eye2.peakR,eye2.errorR,eye2.peakAllR,...
    eye2.baseS,eye2.baseSdS,eye2.peakS,eye2.errorS,eye2.peakAllS]= sigFcmp(eye2.sigF,eye2.win_sig,eye2.matrix);
[eye1.base,eye1.baseSd,eye1.peak,eye1.error,eye1.peakAll]=cal_ER(eye1.sigF,eye1.win_sig);  %peak:1*Var*ncell
[eye2.base,eye2.baseSd,eye2.peak,eye2.error,eye2.peakAll]=cal_ER(eye2.sigF,eye2.win_sig);  %peak:1*Var*ncell

% Baseline adjusted

% [eye1.OSI_r,eye1.gOSI_r,eye1.DSI_r,eye1.gDSI_r,eye1.pref_ori_r] = calOSI(eye1.peakR-eye1.baseR);
% [eye2.OSI_r,eye2.gOSI_r,eye2.DSI_r,eye2.gDSI_r,eye2.pref_ori_r] = calOSI(eye2.peakR-eye2.baseR);
% 
% [eye1.OSI_s,eye1.gOSI_s,eye1.DSI_s,eye1.gDSI_s,eye1.pref_ori_s] = calOSI(eye1.peakS-eye1.baseS);
% [eye2.OSI_s,eye2.gOSI_s,eye2.DSI_s,eye2.gDSI_s,eye2.pref_ori_s] = calOSI(eye2.peakS-eye2.baseS);
% 
% [eye1.OSI,eye1.gOSI,eye1.DSI,eye1.gDSI,eye1.pref_ori] = calOSI(eye1.peak-eye1.base);
% [eye2.OSI,eye2.gOSI,eye2.DSI,eye2.gDSI,eye2.pref_ori] = calOSI(eye2.peak-eye2.base);

% Baseline no adjustment
[eye1.OSI_r,eye1.gOSI_r,eye1.DSI_r,eye1.gDSI_r,eye1.pref_ori_r] = calOSI(eye1.peakR);
[eye2.OSI_r,eye2.gOSI_r,eye2.DSI_r,eye2.gDSI_r,eye2.pref_ori_r] = calOSI(eye2.peakR);

[eye1.OSI_s,eye1.gOSI_s,eye1.DSI_s,eye1.gDSI_s,eye1.pref_ori_s] = calOSI(eye1.peakS);
[eye2.OSI_s,eye2.gOSI_s,eye2.DSI_s,eye2.gDSI_s,eye2.pref_ori_s] = calOSI(eye2.peakS);

[eye1.OSI,eye1.gOSI,eye1.DSI,eye1.gDSI,eye1.pref_ori] = calOSI(eye1.peak);
[eye2.OSI,eye2.gOSI,eye2.DSI,eye2.gDSI,eye2.pref_ori] = calOSI(eye2.peak);

