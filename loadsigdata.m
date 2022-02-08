function [eye1,eye2,sigtype]=loadsigdata

[~,ipsieyepath]=uigetfile('.mat','select PeakSI matfile for ipsieye/pre-expo');
[~,contraeyepath]=uigetfile('.mat','select PeakSI matfile for contraeye/post-expo');

variableInfo = who('-file',fullfile(ipsieyepath,'peakSI.mat'));
pos1 = strfind(ipsieyepath,filesep);
pos2 = strfind(contraeyepath,filesep);

if ismember('sigF',variableInfo)
    sigtype = 'ca';
    eye1=load(fullfile(ipsieyepath,'peakSI.mat'),'matrix','win_sig','sigF');
    eye2=load(fullfile(contraeyepath,'peakSI.mat'),'matrix','win_sig','sigF');

    [eye1.peakR,eye1.sigR,eye1.errorR,eye1.peakS,eye1.sigS,eye1.errorS]= sigFcmp(eye1.sigF,eye1.win_sig,eye1.matrix);
    [eye2.peakR,eye2.sigR,eye2.errorR,eye2.peakS,eye2.sigS,eye2.errorS]= sigFcmp(eye2.sigF,eye2.win_sig,eye2.matrix);
    [eye1.peak,eye1.sig,eye1.error]=cal_ER(eye1.sigF,eye1.win_sig);  %peak:1*Var*ncell
    [eye2.peak,eye2.sig,eye2.error]=cal_ER(eye2.sigF,eye2.win_sig);  %peak:1*Var*ncell
elseif ismember('sigspF',variableInfo)
    sigtype = 'sp';
    eye1=load(fullfile(ipsieyepath,'peakSI.mat'),'matrix','win_sig','sigspF');
    eye2=load(fullfile(contraeyepath,'peakSI.mat'),'matrix','win_sig','sigspF');

    [eye1.peakR,eye1.sigR,eye1.errorR,eye1.peakS,eye1.sigS,eye1.errorS]= sigspFcmp(eye1.sigspF,eye1.win_sig,eye1.matrix);
    [eye2.peakR,eye2.sigR,eye2.errorR,eye2.peakS,eye2.sigS,eye2.errorS]= sigspFcmp(eye2.sigspF,eye2.win_sig,eye2.matrix);
    [eye1.peak,eye1.sig,eye1.error]=cal_spER(eye1.sigspF,eye1.win_sig);  %peak:1*Var*ncell
    [eye2.peak,eye2.sig,eye2.error]=cal_spER(eye2.sigspF,eye2.win_sig);  %peak:1*Var*ncell
elseif ismember('sigMLF',variableInfo)
    sigtype = 'ML';
    eye1=load(fullfile(ipsieyepath,'peakSI.mat'),'matrix','win_sig','sigMLF');
    eye2=load(fullfile(contraeyepath,'peakSI.mat'),'matrix','win_sig','sigMLF');

    [eye1.peakR,eye1.sigR,eye1.errorR,eye1.peakS,eye1.sigS,eye1.errorS]= sigspFcmp(eye1.sigMLF,eye1.win_sig,eye1.matrix);
    [eye2.peakR,eye2.sigR,eye2.errorR,eye2.peakS,eye2.sigS,eye2.errorS]= sigspFcmp(eye2.sigMLF,eye2.win_sig,eye2.matrix);
    [eye1.peak,eye1.sig,eye1.error]=cal_spER(eye1.sigMLF,eye1.win_sig);  %peak:1*Var*ncell
    [eye2.peak,eye2.sig,eye2.error]=cal_spER(eye2.sigMLF,eye2.win_sig);  %peak:1*Var*ncell
end
