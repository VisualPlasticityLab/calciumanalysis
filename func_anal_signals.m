clear all;close all;clc
% func_sync_signals
% func_sync_signals

figname=uigetfile('*.fig','find fig with the labeled cells');
open(figname);

contraeye=uigetfile('*signals.mat','open signals for contraeye');
contra=load(contraeye)

ipsieye=uigetfile('signals.mat','open signals for ipsieye');
ipsi=load(ipsieye)

%% running with orientation
% sigF:seg,rep,Var,ncell
%matrix:rep,Var
sigF=cat(3,contra.sigF,ipsi.sigF);
matrix=cat(2,contra.matrix,ipsi.matrix);
sigwin=intersect(contra.sigwin,ipsi.sigwin);
Cor=contra.Cor;

[hsigF,Gd]=sigFplt(sigF,matrix,sigwin,Cor);
%% calculate ODI using preferred orientation
 
ncells = size(sigF,4);
nSteps = size(sigF,3);

%SI.peakR:1*Var*ncell
finalvalue.R=cat(1,contra.SI.peakR,ipsi.SI.peakR);
finalvalue.S=cat(1,contra.SI.peakS,ipsi.SI.peakS);

if mod(nSteps,2)==0
    response=finalvalue.R;
else
    response=finalvalue.R(:,1:end-1,:);
end

[~,pref_contra]=max(response(1,:,:),[],2); %contralateral
[~,pref_ipsi]=max(response(2,:,:),[],2);
response_comb=(response(1,:,:)+response(2,:,:))/2;
[response_pref,pref]=max(response_comb,[],2);
pref=squeeze(pref);
pref_contra=squeeze(pref_contra);
pref_ipsi=squeeze(pref_ipsi);


contra.SI.peakR=squeeze(contra.SI.peakR);
ipsi.SI.peakR=squeeze(ipsi.SI.peakR);

pos=sub2ind(size(contra.SI.peakR),pref',1:ncells);
ODI_pref=(contra.SI.peakR(pos)-ipsi.SI.peakR(pos))./(ipsi.SI.peakR(pos)+contra.SI.peakR(pos));
ODI_pref_S=(contra.SI.peakS(pos)-ipsi.SI.peakS(pos))./(ipsi.SI.peakS(pos)+contra.SI.peakS(pos));


%% 

figure;hold on;
plot(abs(mod(pref_contra(logical(Gd))-pref(logical(Gd))+1,4)-1),'ro')
plot(abs(mod(pref_ipsi(logical(Gd))-pref(logical(Gd))+1,4)-1),'b+')
ylim([-.5 2.5]);
g=gca;
g.XLabel.String='Neurons';
g.YTick=[0 1 2];
g.YTickLabel={'0','\pi/4','\pi/2'};
g.YLabel.String='Difference of preferred orientation';
legend('contra','ipsi')
%saveas(gca,'preferred orientation.fig')

%% 

h=figure;
hold on;
seg=.3;
histogram(ODI_pref(logical(Gd)),[-1.2:.3:1.2]);
title(sprintf('total cell %d, using preferred orientation', sum(Gd)));
ODIfigname= [figname(1:end-4) 'good' num2str(sum(Gd)) '_ODI_prefori'];
% saveas(h,[ODIfigname '.fig']); saveas(h,[ODIfigname '.png']);
 save([ODIfigname '.mat'] ,'ODI_pref','Gd');
