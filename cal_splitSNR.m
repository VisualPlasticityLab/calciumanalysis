function [SNR,h] = cal_splitSNR

[eye1fn,p1]=uigetfile('.signals','select signal file for eye1eye');
[eye2fn,p2]=uigetfile('.signals','select signal file for eye2eye');

eye1sig = load(fullfile(p1,eye1fn),'-mat');
eye2sig = load(fullfile(p2,eye2fn),'-mat');

pos_ = strfind(eye1fn,'_');
eye1run = load(fullfile(p1,[eye1fn(1:pos_(3)) 'ball.mat']),'velocity');
eye2run = load(fullfile(p2,[eye2fn(1:pos_(3)) 'ball.mat']),'velocity');
eye1stim = load(fullfile(p1,eye1fn(1:pos_(3)-1)),'info');
eye2stim = load(fullfile(p2,eye2fn(1:pos_(3)-1)),'info');

SNR.skew = max( [skewness(eye1sig.sig); skewness(eye2sig.sig)]);

mmpeaks1 = prctile(eye1sig.sig,[5 50 99]); %;eye2sig.sig
mmpeaks2 = prctile(eye2sig.sig,[5 50 99]); %;eye2sig.sig
SNR.ca = max((mmpeaks1(end,:)+1)./(mmpeaks1(1,:)+1),(mmpeaks2(end,:)+1)./(mmpeaks2(1,:)+1));
% SNR = (mmpeaks(3,:)-mmpeaks(1,:))./mmpeaks(2,:);

mmpeaks_sp = prctile([eye1sig.sigsp;eye2sig.sigsp],[5 50 99]); %
SNR.sp = (mmpeaks_sp(end,:)+1)./(mmpeaks_sp(1,:)+1);
% SNR.sp = (mmpeaks_sp(3,:)-mmpeaks_sp(1,:))./mmpeaks_sp(2,:);

[SNRsorted,Cor] = sort(SNR.ca,'descend');
selected = Cor(SNRsorted>=1.5);
h(1)=figure;
hold on;
scatter(SNR.ca,SNR.sp,'ko')
scatter(SNR.ca(selected),SNR.sp(selected),'ro')
xlabel('Calcium signal noise ratio 99%-5%')
ylabel('Spikie signal noise ratio 99.5% - 50%');

h(2) = sigplt_onepage([eye1sig.sigsp(:,selected);eye2sig.sigsp(:,selected)],[eye1run.velocity eye2run.velocity],...
    [eye1stim.info.frame;eye2stim.info.frame+eye1stim.info.max_idx+1],[eye1stim.info.stimtype; eye2stim.info.stimtype],selected)

h(3) = sigplt_onepage([eye1sig.sig(:,selected);eye2sig.sig(:,selected)],[eye1run.velocity eye2run.velocity],...
    [eye1stim.info.frame;eye2stim.info.frame+eye1stim.info.max_idx+1],[eye1stim.info.stimtype; eye2stim.info.stimtype],selected)

 saveas(h(2),'traces_ca.fig')
% saveas(h(3),'traces_sp.fig')

