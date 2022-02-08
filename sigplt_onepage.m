function hsig=sigplt_onepage(sig,velocity,stimtime,stimtype,Cor)


ncell=size(sig,2);
nframe=size(sig,1);
%% conversion of velcotiy and stimtime to be compatiable with sig

totalframes = numel(velocity);
velocity = downsample(velocity,round(totalframes/nframe));
v_cutoff = prctile(velocity,[99.9 99.99]);
if v_cutoff(2)/v_cutoff(1)>2
    velocity(velocity>v_cutoff(1))=v_cutoff(1);
end
%% sort by correlation with running
% if ~exist('Cor')
%     mmpeaks = prctile(sig,[5 50 99]);
% %     [SNR,Cor] = sort((mmpeaks(3,:)-mmpeaks(2,:)+1)./(mmpeaks(2,:)-mmpeaks(1,:)+1));
%     [SNR,Cor] = sort((mmpeaks(3,:)+1)./(mmpeaks(1,:)+1),'descend');
%     Cor = Cor(SNR>1.25);
%     sig = sig(:,Cor);
% else
%     mmpeaks = prctile(sig,[5 50 99]);
% %   [SNR,Cor] = sort((mmpeaks(3,:)-mmpeaks(2,:)+1)./(mmpeaks(2,:)-mmpeaks(1,:)+1));
%     SNR= (mmpeaks(3,:)+1)./(mmpeaks(1,:)+1);
% end
% for nth=1:ncell
%      SNR(nth) = corr(sig(:,nth),velocity(1:nframe)');
% end
    SNR = skewness(sig);

if ~exist('Cor')
    [~,I]= sort((SNR),'descend');
    SNR = SNR(I);
    sig = sig(:,I);
    Cor = I;
end
%% PLOT
hsig=figure('Name',['CA sig w stim_one page' ], 'position',[50 50 1500 50*ncell]);
hold on;

for nth=1:ncell
%         y =  y/(max(y)-min(y));
        plot(1:nframe,sig(:,nth)/max(sig(:,nth))*2+nth,'LineWidth',1,'color','k');   
        text(nframe*1.00,nth,sprintf('#%d SNR=%.2f',Cor(nth),SNR(nth)),...
            'HorizontalAlignment','left','VerticalAlignment','bottom');
end
plot(1:nframe,(velocity(1:nframe)')/max(velocity)+ncell+2,'LineWidth',3,'Color',[0 0 0 .3]);
text(nframe*1.00,ncell+2,['Run speed'],'HorizontalAlignment','left');
% colorscheme = hsv2rgb([(1:nstim)'/nstim, ones(nstim,1), ones(nstim,1)]);
% xaxis = find(nstim>0);
% yaxis = nstim(nstim>0);

%
%% color scheme for stimulus
if exist('stimtime')
    if numel(stimtime) >= 2*numel(stimtype)
        stimend = stimtime(2:2:end);
        stimtime = stimtime(1:2:end);
    else
        stimend = stimtime+3;
    end
    stimtime = round(stimtime/round(totalframes/nframe));
    stimend = round(stimend/round(totalframes/nframe));
    nstim = floor(numel(unique(stimtype))/2);
    clut = 1-hsv2rgb([(1:nstim)'/nstim, ones(nstim,1), ones(nstim,1)]);
    clut(nstim+1:nstim*2,:) = clut(1:nstim,:)/1.2;
    clut(end+1,:) = [.1 .1 .1];

    % if ~exist('sigsp')
    %     for i=1:ncell
    %         [sigs(:,i),sigsp(:,i)]=deconvolveCa(sig(:,i));
    %     end
    % end
    for ii=1:nstim*2+1
    %     plot([stimtime(stimtype==ii) stimtime(stimtype==ii)],[0 nth+2],'Color',[clut(ii,:) .5]);
        stimpos = find(stimtype==ii);
        x = cat(2,stimtime(stimpos),stimend(stimpos),...
            stimend(stimpos),stimtime(stimpos));
        y = repmat([1 1 ncell+6 ncell+6],[numel(stimpos) 1]);
        hp = patch(x',y',clut(ii,:),'EdgeColor','none');
        alpha(hp,.3)
    end
end

%%
axis off
scalebar('Unit','frames','Location','northeast')
