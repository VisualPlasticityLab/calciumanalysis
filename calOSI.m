function [OSI,gOSI,DSI,gDSI,pref_ori1,pref_dir]=calOSI(peak)
%peak size: ori*ncell
peak = squeeze(peak);
if ndims(peak)~=2
 warning('reformat peak data...');
 for ncond=1:size(peak,1)
    [OSI(ncond,:),gOSI(ncond,:),DSI(ncond,:),gDSI(ncond,:),pref_ori1(ncond,:),h]=calOSI(squeeze(peak(ncond,:,:)));
 end
 return;
end
ori=size(peak,1);
ncell=size(peak,2);
% badcells=squeeze(sum(sum(peak,1),2)==0);

%if there is blank get rid of it
ori1=ori-mod(ori,2);
peak = peak(1:ori1,:);
  
% tuning curve fitting
% pref0=nanmax(peak_p,[],1);
% for ii=1:ncell
%     peak_p2 = fit2peakgaussion(pref0,peak(:,ii))
% end
peak_p2=peak;

[pref1,pref_ori1]=nanmax(peak_p2,[],1);
pref_ori2=mod(pref_ori1+ori1/2-1,ori1)+1;
pref_dir=mod(pref_ori1-1,ori1/2)+1;
orth_ori1=mod(pref_ori1+ori1/4-1,ori1)+1;
orth_ori2=mod(pref_ori1-ori1/4-1,ori1)+1;

pref2 = peak_p2(sub2ind([ori1,ncell],pref_ori2,1:ncell));
orth1 = peak_p2(sub2ind([ori1,ncell],orth_ori1,1:ncell));
orth2 = peak_p2(sub2ind([ori1,ncell],orth_ori2,1:ncell));

orth = (orth1+orth2)/2;
OSI = (pref1-orth)./(pref1+orth);
DSI = (pref1-pref2)./(pref1+pref2);

gOSI = abs(exp(2i*(0:180/ori1:180-180/ori1)./180*pi)*peak_p2)./sum(peak_p2);
gDSI = abs(exp(i*(0:180/ori1:180-180/ori1)./180*pi)*peak_p2)./sum(peak_p2);

% OSI(:,:,badcells)=NaN;
% pref_ori1(badcells) =NaN;
% 
% if type ==2 %left and right eye
% h=figure;
% subplot(2,2,1)
% scatter(pref1(:,:,1),pref1(:,:,2))
% xlabel('contra')
% ylabel('ipsi')
% title('Maxium response on contra/ipsi eye')
% axis square
% subplot(2,2,2)
% scatterjitter(pref_ori1(:,:,1),pref_ori1(:,:,2))
% xlabel('contra')
% ylabel('ipsi')
% title('Preferred ori on contra/ipsi eye')
% axis square
% subplot(2,2,3)
% try
% histplt2(pref1(:,:,1),pref1(:,:,2))
% legend('Contra eye','Ipsi eye')
% end
% subplot(2,2,4)
% try
% histplt2(pref_ori1(:,:,1),pref_ori1(:,:,2))
% legend('Contra eye','Ipsi eye')
% end
% else
%     h=[];
% end
