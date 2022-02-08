function m=sbxshownplane(fn)
%%
if nargin ==0
    [fname,p]=uigetfile('.align','load aligned data');
else
    [p,fnam,ext]=fileparts(fn);
    fname=[fnam ext];
end

load(fullfile(p,fname),'mr','mg','-mat');

if ~exist('mr')
    load(fullfile(p,fname),'m','-mat');
    mr = m;
end
if ~exist('mg')
    mg = m;
end

clear m

nplane=size(mg,3);
for i=1:nplane
    redchan = double(mr(:,:,i));
    redchanmax(i) = round(prctile(redchan(:),99));
    redchanmin(i) = round(prctile(redchan(:),10));

    greenchan = double(mg(:,:,i));
    greenchanmax(i) = round(prctile(greenchan(:),99));
    greenchanmin(i) = round(prctile(greenchan(:),10));
    m{i} =(redchan-redchanmin(i))/redchanmax(i);
    m{i}(:,:,2)=(greenchan-greenchanmin(i))/greenchanmax(i);
    m{i}(:,:,3)=0;
%     m{i}(:,:,2)=0;
end

sbxread(fullfile(p,strtok(fname,'.')),1,1);global info;
try
    z=info.config.knobby.pos.z
end

figure('Position',[400 200 1700 700],'Name',fname);hold on
for i=1:nplane
    subplot(1,nplane,i);
    title(sprintf('r=(%d,%d),g=(%d,%d)',redchanmax(i),redchanmin(i),greenchanmax(i),greenchanmin(i)))
    imshow(m{i});
%     figure  % for example case print
%     temp = circshift(m{i},[14 -45 0]); %file1
%     temp = m{i}; %file2
%     imshow(temp(43:491,128:640,:));
%     colormap gray
end