% function fns = findmatchedplane
% % NLW scope Feb 2016
% spatial calibration
%  4X .54micro/pixel  .316micro/pixel  
%  2X 1.073micro/pixel .632micro/pixel
%  1X 2.217micro/pixel(y) 1.36micro/pixel(x)
% sbxshownplane
ranges = 15:-5:-15
startingfolder =100;
%%
global scanbox_h sbconfig
h = guihandles(scanbox_h);
% set animal and unit numbers
h.expt.String = num2str(startingfolder);  
while exist([h.animal.String '_' h.unit.String '_' h.expt.String '.sbx'],'file') 
    startingfolder = startingfolder+100;
    h.expt.String = num2str(startingfolder);  
end
h.expt.Callback(h.expt,[]);
 h.volscan.Value = 1; h.volscan.Callback(h.volscan,[]);  % volumetric mode
h.tfilter.Value = 3; h.tfilter.Callback(h.tfilter,[]);  % accumulate mode
h.savesel.Value = 2; h.savesel.Callback(h.savesel,[]);% save green+red / green/ red channel
step = ranges(2)-ranges(1);
% defult take 80frames for each step
h.frames.String = num2str(120*str2double(h.optoperiod.String));h.frames.Callback(h.frames,[]); 
%  h.frames.String = num2str(180);h.frames.Callback(h.frames,[]); 

%% Grab
fns = [];
tri_send('KBY',0,12,0); % set knobby to super fine
tri_send('KBY',0,0,ranges(1));
for i = 1:numel(ranges)
    fns{i} = [h.animal.String '_' h.unit.String '_' h.expt.String ];
    h.grabb.Callback(h.grabb,[]);
    pause(1);
    tri_send('KBY',0,0,step);
end
tri_send('KBY',0,0,-ranges(end)-step);
%% Analyze
cd([h.dirname.String '\' h.animal.String]);
sbxalignnplanedir(1,fns{:});

%% Previous image
clear m mg
[fname,p]=uigetfile('.align','load previous aligned  images');
load(fullfile(p,fname),'m','mg','-mat');
greenimage = menu('Use which Channels?','green','red')
if greenimage ==1 & ~exist('m','var')
    m = mg ;
end
B = m ; %previous
z0 = squeeze(mean(m,3));
N = min(size(z0));    % leave out margin and focus on center
yidx = round(size(z0,1)/2)-N/2 + 1 : round(size(z0,1)/2)+ N/2;
xidx = round(size(z0,2)/2)-N/2 + 1 : round(size(z0,2)/2)+ N/2;
z0 = z0(yidx,xidx);
%% COMPARE the cells
clear z1 z2 A
xc = zeros(size(B,3),numel(ranges));
f = figure('Position',[100 200 2000 720]);

for i= 1:numel(ranges)
    load([fns{i} '.align'],'m','mg','-mat');
    if greenimage ==1& ~exist('m','var')
     m = mg ;
    end
    z1 = squeeze(mean(m,3));
    z1 = z1(yidx,xidx);
    [u(i),v(i)] = fftalign(z1,z0);
    
    A = circshift(m,[u(i) v(i) 0]);
    z2(:,:,:,i) = A;

    for k=1:size(B,3)
        old = double(B(yidx,xidx,k));
        new = double(A(yidx,xidx,k));
        [rho,~] = corrcoef(old,new); %(old>prctile(old(:),40))
        xc(k,i) = rho(1,2);
        subplot(1,size(A,3),k);
        imshowpair(B(yidx,xidx,k),A(yidx,xidx,k));
        title(sprintf('Old Red,New Green,Moved%d, Cor=%.2f',ranges(i),xc(k,i)));
        drawnow;
    end
    pause(5);
end

%% compute peak locations
[bestmatch, depth]= max(xc,[],2);
% print aligned and movement suggestion
f4=figure;
for k=1:size(A,3)
subplot(1,size(A,3),k);
imshowpair(B(yidx,xidx,k),z2(yidx,xidx,k,depth(k)));
title(sprintf('Old Red,New Green,BestMatch: Move%d,cor=%.2f',ranges(depth(k)),bestmatch(k)));
drawnow;
end
hf=figure;hold on;
plot(ranges,xc','o');
title('correlation for each plane');
legend plane1 plane2;
drawnow;
% set(gca,'XTickLabel',ranges)

load('magnificationlist.mat');
mag = maglist(h.magnification.Value);
x = median(u(depth))*(1.36/mag);
y = median(v(depth))*(-2.217/mag);
z = mean(ranges(depth));
%% Move to better location, take pictures, and come back
tri_send('KBY',0,2,x);
tri_send('KBY',0,1,y);
tri_send('KBY',0,0,z);

fn = [h.animal.String '_' h.unit.String '_' h.expt.String ];
h.grabb.Callback(h.grabb,[]);
sbxalignnplanedir(2 - greenimage ==1,fn);

load([fn '.align'],'m','mg','-mat');
if greenimage ==1& ~exist('m','var')
    m = mg ;
end
z1 = squeeze(mean(m,3));
z1 = z1(yidx,xidx);

[u0,v0] = fftalign(z1,z0);
% A = m;
 A = circshift(m,[u0 v0 0]);
% f = figure('Position',[400 800 1200 400]);
figure(f4);
for k=1:size(A,3)
    old = double(B(:,:,k));
    new = double(A(:,:,k));
    [rho,p] = corrcoef(old,new);
    subplot(1,size(A,3),k);
    imshowpair(B(:,:,k),A(:,:,k));
    title(sprintf('Old Red,New Green,Corrected Plane:Cor=%.2f',rho(1,2)));
    drawnow;
end
saveas(f4,[fn 'Correlatewith' strtok(fname,'.') '.fig'])

tri_send('KBY',0,2,-x);
tri_send('KBY',0,1,-y);
tri_send('KBY',0,0,-z);
x = x+ u0*(1.36/mag);
y = y+ v0*(-2.217/mag);
disp([x y z])
% fclose(info.fid);
% info = [];
figure(hf);
axis on