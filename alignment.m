function [x,y]=alignment(day0,day1,mag)
%default widefield , if mag exist, then it's 2p mode
if nargin<3
%     micronxy=[.55 .32];
    micronxy=[3.4 -2];
else
    micronxy=[1.36 -2.217]/mag;%  1X 2.217micro/pixel(y) 1.36micro/pixel(x)

end

if nargin<2
    day0f = uigetfile('*.mat','select snapshot before stroke');
    day1f = uigetfile('*.mat','select snapshot after stroke');
    load(day0f,'img'); day0=img; clear img;
    load(day1f,'img'); day1=img; clear img;
end
    
fixed=double(day0);%fixed=img(:,:,2);
moving=double(day1);%moving=img(:,:,2);
moving=moving*(max(fixed(:))/max(moving(:)));
rho = corrcoef(moving(:),fixed(:));


figure;
fig2=subplot(1,2,1);
imshowpair(moving,fixed); %first image day1 magenta, second green
title(sprintf('day1 red,day0 green: Correlation=%.2f ',rho(1,2)));
drawnow;
[optimizer, metric] = imregconfig('multimodal');
[moving_reg,R_reg] =  imregister(moving, fixed, 'translation', optimizer, metric);
tform = imregtform(moving, fixed, 'translation', optimizer, metric);
disp(tform.T);
%have to set 'OutputView' mode
movingRegistered= imwarp(moving,tform,'OutputView',imref2d(size(moving)));
rho2 = corrcoef(movingRegistered(:),fixed(:));

subplot(1,2,2);
imshowpair(movingRegistered,fixed)
title(sprintf('Corrected day1 red,day0 green:Correlation=%.2f',rho2(1,2)));
%             v(i)=round(tform.T(3,1));
%             u(i)=round(tform.T(3,2));
    x = tform.T(3,2)*micronxy(1);%3.4;
    y = tform.T(3,1)*micronxy(2);%-2;

disp(sprintf('Method1:move[x y]=[%.1f %.1f]',x,y))


disp('select area on the map for method 2')
fig1img=moving;
fig2img=fixed;
 rect=round(getrect(fig2));
    img1=fig1img(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3));
    img2=fig2img(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3));
[u v] = fftalign(img1,img2);
fig1imgnew = circshift(fig1img,[u,v]);
figure;imshowpair(fig1imgnew,fig2img)
    x2 = u*micronxy(1);%3.4;
    y2 = v*micronxy(2);%-2;
disp(sprintf('Method2:move[x y]=[%.1f %.1f]',x2,y2))
