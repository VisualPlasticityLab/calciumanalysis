function plotROI

yr=423:722
xr=46:346
xlim([xr(1) xr(end)])
ylim([yr(1) yr(end)])
ROI=fig1imgnew(xr,yr);
ROImin=min(ROI(:));
ROImax=max(ROI(:));
cdata1=(fig1imgnew-ROImin)/(ROImax-ROImin);
cdata1(:,:,2)=0;
cdata1(:,:,3)=0;
figure;imagesc('CData',cdata1)
axis tight;axis off
clc

ROI2=fig2img(xr,yr);
ROI2min=min(ROI(:));
ROI2max=max(ROI(:));
cdata2(:,:,2)=(fig2img-ROI2min)/(ROI2max-ROI2min);
cdata2(:,:,1)=0;
cdata2(:,:,3)=0;
figure;imagesc('CData',cdata2)
axis off
axis tight

cdata=cdata1;
cdata(:,:,2)=cdata2(:,:,2);
figure;hold on
imagesc('CData',cdata)
axis off
axis tight
for i=1+cellnum1-pair1
    j=i+cellnum1;
%     text(kC1(i).Position(1)+v,kC1(i).Position(2)+u,kC1(i).String,'Color',[ 1 0 0])
    contour(kC1(j).XData,kC1(j).YData,circshift(kC1(j).ZData,[u v]),[0.01 1],'LineColor',[1 0 0])
end
%2nd image, green
for i=1+cellnum2-pair2
    j=i+cellnum2;
%     text(kC2(i).Position(1)-10,kC2(i).Position(2)-10,kC2(i).String,'Color',[ 0 1 0])
    contour(kC2(j).XData,kC2(j).YData,kC2(j).ZData,[0.01 1],'LineColor',[0 1 0])
end

sigONEplt(eye2.SI,eye2.sigF,eye2.matrix,win_sig2,pair2(Cor2))
sigONEplt(eye1.SI,eye1.sigF,eye1.matrix,win_sig1,pair1(Cor2))