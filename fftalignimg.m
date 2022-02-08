function  [u v]=fftalignimg(align1img,align2img)

fig3 = figure();
while menu('Select areas for image alignment?','Yes','Done')==1
    imshowpair(align1img,align2img);
    rect=round(getrect(fig3));
    img1=align1img(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3));
    img2=align2img(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3));
    [u v] = fftalign(img1,img2);
    align1imgnew = circshift(align1img,[u,v]);
    imshowpair(align1imgnew,align2img);
    imgcorr = corr(double(align1imgnew(:)),double(align2img(:)));
    title(sprintf('u=%d,v=%d,corr=%.2f',u,v,imgcorr))
    %     set(fig3,'Position',[0 0 2000 1200]);
end