function make360movie(fname)
figure;
colormap gray;
writerObj = VideoWriter([strtok(fname,'.') '.mp4'],'MPEG-4');
writerObj.Quality = 96;
writerObj.FrameRate = 30;
open(writerObj);
for ii = 91 : 180
temp_tiff = imread(fname, ii);
imshow(temp_tiff);axis off
drawnow;
frame = getframe;
writeVideo(writerObj,frame);
end
close(writerObj);
