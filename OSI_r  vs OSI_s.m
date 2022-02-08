figure;
for i=1:46 
    plot(1:6,peakR(1,:,i),'r*-',1:6,peakS(1,:,i),'b*-'); 
    title([num2str(i) 'gOSI_r:gOSI_s ' num2str(gOSI_R(:,:,i)) ':'  num2str(gOSI_S(:,:,i))])
    pause
end