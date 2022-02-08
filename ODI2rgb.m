function  RGB = ODI2rgb(ODI,peak)
if isnan(ODI)
    RGB = [1 1 1];
else
    h=(ODI+1)/2*.75;  % make ODI[-1:1] be a hue that goes from 0 to 0.9
    s = log(peak/.05); % make peak [0.05 .5] be a value that goes from 0 to 1
    s(s<0) = 0 ;
    s(s>1) = 1;
    v = sqrt(peak/.05)/sqrt(10);% make peak [0.05 .5] be a saturation that goes from 0 to 1
    v(v>1) = 1;
    RGB = squeeze(hsv2rgb(h,s,v));
end

        
