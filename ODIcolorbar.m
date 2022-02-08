function ODIcolorbar
odi_range = (-1:.01:1);
peak_range = (0:.05:.5)';
% peak_range(:) =1;
odi_bar= repmat(odi_range,numel(peak_range),1);
peak_column = repmat(peak_range,1,numel(odi_range));
ODIcolorbar = ODI2rgb(odi_bar,peak_column);

imshow(ODIcolorbar,'InitialMagnification',400)
axis off
% h=gca;
% h.YLabel.String = 'Peak(df/F)';
% % h.YTick=[ .5 11 21.5];
% % h.YTickLabel= {'0','0.5','1'}
% 
% h.XLabel.String = 'ODI';
% h.XTick=[ .5 22 41.5];
% h.XTickLabel= {'-1','0','1'}

