function sbxeyemotiondir(s)

if(nargin>=1)
    d=dir([s '_eye.mat']);
else
d = dir('*_eye.mat');
end
for(i=1:length(d))
    try
        fn = d(i).name;
        vars = whos('-file',fn);

        if ((exist([fn(1:strfind(fn,'_eye')-1) '.eye'])==0) && ~ismember('eye', {vars.name}))
           display(sprintf('extract eye motion for %s',fn))
           sbxeyemotion(fn);
        else
           sprintf('Eye motion %s is already extracted',fn)
        end
    catch
        display(sprintf('Could not do eyemotion for %s',fn))
    end
end
