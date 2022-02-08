function sbxballmotiondir(s)

if(nargin>=1)
    d=dir([s '_ball.mat']);
else
    d = dir('*_ball.mat');
end

for i=1:length(d)
    try
        fn = d(i).name;
        vars = whos('-file',fn);
         if ~ismember('ball', {vars.name})
             fprintf('extract ball motion for %s\n',fn)
             sbxballmotion(strtok(fn,'.'));
         else
            load(fn,'ball','-mat');
            if prctile(abs(ball),5)>106
             fprintf('reextract ball motion for %s\n',fn)
             sbxballmotion(strtok(fn,'.'));
            else
                sprintf('Ball motion %s is already extracted',fn)
            end
        end
    catch
        fprintf('Could not do ballmotion for %s\n',fn)
    end
end

allfiles=dir;
for i=3:size(allfiles,1)
    if isdir(allfiles(i).name)
        disp(pwd);
        cd(allfiles(i).name);sbxballmotiondir; cd ..
    end
end

