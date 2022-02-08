function sbx2tifdir
% nbinning =1;

d=dir('*.sbx');
% temp = pwd;
% pos = strfind(temp,filesep);
% foldern = [temp(pos(end-1)+1 :end)];
% mkdir(['/nas/data1/mdadarla/CA2p/J1/' foldern])
if numel(d)>0
    for i=1:numel(d)
        fn = strtok(d(i).name,'.');
        if exist('*_*_*.tif','file')
            delete('*_*_*.tif');
        end
        sbx2tif_green_and_red(fn);
    %     mkdir(['/nas/data1/mdadarla/CA2p/' foldern filesep num2str(i)]);
    %     movefile([fn '*.tif'],['/nas/data1/mdadarla/CA2p/J1/' foldern '/' num2str(i)]);
    end
end

allfiles=dir;
for i=3:size(allfiles,1)
    if isdir(allfiles(i).name)
        disp(pwd);
        cd(allfiles(i).name);sbx2grtifdir; cd ..
    end
end

clear