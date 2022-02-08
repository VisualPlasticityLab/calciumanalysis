function caanalysisdir(fast)

if nargin<1
    fast =.5; % means using default trace from suite2p
end

files=dir('*cell*.signals');
% files=dir('*cell_tau6_*.signals')
for i=1:numel(files)
    fn = files(i).name;
    if isempty(strfind(fn,'&'))
        try
            caanalysis([strtok(fn,'.') '.signals'],fast)
        catch
            disp(['cant process' strtok(fn,'.') '.signals'])
        end
    end
    
end
% 
% allfiles=dir;
% for i=3:size(allfiles,1)
%     if isdir(allfiles(i).name)
%         disp(pwd);
%         cd(allfiles(i).name);caanalysisdir; cd ..
%     end
% end
