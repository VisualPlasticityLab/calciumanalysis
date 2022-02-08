
function fname=memmapdir(type)
%type, predefine memmapchoice type: single file,multiple plane,Multiple
%files,mutiple planesmultiplefiles
%string, specify filenames to be aligned

setpref('Internet','SMTP_Server','keck.ucsf.edu')
setpref('Internet','E_mail','jsun@phy.ucsf.edu')

if ~nargin   %
    memmapchoice=menu('Which type?','Single file','Multiple planes','Multiple files','Multiple planes multiple files')
else
    memmapchoice=type;
end
d = dir('*.sbx');

sbxalignnplanedir;
% try
switch memmapchoice
    case 2 %Single file
        for i=1:numel(d)
            fns{1} = strtok(d(i).name,'.');
            newnames= makememmapnfilesnplane(fns); %makememmap2plane(fn);
            for j=1:numel(newnames)
                pos = strfind(newnames{j},'_memmap.mat');
                fname= memmap(newnames{j}(1:pos-1));
            end
        end
        %     case 3 %Multiple files
        %         newname=makememmapnfiles;
        %         pos = strfind(newname,'_memmap.mat');
        %         fname=memmap(newname(1:pos-1));
        %         caanalysisdir
        %         disp(sprintf('finished %s in %d mins',fname,round(toc/60)));
        
    case 4 %Multiple files
        for i=1:numel(d)
            fns{i} = strtok(d(i).name,'.');
        end
        newnames=makememmapnfilesnplane(fns);
        for i=1:numel(newnames)
            pos = strfind(newnames{i},'_memmap.mat');
            fname=memmap(newnames{i}(1:pos-1));
        end
end
close all;
clear;
sbxballmotiondir;
sbxeyemotiondir;
try
    caanalysisdir;
end
close all;
setpref('Internet','SMTP_Server','keck.ucsf.edu')
setpref('Internet','E_mail','jsun@phy.ucsf.edu')
% sendmail('mcdadarlat@gmail.com','Memmap done,continue with caanalysis', ...
%          'Memmap finished');
%


%     sendmail({'j.suninchina@gmail.com','anuta278@yandex.ru'},'Finished Memmapdir', ...
%         [fns 'projects finished']);

% catch
%     sendmail({'j.suninchina@gmail.com','anuta278@yandex.ru'},'Error Memmap', ...
%         [fns 'got error']);
% end


% global info_loaded info
% if(~isempty(info_loaded))   % try closing previous...
%     try
%         fclose(info.fid);
%     catch
%     end
% end
% switch memmapchoice
%     case 1  %single file
%         if exist('s','var')
%             d = dir([s '.sbx']);
%         else
%             d = dir('*.sbx');
%         end
%         for i=1:numel(d)
%             fn = strtok(d(i).name,'.');
%             newname = makememmap1file(fn);
%             pos = strfind(newname,'_memmap.mat');
%             fname=memmap(newname(1:pos-1));
%             disp(sprintf('finished cell extraction for %s in %d mins',fname,round(toc/60)));
%         end
%     case 2 %Multiple plane
%         if exist('s','var')
%             d = dir([s '.sbx']);
%         else
%             d = dir('*.sbx');
%         end
%         for i=1:numel(d)
%             fns = {strtok(d(i).name,'.')};
%              newnames= makememmapnfilesnplane(fns); %makememmap2plane(fn);
%             for j=1:numel(newnames)
%                 pos = strfind(newnames{j},'_memmap.mat');
%                 fname= memmap(newnames{j}(1:pos-1));
%                 disp(sprintf('finished cell extraction for %s in %d mins',fname,round(toc/60)));
%             end
%         end
%     case 3 %Multiple files
%         newname=makememmapnfiles;
%         pos = strfind(newname,'_memmap.mat');
%         fname=memmap(newname(1:pos-1));
%         caanalysisdir
%         disp(sprintf('finished %s in %d mins',fname,round(toc/60)));
%     case 4 %Multiple files with multiple planes
%         if exist('s','var')
%             d = dir([s '.sbx']);
%         else
% %             d = uigetfile('*.sbx','Select the files to combine','MultiSelect', 'on');
%             d = dir('*.sbx');
%         end
%         for i=1:numel(d)
%             fns{i} = strtok(d(i).name,'.');
%         end
%         newnames=makememmapnfilesnplane(fns);
%         for i=1:numel(newnames)
%             pos = strfind(newnames{i},'_memmap.mat');
%             fname=memmap(newnames{i}(1:pos-1));
% %             for k=1:numel(fname)
% %              caanalysis(fname{k});
% %             end
%             disp(sprintf('finished cell extraction for %s in %d mins',cell2str(fname),round(toc/60)));
%         end
% end
%
%
% setpref('Internet','SMTP_Server','keck.ucsf.edu')
% setpref('Internet','E_mail','jsun@phy.ucsf.edu')
% sendmail('j.suninchina@gmail.com','Memmap done,continue with caanalysis', ...
%          'Memmap finished');
%
% sbxballmotiondir;
% sbxeyemotiondir;
% caanalysisdir;
%
% sendmail('anuta278@yandex.ru','Finished Memmapdir', ...
%          'All projects finished');
%
% sendmail('j.suninchina@gmail.com','Finished Memmapdir', ...
%          'All projects finished');
%