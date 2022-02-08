function [names,names_dir] = num2ori(nSteps)

nSteps1 = nSteps - mod(nSteps,2);

for i=1:nSteps1/2
    if mod((i-1)*360,nSteps1)==0
        names{i} = sprintf('%d/%d',(i-1)*360/nSteps1,(i-1)*360/nSteps1+180);
    else
        names{i} = sprintf('%.1f/%.1f',(i-1)*360/nSteps1,(i-1)*360/nSteps1+180);
    end
end
if nSteps>nSteps1
    names{end+1} = 'blk';
end

for i=1:nSteps1
    if mod((i-1)*360,nSteps1)~=0
        names_dir{i} = sprintf('%.1f',(i-1)*360/nSteps1);
    else
        names_dir{i} = sprintf('%d',(i-1)*360/nSteps1);
    end
end
if nSteps>nSteps1
    names_dir{nSteps} = 'blk';
end