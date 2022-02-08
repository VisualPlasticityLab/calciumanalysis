function [sig, dsig]= func_preprocess_data(sig_chunk,startfram,endfram)

%%
debug =1;
numofcondition=numel(startfram);
if size(sig_chunk,1)==1
    sig_chunk=sig_chunk';
end
sig= sig_chunk;
dsig =sig_chunk;

%%
if debug ==1
    figure('Name','whole trace');    hold on;
    plot([startfram;startfram],[0;1]*ones(1,numofcondition),'k-')
    plot([endfram;endfram],[0;1]*ones(size(1,numofcondition)),'k--')
    beginning= 1:startfram(1)-1;
    plot(beginning,sig_chunk(beginning),'b.')
    xlim([0 beginning(end)]);
end

[decay_each,decay]=fitexp1(sig_chunk,100,0);

seg=min(endfram-startfram);
startfram(numofcondition+1)=min(startfram(end)+min(diff(startfram)),size(sig,2));
for i=1:numofcondition;
    if i==1
        prestim = 1:startfram(i)-1;
    else
        prestim = endfram(i-1):startfram(i)-1;
    end
    current = startfram(i):startfram(i+1)-1;
    future = startfram(i+1)-1:size(sig,2);
    if debug ==1
        plot(current,sig(current),'r--')
    end
    [~, x1] = findpeaks(sig(prestim));
    if isempty(x1)
        baseline =prestim;
    else
        baseline =prestim(x1(end):end);
    end
    try
        f = fit(baseline',sig(baseline),'exp1');
        if abs(f.b-decay)<std(decay_each)
            sig([current,future]) = sig([current,future])- f([current,future])';
        end
        sig([current,future])= sig([current,future])-min(sig(current));
        
    end
%     temp = detrend(sig([prestim,current]));
%     dsig (current) = temp (numel(prestim)+1:end);
%     try
%         dsig (current) = dsig (current) - dsig(current(1));
%     catch
%         disp(i)
%         disp(size(dsig(current)))
%     end
    if debug ==1
        plot(current,sig_chunk(current),'b.')
        plot(current,sig(current),'r')
%         plot(current,dsig (current),'k-')
        xlim([0 current(end)]);
        drawnow;
        waitforbuttonpress;
    end
end

% beginning= startfram(1):startfram(2)-1;
% sig(beginning)=sig(beginning)-min(sig(beginning));
% dsig(beginning) =sig(beginning);
% plot(beginning,sig_chunk(beginning),'b.')
% plot(beginning,sig(beginning),'r')

