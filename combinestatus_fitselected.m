function fit_selected = combinestatus_fitselected(fit_selected,day1matrix,day2matrix)

runtrials1 = sum(day1matrix)/size(day1matrix,1);
stilltrials1 = 1 - runtrials1;

runtrials2 = sum(day2matrix)/size(day2matrix,1);
stilltrials2 = 1 - runtrials2;

for i=1:numel(fit_selected)
    fit_selected{i}.peak1 = nanmean(cat(1,stilltrials1.*fit_selected{i}.peak1S,...
    runtrials1.*fit_selected{i}.peak1R))./nanmean(cat(1,stilltrials1.*~isnan(fit_selected{i}.peak1S),...
    runtrials1.*~isnan(fit_selected{i}.peak1R)));
    
    [fit_selected{i}.a1,fit_selected{i}.ori1] = nanmax(fit_selected{i}.peak1);

    fit_selected{i}.peak2 = nanmean(cat(1,stilltrials2.*fit_selected{i}.peak2S,...
    runtrials2.*fit_selected{i}.peak2R))./nanmean(cat(1,stilltrials2.*~isnan(fit_selected{i}.peak2S),...
    runtrials2.*~isnan(fit_selected{i}.peak2R)));

    [fit_selected{i}.a2,fit_selected{i}.ori2] = nanmax(fit_selected{i}.peak2);
end
