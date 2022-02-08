function peak_selectedday = combinestatus(peak_selected,day1matrix,day2matrix)

runtrials1 = sum(day1matrix)/size(day1matrix,1);
stilltrials1 = 1 - runtrials1;

peak_selectedday1 = nanmean(cat(3,diag(stilltrials1)*squeeze(peak_selected(1,:,:)),...
    diag(runtrials1)*squeeze(peak_selected(2,:,:))),3);


runtrials2 = sum(day2matrix)/size(day2matrix,1);
stilltrials2 = 1 - runtrials2;

peak_selectedday2 = nanmean(cat(3,diag(stilltrials2)*squeeze(peak_selected(3,:,:)),...
    diag(runtrials2)*squeeze(peak_selected(4,:,:))),3);


peak_selectedday = permute(cat(3,squeeze(peak_selectedday1),squeeze(peak_selectedday2)),[3 1 2]);