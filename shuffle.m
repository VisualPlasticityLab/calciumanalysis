 function v=shuffle(v)
     v=v(randperm(length(v)));
     %shuffle = @(v)v(randperm(numel(v)));

 end