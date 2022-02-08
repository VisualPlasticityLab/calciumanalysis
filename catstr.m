function S1 = catstr(S1,S2)

 f = fieldnames(S2);
 for i = 1:numel(f)
     n = ndims(S1.(f{i}));
    S1.(f{i}) = cat(n,S1.(f{i}),S2.(f{i}))
 end