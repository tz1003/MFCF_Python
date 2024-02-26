function [values, sep_ranked] = greedy_sortsep(vtx, sep, W)
  sep_ranked = sep(~isnan(sep));
  pad = numel(sep) - numel(sep_ranked);
  [values, ranks] = sort(W(vtx, sep_ranked), 'descend');
  sep_ranked = [sep_ranked(ranks) NaN(1,pad)];
  values = [values zeros(1, pad)];
end