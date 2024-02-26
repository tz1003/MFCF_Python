function [ranked_values, ranked_sep] = greedy_sortsep_v( vertices, sets, W )
  sort_fun = @(x) greedy_sortsep(vertices(x), sets(x,:), W);
  [val_cell, rank_cell] = arrayfun(sort_fun,(1:numel(vertices)).', 'UniformOutput', false);
  ranked_values = cell2mat(val_cell);
  ranked_sep = cell2mat(rank_cell);
end