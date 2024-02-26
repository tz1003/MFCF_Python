function [ranked_values, ranked_seps] = apply_threshold_v (ranked_values, ranked_seps, mincsize, threshold)
[vsz, ~] = size(ranked_values);
threshold_fn = @(rowidx) apply_threshold(ranked_values(rowidx,:), ranked_seps(rowidx,:), ...
mincsize, threshold);
[val_cell, rank_cell] = arrayfun(threshold_fn,(1:vsz).', 'UniformOutput', false);
ranked_values = cell2mat(val_cell);
ranked_seps = cell2mat(rank_cell);
end