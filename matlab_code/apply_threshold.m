function [val, sep] = apply_threshold(val, sep, mincsize, threshold)
  idx = find(val < threshold);
  idx = idx(idx > mincsize - 1);
  val(idx) = 0;
  sep(idx) = NaN;
 end