function [nodes gains sepas] = gf_sum_squares( M, clq, v, ct_control )
%GF_SUM_SQUARES2 Returns the sum of the squares of the elements of the
%matrix that are above a given threshold

% Guido Previde Massara 20-04-2020

% remove NaN
clq = clq(~isnan(clq));
sz = numel(clq);
vn = numel(v);

% extract threshold
if  isfield(ct_control, 'threshold')
    threshold = ct_control.threshold;
else
    threshold = 0.0;
end;

if  isfield(ct_control, 'cachesize')
    cachesize = ct_control.cachesize;
else
    cachesize = min(4, ct_control.min_clique_size);
end;

% generate facets
facets = nchoosek(clq, sz-1);
[block_rows, ~] = size(facets);
the_vs = reshape(repmat(v, block_rows, 1), block_rows*vn,1);
the_fs = repmat(facets, vn, 1);
gainfun = @(rowidx) sum(M(the_vs(rowidx),the_fs(rowidx,:)).^2);
% see which gains are above threshold
gains = arrayfun(gainfun,1:numel(the_vs));
selector = repmat(transpose(1:block_rows),vn,1);
the_table = [the_vs, gains', the_fs];
the_table = sortrows(the_table, [1,-2] );
%compress the table by keeping only -cachesize- entries per node;
the_table = the_table(selector <= cachesize, :);

nodes = the_table(:,1);
gains = the_table(:,2);
sepas = the_table(:,3:end);

end

