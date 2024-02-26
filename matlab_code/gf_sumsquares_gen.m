% Copyright (C) 2021 XM78AT
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

% -*- texinfo -*-
% @deftypefn {} {@var{retval} =} gf_sumsquares_gen (@var{input1}, @var{input2})
%
% @seealso{}
% @end deftypefn

% Author: XM78AT <XM78AT@WPU8L0063505>
% Created: 2021-03-29

function [nodes gains seps] = gf_sumsquares_gen (M, clq, v, ct_control)
%GF_SUMSQUARES_GEN returns the sum of squares subjet to max_clique_size, 
% min_clique_size, threshold constraints


% remove NaN and get sizes
clq = clq(~isnan(clq));
csz = numel(clq);
vn = numel(v);

% generate weight matrix to be used for sorting the elements of the separtors
W = M.^2;

% extract the threshold ...
if  isfield(ct_control, 'threshold')
    threshold = ct_control.threshold;
else
    threshold = 0.0;
end;

% ... and size of the cache. By default we keep as many as the TMFG
if  isfield(ct_control, 'cachesize')
    cachesize = ct_control.cachesize;
else
    cachesize = min(4, ct_control.min_clique_size);
end;

% ... and min size of cliques. By default we keep as much as the TMFG
if  isfield(ct_control, 'min_clique_size')
    mincsize = ct_control.min_clique_size;
else
    mincsize = 4;
end;

% ... and max size of cliques. By default we keep as much as the TMFG
if  isfield(ct_control, 'max_clique_size')
    maxcsize = ct_control.max_clique_size;
else
    maxcsize = 4;
end;

% if the clique is less than max_clique_size then we keep the clique itself ...
if csz < maxcsize 
  facets = clq;
 else
% otherwise we generate the facets of the clique
  facets = nchoosek(clq, csz-1);
end

% now let's generate the table for the calculation of gains
% get number of facets
[block_rows, ncols] = size(facets);
% repeat the outstanding vertices times the number of facets
the_vs = reshape(repmat(v, block_rows, 1), block_rows*vn,1);
% repeat facets times the number of outstanding vertices
the_fs = repmat(facets, vn, 1);

% now rank the facets (separators) greedily
[ranked_values, ranked_seps] = greedy_sortsep_v(the_vs, the_fs, W);

% apply thresholding to gains of elements between min_clique_size and max_clique_size
[ranked_values_thr, ranked_seps_thr] = apply_threshold_v(ranked_values, ranked_seps, mincsize, threshold);


% apply gain function
gainfun = @(rowidx) sum(ranked_values_thr(rowidx,:));
% see which gains are above threshold
gains = arrayfun(gainfun,1:numel(the_vs));

selector = repmat(transpose(1:block_rows),vn,1);
the_table = [the_vs, gains', ranked_seps_thr];
the_table = sortrows(the_table, [1,-2] );
%compress the table by keeping only -cachesize- entries per node;
the_table = the_table(selector <= cachesize, :);

nodes = the_table(:,1);
gains = the_table(:,2);
seps = [the_table(:,3:end), NaN(numel(nodes),maxcsize - 1- ncols)];


end


