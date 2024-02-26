function [ cliques, separators, peo, tree ] = MFCF_Forest( X, ct_control, gain_function )
%MFCF Builds a clique forest from data
% See:
% Massara, G. P., & Aste, T. (2019). Learning Clique Forests. 
% arXiv preprint arXiv:1905.02266.


% Analyse X
sizeX = size(X);
isSquare = sizeX(1) == sizeX(2);
p = sizeX(2);
max_cliques = p - 1;

% Preallocate output
cliques = NaN(max_cliques, ct_control.max_clique_size);
separators = NaN(max_cliques - 1, ct_control.max_clique_size - 1);
peo = [];
tree = sparse(max_cliques, max_cliques);
GT = init_gain_table(ct_control.max_clique_size);

outstanding_nodes = 1:p;
clique_no = 0;
sep_no = 0;

% If matrix is not square make kernel matrix
if ~isSquare 
    kernel_func = @(r,c) gain_function(X, r, c, ct_control);
    ij = nchoosek(1:p,2);
    vals = arrayfun(kernel_func, ij(:,1), ij(:,2));
    K = zeros(p,p);
    idx = sub2ind(size(K), ij(:,1), ij(:,2));
    K(idx) = vals;
    K = K + K.';
    sums = sum(K.*(K>mean(K(:))),2);
else
    sums = sum(X.*(X>mean(X(:))),2);
end

% get first ct_control.min_clique nodes as the first clique
[~,j]=sort(sums,'descend');
first_clique = j(1:ct_control.min_clique_size).';
clique_no = 1;
cliques(clique_no, 1:numel(first_clique)) = first_clique;
% remove first clique form outstanding_nodes
outstanding_nodes = outstanding_nodes(~ismember(outstanding_nodes, first_clique));
peo = first_clique.';
% calculate gains for this clique
%gf = @(v) gain_function(X, first_clique, v, ct_control);
%[nodes gains seps] = arrayfun(gf, outstanding_nodes.', 'UniformOutput', false);
[nodes gains seps] = gain_function(X, first_clique, outstanding_nodes, ct_control);
new_gains = numel(gains);
from = GT.tot_records + 1;
to   = from + new_gains - 1;
GT.tot_records = GT.tot_records + new_gains;
idx = (from:to).';
clq_len = numel(first_clique);
GT.cliques(idx, 1:clq_len) = repmat(first_clique,new_gains,1);
if clq_len < ct_control.max_clique_size
    GT.cliques(idx, (clq_len+1):end) = NaN;
end
GT.gains(idx) = gains;
GT.separators(idx, :) = seps;
GT.nodes = nodes;


while(~isempty(outstanding_nodes))

   [the_gain, idx] = max(GT.gains, [], 'omitnan');

   % case no gain, we must add an isolated clique, let us pick up the first 
   % outstanding node
   if isnan(the_gain)
     the_node = outstanding_nodes(1);
     the_sep  = [];
     parent_clique = []; 
     parent_clique_id = NaN;
   else    
     idx = idx(1); % keep only first match
     the_node = GT.nodes(idx);
     the_sep  = GT.separators(idx,:);
     the_sep = the_sep(~isnan(the_sep));
     parent_clique = GT.cliques(idx,:);
     parent_clique_id = id_from_set(cliques, parent_clique);
   end;
   
   new_clique = [the_sep the_node];
   outstanding_nodes = outstanding_nodes(outstanding_nodes ~= the_node);
   peo = [peo; the_node];
   
   %togo = numel(outstanding_nodes)
   
   % track if it is a clique extension or a new clique
   clique_extension = 0;
   
   % in case of no gain add an isolated clique
   if isnan(the_gain)
     clique_no = clique_no +1;
     clique_to_update = clique_no;
     cliques(clique_to_update, 1:numel(new_clique)) = new_clique;
   % case new clique with existing intersection
   elseif numel(new_clique) <= sum(~isnan(parent_clique))
       clique_no = clique_no + 1;
       clique_to_update = clique_no;
       sep_no = sep_no + 1;
       cliques(clique_to_update, 1:numel(new_clique)) = new_clique;
       separators(sep_no, 1:numel(the_sep)) = the_sep;
       tree(clique_to_update, parent_clique_id) = 1;
   % case extension of existing clique
   else
       clique_to_update = parent_clique_id;
       cliques(clique_to_update, 1:numel(new_clique)) = new_clique;
       clique_extension = 1;
       % parent clique is the same as the clique without the last estension
       old_clique_idx = id_from_set(GT.cliques, parent_clique);
   end
   
   
   % don't update the gain table when nodes finished
   if isempty(outstanding_nodes) 
       break;
   end;
    
   % finally update gain table
   [nodes gains seps] = gain_function(X, new_clique, outstanding_nodes, ct_control);
   new_gains = numel(gains);
   from = GT.tot_records + 1;
   to   = from + new_gains - 1;
   GT.tot_records = GT.tot_records + new_gains;
   idx = (from:to).';
   new_clique = new_clique(~isnan(new_clique));
   clq_len = numel(new_clique);
   GT.cliques(idx, 1:clq_len) = repmat(new_clique,new_gains,1);
   if clq_len < ct_control.max_clique_size
      GT.cliques(idx, (clq_len+1):end) = NaN;
   end
   GT.gains(idx) = gains;
   GT.separators(idx, :) = seps;
   GT.nodes(idx)= nodes;
   
   % remove from gain table where node is the new one
   idx = find(GT.nodes == the_node);
   GT.gains(idx) = NaN;
   
   % if the clique was expanded remove the records with the old clique from the 
   % gain table
   if clique_extension == 1
     GT.gains(old_clique_idx) = NaN;
   end

   % remove separators who have reached max cooridination number
   if ct_control.coordination_num > 0
       if (~isempty(the_sep))
           % count appearance of separator just used
           idx = id_from_set(separators, the_sep);
           if numel(idx) >= ct_control.coordination_num
              idx = id_from_set(GT.separators, the_sep);
              GT.gains(idx) = NaN;
           end;
       end;
   end


   % if drop sep remove also the separator just used
   if ct_control.drop_sep == true
       idx = id_from_set(GT.separators, the_sep);
       GT.gains(idx) = NaN; 
   end
   

end
   peo = flipud(peo);

end
