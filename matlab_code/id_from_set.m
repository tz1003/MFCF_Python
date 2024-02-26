function idx = id_from_set( set, an_element )
%ID_FROM_CLIQUE Find row ids of gain table from clique
partial = sum(ismember(set, an_element), 2);
idx = find(partial == sum(~isnan(an_element)));
end

