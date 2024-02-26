function [ JS ] = j_LoGo( S, cliques, separators)
%LOGO assembles the direct or inverse (recommended) matrix from cliques and
%separators
N = size(S,1);
JS = sparse(N,N);

for c=cliques.'
    clique = c(~isnan(c));
    JS(clique, clique) = JS(clique, clique) + inv(S(clique, clique));

end

for s=separators.'
    separator = s(~isnan(s));
    S(separator, separator) = JS(separator, separator) - inv(S(separator, separator));

end

end



