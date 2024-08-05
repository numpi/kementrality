function A = adjacency2(G, weights)
% replacement for Matlab's adjacency(G, weights), since it does not exist
% in 2017a.

if any(G.Edges.EndNodes(:,1) == G.Edges.EndNodes(:,2))
    error('Graph must not have self-loops')
end

A = sparse([G.Edges.EndNodes(:,1) G.Edges.EndNodes(:,2)], [G.Edges.EndNodes(:,2) G.Edges.EndNodes(:,1)], [weights weights], numnodes(G), numnodes(G));
