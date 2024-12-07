function [GR, mapback, weights] = reduce_graph(G, weights)
% Create a reduced graph GR that removes all degree-2 edges.
% 
% If the original (undirected) graph contained a path i -> j1 -> ... -> jl -> k, 
% where j1,...jl are not connected to any other node, we replace j1,...,jl 
% with a single edge from i to k.
% This procedure may create multiple edges between the same two nodes i and
% k. If this happens, we replace them with a single edge.
% Weights for the reduced graph weightsR are generated from weights as if
% weights were conductances, using parallel/series resistance rules:
% * if a path is replaced by a single edge, its weight satisfies
% 1/w = 1/w1 + 1/w2 + ... + 1/wl.
% * if multiple edges between the same nodes i,k are replaced by a single 
% edge, its weight is the sum w1 + ... + wm.
%
% Also returns a vector "mapback" such that mapback(l) = k whenever 
% edge l in G is a section of edge k in GR. In this way, once computed a
% measure on the graph kementralityR, you can obtain the corresponding
% measure on G via kementrality = kementralityR(mapback)

if isa(G, "char")
    % supports a basename as first input
    G = convert_graphs(G, true);
end
if not(exist('weights', 'var')) || isempty(weights)
    G.Edges.Length = hypot(G.Edges.x1-G.Edges.x2, G.Edges.y1-G.Edges.y2);
    weights = exp(-G.Edges.Length/max(G.Edges.Length));
end

% We call a "bucket" a sequence of edges to reduce, i.e., 
%  i -> j1 -> ... -> jl -> k, where deg(j1)==...==deg(jl)==2.
% A bucket contains two "extremal edges" (i,j1) and (jl,k)
% and zero or more "internal edges" (j1,j2), ... (j{l-1},jl)

% We construct a subgraph G2 that contains only the j-edges. 
% In this way, each bucket is a connected component of G2.
% We number them using conncomp()

jnodes_logical = degree(G)==2;
G2 = subgraph(G, jnodes_logical);
% find the bucket and weight of each edge in G2
buckets2_nodes = transpose(conncomp(G2));
buckets2_edges = buckets2_nodes(G2.Edges.EndNodes(:,1));

internal_edges_logical = sum(jnodes_logical(G.Edges.EndNodes),2)==2;
weights2 = weights(internal_edges_logical);

% vector that tells which bucket a node in G belongs to, or 0
% if it does not belong to G2
nodeG2comp = zeros(size(G.Nodes,1), 1);
nodeG2comp(jnodes_logical) = buckets2_nodes;
internal_buckets = nodeG2comp(G.Edges.EndNodes);

% find extremal edges in G
extremal_edges = sum(jnodes_logical(G.Edges.EndNodes),2) == 1;
% and their corresponding bucket
% since each row in internal_buckets(boundary_edges,:) contains exactly one
% zero, we can just sum the values.
extremal_buckets = sum(internal_buckets(extremal_edges,:),2);

% compute the conductance of each bucket
resistance_extremal_edges = accumarray(extremal_buckets, 1./weights(extremal_edges));
resistance_internal_edges = accumarray(buckets2_edges, 1./weights2, [length(resistance_extremal_edges) 1]);

weightsBuckets = 1./(resistance_extremal_edges + resistance_internal_edges);

% identify nodes i,k corresponding to each bucket.

% endpoint_node contains for each extremal_edge the index of its endpoint
% node i or k.
idx = transpose(1:numnodes(G)); idx(jnodes_logical) = 0;
endpoint_node = sum(idx(G.Edges.EndNodes(extremal_edges,:)),2);

% we sort endpoint_node so that its first two entries correspond to the
% endpoint nodes of bucket 1, the next ones of bucket 2, and so on.
[~, perm] = sort(extremal_buckets); % the first output must be [1 1 2 2 3 3 ...]
endpoint_node_sorted = endpoint_node(perm);
newEndNodes = transpose(reshape(endpoint_node_sorted, 2, []));

% create a mapback vector. It will have to be reordered, following the
% ordering of edges in GR, which is still to determine
% for now we identify edges in GR with a GR_index that is defined as
% the edge index in G for edges that were in G, and
% numedges(G)+bucket_number for edges that come from the reduction.
mapbacktemp = transpose(1:numedges(G));
mapbacktemp(extremal_edges) = numedges(G) + extremal_buckets;
mapbacktemp(internal_edges_logical) = numedges(G) + internal_buckets(internal_edges_logical);

% we are ready to change G
GR = G;
% columns in G.Edges other than EndNodes carry little meaning
% since we will not be able to define them for the new edges, so we
% remove them.
GR.Edges(:,2:end) = [];
GR.Edges.weight = weights;
GR.Edges.original_index = transpose(1:size(GR.Edges,1));
E = table();
E.EndNodes = newEndNodes;
E.weight = weightsBuckets;
E.original_index = transpose(size(GR.Edges,1)+1:size(GR.Edges,1)+size(newEndNodes,1));
GR = addedge(GR, E);
GR = rmnode(GR, find(jnodes_logical));

% mapback contains original_indices, we must convert them into positions in
% GR.Edges
[~, mapback] = ismember(mapbacktemp, GR.Edges.original_index);
GR.Edges.original_index = [];

[GR,map] = simplify(GR, 'sum', 'keepselfloops', 'AggregationVariables', 'weight');
mapback = map(mapback);

weights = GR.Edges.weight;
GR.Edges.weight = [];
