function G = kementrality_der(basename, weights, relative, parallel)
% compute the derivative variant of the kemeny-based centrality from a csv file
%
% G = kementrality_der(basename, weights, relative, parallel)
%
% Compute the Kemeny-based centrality described in 
% https://doi.org/10.1137/22M1486728
%
% This is a single unified script meant to simplify the process.
%
% Usage:
%
% >> G = kementrality_der("map")
%
% or:
%
% >> G = kementrality_der("map", [], false)
%
% (the second form does not scale the measure by edge weight)
%
% The command reads the network from a file called "map.csv", which must contain one 
% line for each edge with columns x1, y1, x2, y2, and writes the result on 
% a file called "map_kementrality.csv". 
%
% A point is considered an intersection only if (x1,y1) or (x2,y2) match 
% exactly those in another road segment. Take care with approximations, e.g.,
% one coordinate being 1.5 on a segment and 1.49999999 on another.
%
% If the network is disconnected, the centrality is computed only on the
% largest connected component; on disconnected roads the column contains 
% NaN.
%
% Additional columns in map.csv are copied to this new file, and a column
% called kementrality is added.
% Auxiliary files map_nodes.csv, map_edges.csv, and map_dual_edges.csv are
% produced.
%
% Alternatively, the first argument can be a Matlab graph object.
%
% Additional parameters:
%   * weights: weight to use for each edge; it can be:
%     - a vector with the edges ordered as in G.Edges, 
%     - 'inv' (default), which is (1/length), where the (Euclidean) length
%       is computed from G.Nodes.x and G.Nodes.y
%     - 'exp', i.e., exp(edge_length/max_length), as used in our first papers.
%   * relative: wheter to scale with edge weight; defaults to true.
%   * parallel: whether to use parallel computation. You may set it to false
%     if you don't have the Matlab parallel toolbox installed, or if you want
%     a faster result for a small network. Otherwise just use the default
%     true.
%   
% To set parallel to false without touching the other defaults, use
% >> G = kementrality("map", [], [], false);
%
% Return value: a Matlab graph object, which includes G.Edges.kementrality.
%
% To plot the map inside Matlab, you can use
%
% plotmeasure(G, kementrality_der);

if not(exist('parallel', 'var')) || isempty(parallel)
    parallel = true;
end
if not(exist('relative', 'var')) || isempty(relative)
    relative = true;
end

if isa(basename, "char")
    G = convert_graphs(basename, true);
else
    G = basename;
end

if exist('weights', 'var') && strcmp(weights, 'exp')
    G.Edges.Length = hypot(diff(G.Nodes.x(G.Edges.EndNodes)'), diff(G.Nodes.y(G.Edges.EndNodes)'))';
    weights = exp(-G.Edges.Length/max(G.Edges.Length));
end

if not(exist('weights', 'var')) || isempty(weights) || strcmp(weights, 'inv')
    G.Edges.Length = hypot(diff(G.Nodes.x(G.Edges.EndNodes)'), diff(G.Nodes.y(G.Edges.EndNodes)'))';
    weights = 1 ./ G.Edges.Length;
end


kementrality_der = kderivative(G, weights, parallel, relative);

G.Edges.kementrality_der = kementrality_der;

if isa(basename, "char")
    % creates a reduced Edges table that undoes the switching in
    % convert_graphs
    x1 = G.Edges.x1;
    x2 = G.Edges.x2;
    y1 = G.Edges.y1;
    y2 = G.Edges.y2;
    [x1(G.Edges.switched), x2(G.Edges.switched)] = deal(x2(G.Edges.switched), x1(G.Edges.switched));
    [y1(G.Edges.switched), y2(G.Edges.switched)] = deal(y2(G.Edges.switched), y1(G.Edges.switched));
    T2 = table(x1, y1, x2, y2, kementrality_der);

    T1 = readtable(strcat(basename, '.csv'));
    T1 = outerjoin(T1, T2, 'Type', 'left', 'MergeKeys', true);
    writetable(T1, strcat(basename, '_kementrality_der.csv'));
end

plotmeasure(G, kementrality_der);
