function G = kementrality_old_negative(basename, weights, parallel)
% compute kemeny-based centrality from a csv file
%
% G = kementrality(basename, weights, parallel)
%
% Compute the old version of the Kemeny-based centrality used in 
% https://www.hvl.no/globalassets/hvl-internett/arrangement/2022/13sss/522altafini.pdf. 
% Note that this measure can take negative values.
%
% This is a single unified script meant to simplify the process.
%
% Usage:
%
% >> G = kementrality("map")
%
% reads the network from a file called "map.csv", which must contain one 
% line for each edge, with columns x1, y1, x2, y2, and writes the result on 
% a file called "map_kementrality.csv". 
%
% A point is considered an intersection only if (x1,y1) or (x2,y2) match 
% exactly those in another road segment. Take care with approximations, e.g.,
% one coordinate being 1.5 on a segment and 1.49999999 on another.
%
% If the network is disconnected, the centrality is computed only on the
% largest connected component; ; on disconnected roads the column contains
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
%   * reg: regularization parameter as described in the paper; defaults 
%     to 10^(-8).
%   * weights: weight to use for each edge; defaults to 
%     exp(edge_length/max_length), as used in the paper
%   * filtered: wheter to apply filtering; defaults to true.
%   * parallel: whether to use parallel computation. You may set it to false
%     if you don't have the Matlab parallel toolbox installed, or if you want
%     a faster result for a small network. Otherwise just use the default
%     true.
%   
% To set parallel to false without touching the other defaults, use
% >> G = kementrality("map", [], [], [], false);
%
% Return value: a Matlab graph object, which includes G.Edges.kementrality.
%
% To plot the map inside Matlab, you can use
%
% >> plot(G, "XData", G.Nodes.x, "YData", G.Nodes.y, "EdgeCData", G.Edges.kementrality);

if not(exist('parallel', 'var')) || isempty(parallel)
    parallel = true;
end

if isa(basename, "char")
    G = convert_graphs(basename, true);
else
    G = basename;
end

if not(exist('weights', 'var')) || isempty(weights)
    G.Edges.Length = hypot(diff(G.Nodes.x(G.Edges.EndNodes)'), diff(G.Nodes.y(G.Edges.EndNodes)'))';
    weights = exp(-G.Edges.Length/max(G.Edges.Length));
end

kementrality = kementrality_old_negative_ldl(G, weights, parallel);

G.Edges.kementrality = kementrality;

if isa(basename, 'char')
    % creates a reduced Edges table that undoes the switching in
    % convert_graphs
    x1 = G.Edges.x1;
    x2 = G.Edges.x2;
    y1 = G.Edges.y1;
    y2 = G.Edges.y2;
    [x1(G.Edges.switched), x2(G.Edges.switched)] = deal(x2(G.Edges.switched), x1(G.Edges.switched));
    [y1(G.Edges.switched), y2(G.Edges.switched)] = deal(y2(G.Edges.switched), y1(G.Edges.switched));
    T2 = table(x1, y1, x2, y2, kementrality);

    T1 = readtable(strcat(basename, '.csv'));
    T1 = outerjoin(T1, T2, 'Type', 'left', 'MergeKeys', true);
    writetable(T1, strcat(basename, '_kementrality.csv'));
end

plotmeasure(G, kementrality);
