function G = kementrality_der(basename, reg, weights, relative, parallel)
% compute the derivative variant of the kemeny-based centrality from a csv file
%
% G = kementrality_der(basename, reg, weights, relative, parallel)
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
% >> G = kementrality_der("map", [], [], false)
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
% Additional parameters:
%   * reg: regularization parameter as described in the paper; defaults 
%     to 10^(-8).
%   * weights: weight to use for each edge; defaults to 
%     exp(edge_length/max_length), as used in the paper
%   * relative: wheter to scale with edge weight; defaults to true.
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

if not(exist('reg', 'var')) || isempty(reg)
    reg = 1e-8;
end
if not(exist('parallel', 'var')) || isempty(parallel)
    parallel = true;
end
if not(exist('relative', 'var')) || isempty(relative)
    relative = true;
end


G = convert_graphs(basename, true);

if not(exist('weights', 'var')) || isempty(weights)
    G.Edges.Length = hypot(G.Edges.x1-G.Edges.x2, G.Edges.y1-G.Edges.y2);
    weights = exp(-G.Edges.Length/max(G.Edges.Length));
end

[kementrality_der, ~] = kementrality_chol_der(G, reg, weights, parallel, relative);

G.Edges.kementrality_der = kementrality_der;
% sets up things 

% creates a reduced Edges table that undoes the switching
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

rescale = (kementrality_der - min(kementrality_der));
rescale = rescale / max(rescale);
rescale = 0.01 + 0.99 * rescale;
rescale = log(rescale);

% "reverse parula" seems the colormap that works better visually
cmap = colormap('parula');
colormap(cmap(end:-1:1,:));
plot(G, "EdgeCData", rescale, "XData", G.Nodes.x, "YData", G.Nodes.y, 'Marker', 'none');
colorbar;
