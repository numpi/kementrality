function [G, H] = convert_graphs(basename, check_connectivity)
% read basename.csv and produce auxiliary files basename_nodes.csv,
% basename_edges.csv, basename_dual_edges.csv

% if check_connectivity==false, the graph is assumed to be connected.

if not(exist('check_connectivity', 'var'))
    check_connectivity = true;
end

f = strcat(basename, '.csv');
Edges = readtable(f);
G = makegraph(Edges);
if check_connectivity
    c = conncomp(G);
    G = subgraph(G, c==mode(c));
end
H = dualgraph(G);
writetable(H.Edges, strcat(basename, '_dual_edges.csv'));
writetable(H.Nodes, strcat(basename, '_edges.csv'));
writetable(G.Nodes, strcat(basename, '_nodes.csv'));
