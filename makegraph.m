function G = makegraph(Edges)
% Create a Matlab graph from the given table
%
% Input: a table with columns x1, y1, x2, y2, and (optionally) some extra edge-specific fields
%
% This is the "primal" graph with roads as its edges, not as its nodes.
%
% Includes some cleanup: removes duplicate edges and loops.
%
% You may plot it with this command
% plot(G, 'XData', G.Nodes.x, 'YData', G.Nodes.y)
%
% or compute centralities with, like,
% centrality(G, 'pagerank');
%
% (keep in mind that centrality(Toscana, 'betweenness' / 'closeness') may
% take a lot of time, though.)
%
% Remark: we keep the field 'connectivity' in the EdgeTable, but it contains
% some wrong entries since we removed a few edges.

Nodes_array = unique([Edges{:,{'x1','y1'}}; Edges{:,{'x2','y2'}}], 'rows');
Nodes = table(Nodes_array(:,1), Nodes_array(:,2), 'VariableNames', {'x2', 'y2'});

% we add integer ids to use for nodes, instead of using coordinate pairs, 
% for ease of use and probably also performance reasons.

% add id2 (we add id2 before id1 so that they come up in the right order in
% the final result)
Nodes.id2 = uint32(1:height(Nodes))';
Edges = innerjoin(Nodes,Edges);

% add id1
Nodes.Properties.VariableNames = {'x1', 'y1', 'id1'};
Edges = innerjoin(Nodes,Edges);

% remove loops -- moved below inside simplify()
%Edges(Edges.id1 == Edges.id2,:) = [];

G = graph(Edges.id1, Edges.id2, Edges);
Nodes.Properties.VariableNames = {'x', 'y', 'id'};
Nodes = movevars(Nodes, 'id', 'before', 1);
G.Nodes = Nodes;

% Matlab reorders edges in increasing order, so building the graph
% may have switched id1 and id2. In this case, we need to switch
% x1, x2, y1, y2. We will also add a column to keep track of these switches

G.Edges.switched = (G.Edges.id1 ~= G.Edges.EndNodes(:,1));
[G.Edges.x1(G.Edges.switched), G.Edges.x2(G.Edges.switched)] = deal(G.Edges.x2(G.Edges.switched), G.Edges.x1(G.Edges.switched));
[G.Edges.y1(G.Edges.switched), G.Edges.y2(G.Edges.switched)] = deal(G.Edges.y2(G.Edges.switched), G.Edges.y1(G.Edges.switched));

% id1 and id2 should now be switched, too; then they become a no longer
% necessary copy of EndNodes, so we just delete them.

G.Edges.id1 = [];
G.Edges.id2 = [];

% We delete id as well, since we won't need it and it may cause confusion
% as it is almost a copy of the node number
G.Nodes.id = [];

G = simplify(G, 'first', 'omitselfloops');
