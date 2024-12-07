function H = dualgraph(G)
% construct edge-graph
%
% Given a graph with roads as its edges, construct a graph with roads as
% its points.
%
% You may plot the graph with 
% plot(H, 'XData', (H.Nodes.x1+H.Nodes.x2)/2, 'YData', (H.Nodes.y1+H.Nodes.y2)/2)
%
% Centrality of the dual graph H plotted on the primal:
% H = dualgraph(G);
% C = centrality(H, 'pagerank');
% plot(G, 'XData', G.Nodes.x, 'YData', G.Nodes.y, 'EdgeCData', C, 'Marker', 'none', 'LineWidth', 3);

A = incidence(G);
A = A'; %it is much faster to search a sparse matrix column-wise than row-wise, so we transpose it

fprintf('Created incidence matrix.\n');

AnglesWithHorizontal = atan2(G.Edges.y2 - G.Edges.y1, G.Edges.x2 - G.Edges.x1);

dualEdges = [];
block = [];
blocksize = 10000; % we grow dualEdges block by block, to avoid resizing it continuously. This has a huge impact on performance.
angles = [];
angleBlock = [];

progress = max(size(A,2) / 100, 1000);

for i = 1:size(A,2)
    % pushes into dualEdges all pairs of indices touching node i
    [edges_touching_i, ~, signs] = find(A(:,i));
    angles_in_i = AnglesWithHorizontal(edges_touching_i) + pi*(signs==1); % the angle becomes its opposite if we consider it wrt (x2,y2) rather than (x1,y1)    
    k = length(edges_touching_i);
    for h = 1:k-1
        block = [block; [repmat(edges_touching_i(h),k-h,1) edges_touching_i(h+1:k)]];
        angleBlock = [angleBlock; mod(angles_in_i(h+1:k) - angles_in_i(h), 2*pi) - pi]; % reducing angles_in_i(h) + angles_in_i(h+1:k) in [-pi,pi]
    end
    if mod(i, progress) == 0
        fprintf('%d/%d ', i, size(A,2));
        if mod(i, 10*progress) == 0
            fprintf('\n');
        end
    end
    if size(block,1) > blocksize
        dualEdges = [dualEdges; block];
        angles = [angles; angleBlock];
        block = [];
        angleBlock = [];
    end
end
dualEdges = [dualEdges; block];
angles = [angles; angleBlock];
fprintf('\n');
dualEdges2 = table(dualEdges, angles, 'VariableNames', {'EndNodes', 'ChangeOfDirection'});

H = graph(dualEdges2, G.Edges);
