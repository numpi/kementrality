function plotmeasure(G, measure)
% plots an edge centrality measure on graph G
%
% measure should be a vector with length numedges(G). 
% It will get automatically rescaled.
%
% G should contain coordinates G.Nodes.x and G.Nodes.y.

rescale = (measure - min(measure));
rescale = rescale / max(rescale);
rescale = 0.01 + 0.99 * rescale;
rescale = log(rescale);

% "reverse parula" seems the colormap that works better visually
cmap = colormap('parula');
colormap(cmap(end:-1:1,:));
plot(G, "EdgeCData", rescale, "XData", G.Nodes.x, "YData", G.Nodes.y, 'Marker', 'none');
colorbar;
