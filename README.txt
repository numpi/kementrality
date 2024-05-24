Matlab code to compute the Kemeny-based centrality described in https://doi.org/10.1137/22M1486728

Usage:

1. prepare (using QGIS or some other tool) a file `somename.csv` file with at least 4 columns `x1`, `y1`, `x2`, `y2`, which identify the beginning / end of each road.
2. run in Matlab 

    >> G = kementrality('somename')

3. a file `somename_kementrality.csv` is created; this contains the smae columns as `somename.csv` plus one column named `kementrality`.
4. Optionally, you may plot the result directly in Matlab with

    >> plot(G, "XData", G.Nodes.x, "YData", G.Nodes.y, "EdgeCData", G.Edges.kementrality)


For more detailed information, type

    >> doc kementrality
