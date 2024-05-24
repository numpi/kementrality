Matlab code to compute the Kemeny-based centrality described in https://doi.org/10.1137/22M1486728

Usage:

1. prepare (using QGis or some other tool) a file `somename.csv` (where you can put any name of your choice instead of `somename`) with at least 4 columns `x1`, `y1`, `x2`, `y2`, which identify the beginning / end of each road.
2. run in Matlab 

    G = kementrality('somename')

3. a file `somename_kementrality.csv` is created; this contains the smae columns as `somename.csv` plus one column named `kementrality`, which contains the values of the measure. This column may contain NaN if some values could not be computed; this happens for instance if the network is disconnected or has self-loop edges (zero-length roads).
4. You may then re-import `somename-kementrality.csv` into QGis or your tool of choice.
5. Optionally, you may plot the result directly in Matlab with

    plot(G, "XData", G.Nodes.x, "YData", G.Nodes.y, "EdgeCData", G.Edges.kementrality)

or interact with it; the function returns a Matlab `Graph` object.

For more detailed information, type into Matlab

    doc kementrality
