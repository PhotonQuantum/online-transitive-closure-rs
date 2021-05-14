# Online Transitive Closure (SCC)

Let G=<V, E> be a directed graph, G*=<V, E*> its transitive closure. Let E* be represented by incidence matrices.
Suppose edges are inserted in G one at a time. We consider the problem of efficiently updating G* each time an edge
is inserted.

This algorithm for updating G* in the case of edge insertions requires O(|Enew||V|) time for 
q consecutive insertions.

The algorithm especially yield better time complexity for graphs with Eold << |E*| and for graphs
with relatively small components.

This crate only implements adding edges and G_c & G*.

Please refer to the original paper for details:

La, J. & Leeuwen, Poutr & Informatica, Vakgroep & Nthrlan, The. (2001). Maintenance Of Transitive Closures And Transitive Reductions Of Graphs. 