# {#mainpage}

Nerstrand
=============================

Nerstrand is a multi-threaded multilevel graph clustering tool for generating
clusterings with high modularity. It supports both finding a specified number
of clusters/communities as well as detecting the number of
clusters/communities. 


The 'nerstrand' Executable
-----------------------------

Nerstrand can then be used to cluster a graph as follows:

    nerstrand test.graph


Most users will actually want the cluster assignment and not just runtime
statistics, and should append the name of the output file.

    nerstrand test.graph test.cluster 


The number of clusters to generate can be specified with the '-c' or 
'--clusters' option.

    nerstrand -c16 test.graph test.cluster


Many other runtime parameters (including experimental ones) are available. Use 
the '-h' option to view these.

    nerstrand -h
    


The Nerstrand API
-------------------------------

The file [nerstrand.h](@ref nerstrand.h) is the header that should be included
by external programs wishing link to Nerstrand. It contains three functions for
creating cluterings.


    nerstrand_cluster_anyway()

    nerstrand_cluster_kway()   

    nerstrand_cluster_explicit()

    nerstrand_init_options()


These functions are documented [here](@ref nerstrand.h).

