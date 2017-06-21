Nerstrand
=========

Nerstrand is a multi-threaded multilevel graph clustering tool for generating
clusterings with high modularity. It supports both finding a specified number
of clusters/communities as well as detecting the number of
clusters/communities. 


Packages
--------

Built using https://bytepackager.com.

<a href="https://github.com/dlasalle/nerstrand/releases">
  <img src="https://bytepackager.com/badge/xoguahkIzHq9O1OgT3opySIVJ3A/status.svg"/>
</a>

<a href="https://github.com/dlasalle/nerstrand/releases">
  <img src="https://bytepackager.com/badge/YgTwmFEn4A5u2i8Bd-kGeuJu86c/status.svg"/>
</a>

<a href="https://github.com/dlasalle/nerstrand/releases">
  <img src="https://bytepackager.com/badge/Tx36QrcDiK4cD6LyN2yR25YhwBU/status.svg"/>
</a>

<a href="https://github.com/dlasalle/nerstrand/releases">
  <img src="https://bytepackager.com/badge/IrUcbkkLLSx-5XLjyU3Byr_YNIE/status.svg"/>
</a>


The 'nerstrand' Executable
--------------------------

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
    


Unit Tests
----------

<a href="https://travis-ci.org/dlasalle/nerstrand">
  <img src="https://travis-ci.org/dlasalle/nerstrand.svg?branch=master"/>
</a>



The Nerstrand API
-----------------

The file [nerstrand.h](@ref nerstrand.h) is the header that should be included
by external programs wishing link to Nerstrand. It contains three functions for
creating clusterings.


    nerstrand_cluster_anyway()

    nerstrand_cluster_kway()   

    nerstrand_cluster_explicit()

    nerstrand_init_options()


These functions are documented [here](@ref nerstrand.h).

