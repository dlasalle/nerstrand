/**
 * @file newman.h
 * @brief Function prototypes for agglormerative clustering
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-10-24
 */




#ifndef NERSTRAND_NEWMAN_H
#define NERSTRAND_NEWMAN_H




#include "cluster.h"
#include "graph.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define generate_newman_clustering __nerstrand_generate_newman_clustering
/**
 * @brief Generate a new clustering using CNM.
 *
 * @param graph The graph to cluster.
 *
 * @return The generated clustering.
 */
clustering_t * generate_newman_clustering(
    graph_t const * graph);




#endif
