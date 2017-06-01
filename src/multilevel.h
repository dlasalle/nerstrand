/**
 * @file multilevel.h
 * @brief Functions for multilevel graph clustering
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2014-03-03
 */




#ifndef NERSTRAND_MULTILEVEL_H
#define NERSTRAND_MULTILEVEL_H




#include "base.h"
#include "objective.h"
#include "graph.h"
#include "coarsen.h"
#include "uncoarsen.h"
#include "initialclustering.h"




/******************************************************************************
* SERIAL FUNCTION PROTOTYPES **************************************************
******************************************************************************/


#define cluster_graph __nerstrand_cluster_graph
/**
 * @brief Generate a clustering of a given graph 
 *
 * @param objective The objective containing parameters for the clustering
 * @param graph The graph to cluster
 *
 * @return The generated clustering
 */
clustering_t * cluster_graph(
    objective_t * objective, 
    graph_t const * graph);




/******************************************************************************
* PARALLEL FUNCTION PROTOTYPES ************************************************
******************************************************************************/


#define par_cluster_graph __nerstrand_par_cluster_graph
/**
 * @brief Generate a clustering of a given graph in parallel
 *
 * @param objective The objective containing parameters for the clustering 
 * @param graph The graph to cluster
 *
 * @return The clustered graph
 */
clustering_t * par_cluster_graph(
    objective_t * objective, 
    graph_t const * graph);




#endif
