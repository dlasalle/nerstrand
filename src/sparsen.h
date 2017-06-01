/**
 * @file sparsen.h
 * @brief Edge sparsening functions
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-07-25
 */




#ifndef NERSTRAND_SPARSEN_H
#define NERSTRAND_SPARSEN_H




#include "base.h"
#include "objective.h"
#include "graph.h"




/******************************************************************************
* SERIAL FUNCTION PROTOTYPES **************************************************
******************************************************************************/


#define sparsen_graph __nerstrand_sparsen_graph
/**
 * @brief Lower the edge density of a graph
 *
 * @param objective The objective containing sparsen parameters
 * @param graph The graph to sparsen
 * @param targetnedges The target number of remaining edges
 *
 * @return NERSTRAND_SUCCESS, or an error code if unsuccessful
 */
int sparsen_graph(
    objective_t * objective, 
    graph_t * graph,
    adj_t targetnedges);




#endif
