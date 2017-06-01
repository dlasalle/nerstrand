/**
 * @file coarsen.h
 * @brief Function prototypes for coarsening a graph
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2014-03-03
 */




#ifndef NERSTRAND_COARSEN_H
#define NERSTRAND_COARSEN_H




#include "base.h"
#include "mgraph.h"
#include "objective.h"




/******************************************************************************
* SERIAL FUNCTION PROTOTYPES **************************************************
******************************************************************************/


#define coarsen_graph __nerstrand_coarsen_graph
/**
 * @brief Generate a coarser graph using the methods specified by the objective
 *
 * @param objective The objective structure specifying what methods to use for
 * coarsening
 * @param mgraph The multi-level graph to coarsen
 *
 * @return The coarser graph
 */
mgraph_t * coarsen_graph(
    objective_t * objective, 
    mgraph_t * mgraph);




/******************************************************************************
* PARALLEL FUNCTION PROTOTYPES ************************************************
******************************************************************************/


#define par_coarsen_graph __nerstrand_par_coarsen_graph
/**
 * @brief Generate a coarser graph using the methods specified by the objective
 * -- to be called by each thread in a parallel region
 *
 * @param objective The objective structure specifying what methods to use
 * @param mgraph The multi-level graph to coarsen
 *
 * @return The coarser graph
 */
mgraph_t * par_coarsen_graph(
    objective_t * objective, 
    mgraph_t * mgraph);




#endif
