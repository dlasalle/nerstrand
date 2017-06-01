/**
 * @file initialclustering.h
 * @brief Function prototypes for initial clustering
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-10-08
 */




#ifndef NERSTRAND_INITIALCLUSTERING_H
#define NERSTRAND_INITIALCLUSTERING_H




#include "base.h"
#include "cluster.h"
#include "mgraph.h"
#include "objective.h"




/******************************************************************************
* SERIAL FUNCTION PROTOTYPES **************************************************
******************************************************************************/


#define generate_kclustering __nerstrand_generate_kclustering
/**
 * @brief Generate a k-way clustering serially.
 *
 * @param objective The objective specifiying clustering parameters.
 * @param mgraph The graph to cluster.
 *
 * @return A clustering of the input graph. 
 */
clustering_t * generate_kclustering(
    objective_t * objective,
    mgraph_t * mgraph);


#define generate_aclustering __nerstrand_generate_aclustering
/**
 * @brief Generate an any-way clustering serially.
 *
 * @param objective The objective specifiying clustering parameters.
 * @param mgraph The graph to cluster.
 *
 * @return A clustering of the input graph. 
 */
clustering_t * generate_aclustering(
    objective_t * objective,
    mgraph_t * mgraph);




/******************************************************************************
* PARALLEL FUNCTION PROTOTYPES ************************************************
******************************************************************************/


#define par_generate_kclustering __nerstrand_par_generate_kclustering
/**
 * @brief Generate an k-way clustering in parallel.
 *
 * @param objective The objective specifiying clustering parameters.
 * @param mgraph The graph to cluster.
 *
 * @return A clustering of the input graph. 
 */
clustering_t * par_generate_kclustering(
    objective_t * objective,
    mgraph_t * mgraph);


#define par_generate_aclustering __nerstrand_par_generate_aclustering
/**
 * @brief Generate an any-way clustering in parallel.
 *
 * @param objective The objective specifiying clustering parameters.
 * @param mgraph The graph to cluster.
 *
 * @return A clustering of the input graph. 
 */
clustering_t * par_generate_aclustering(
    objective_t * objective,
    mgraph_t * mgraph);




#endif
