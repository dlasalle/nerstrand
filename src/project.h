/**
 * @file project.h
 * @brief Function prototypes for projecting a clustering
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2014-03-03
 */




#ifndef NERSTRAND_PROJECT_H
#define NERSTRAND_PROJECT_H




#include "base.h"
#include "cluster.h"
#include "mgraph.h"
#include "objective.h"




/******************************************************************************
* SERIAL FUNCTION PROTOTYPES **************************************************
******************************************************************************/


#define project_clustering __nerstrand_project_clustering
/**
 * @brief Project a clustering from a coarse graph to the fine graph one level
 * down
 *
 * @param objective The objective containing parameters for the projection
 * @param clustering The clustering to project
 * @param fmgraph The fine graph
 * @param cmgraph The coarse graph
 *
 * @return The updated clustering
 */
clustering_t * project_clustering(
    objective_t * objective, 
    clustering_t * clustering, 
    mgraph_t * fmgraph, 
    mgraph_t const * cmgraph);




/******************************************************************************
* PARALLEL FUNCTION PROTOTYPES ************************************************
******************************************************************************/


#define par_project_clustering __nerstrand_par_project_clustering
/**
 * @brief Project a clustering from a coarse graph to the fine graph one level
 * down -- to be called by each thread from within a parallel region
 *
 * @param objective The objecitve containing parameters for the projection
 * @param clustering The clustering to project
 * @param fmgraph The fine graph
 * @param cmgraph The coarse graph
 *
 * @return The updated clustering 
 */
clustering_t * par_project_clustering(
    objective_t * objective, 
    clustering_t * clustering, 
    mgraph_t * fmgraph, 
    mgraph_t const * cmgraph);




#endif
