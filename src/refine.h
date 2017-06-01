/**
 * @file refine.h
 * @brief Function prototypes and types for refinement
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-10-08
 */




#ifndef NERSTRAND_REFINE_H
#define NERSTRAND_REFINE_H




#include "base.h"
#include "mgraph.h"
#include "cluster.h"
#include "objective.h"




/******************************************************************************
* STRUCTS *********************************************************************
******************************************************************************/


typedef struct move_t {
  vtx_t vid;
  cid_t to, from;
} move_t;


typedef struct update_t {
  cid_t to, from;
  wgt_t wgt;
  vtx_t vid;
} update_t;




/******************************************************************************
* SERIAL FUNCTION PROTOTYPES **************************************************
******************************************************************************/


#define explicit_refine_graph __nerstrand_explicit_refine_graph
/**
 * @brief Refine a clustering on a graph.
 *
 * @param clustering The clustering to refine.
 * @param mgraph The clustered graph.
 * @param type The type of refinement to perform.
 * @param npasses The number of refinement passes to perform.
 * @param empty Whether to allow clusters to be emptied (removed).
 * @param seed The random seed to use.
 * @param refstats A structure for accumulating refinement statistics.
 *
 * @return The updated clsutering
 */
clustering_t * explicit_refine_graph(
    clustering_t * clustering, 
    mgraph_t const * mgraph, 
    nerstrand_reftype_t type,
    size_t npasses,
    int empty,
    unsigned int seed,
    refinement_stats_t * refstats);


#define refine_graph __nerstrand_refine_graph
/**
 * @brief Refine a clustering on a graph
 *
 * @param objective The objective with parameters to be used for refinement.
 * @param clustering The clustering to refine.
 * @param mgraph The clustered graph.
 * @return The updated clsutering.
 */
clustering_t * refine_graph(
    objective_t * objective, 
    clustering_t * clustering,
    mgraph_t const * mgraph);




/******************************************************************************
* PARALLEL FUNCTION PROTOTYPES ************************************************
******************************************************************************/


#define par_explicit_refine_graph __nerstrand_par_explict_refine_graph
/**
 * @brief Refine a clustering on a graph -- to be called by each thread in a
 * a parallel region.
 *
 * @param clustering The clustering to refine.
 * @param mgraph The clustered graph.
 * @param type The type of refinement to perform.
 * @param npasses The number of refinement passes to perform.
 * @param empty Whether to allow clusters to be emptied (removed).
 * @param seed The random seed to use.
 * @param refstats A structure for accumulating refinement statistics.
 * @param comm The thread communication structure.
 *
 * @return The updated clsutering.
 */
clustering_t * par_explicit_refine_graph(
    clustering_t * clustering, 
    mgraph_t const * mgraph, 
    nerstrand_reftype_t type,
    size_t npasses,
    int empty,
    unsigned int seed,
    refinement_stats_t * refstats,
    dlthread_comm_t const comm);


#define par_refine_graph __nerstrand_par_refine_graph
/**
 * @brief Refine a clustering on a graph -- to be called by each thread in
 * a parallel region.
 *
 * @param objective The objective containing refinement parameters.
 * @param clustering The clustering to refine.
 * @param mgraph The clustered graph.
 *
 * @return The updated graph.
 */
clustering_t * par_refine_graph(
    objective_t * objective, 
    clustering_t * clustering, 
    mgraph_t const * mgraph);




#endif
