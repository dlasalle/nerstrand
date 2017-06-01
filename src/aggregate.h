/**
 * @file aggregate.h
 * @brief Functions revolving around vertex aggregation
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-07-23
 */




#ifndef NERSTRAND_AGGREGATE_H
#define NERSTRAND_AGGREGATE_H




#include "base.h"
#include "mgraph.h"
#include "objective.h"




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef struct aggregate_t {
  vtx_t * nvtxs;
  vtx_t ** match;
  vtx_t ** cmap;
  vtx_t ** fmap;
  int free_match; 
  int free_cmap; 
  int free_fmap; 
  int free_nvtxs;
  tid_t npar;
  dlthread_comm_t comm;
} aggregate_t;




/******************************************************************************
* SERIAL FUNCTION PROTOTYPES **************************************************
******************************************************************************/


#define init_aggregate __nerstrand_init_aggregate
/**
 * @brief Initialize an aggregate structure
 *
 * @param aggregate The structure to initialize
 *
 * @return A pointer to the initialized structure
 */
aggregate_t * init_aggregate(
    aggregate_t * aggregate);


#define setup_aggregate __nerstrand_setup_aggregate
/**
 * @brief Setup an aggregate structure
 *
 * @param nvtxs The number of vertices per thread (optionally leave NULL if all
 * vectors are supplied)
 * @param cnvtxs The number of coarse vertices per thread (optionally leave
 * NULL and nvtxs will be used)
 * @param match The match vectors (optionally leave NULL and they will be
 * allocated) 
 * @param cmap The cmap vectors (optionally leave NULL and they will be
 * allocated)
 * @param fmap The fmap vectors (optionally leave NULL and they will be
 * allocated at size of cnvtxs, or nvtxs if cnvtxs is NULL)
 * @param nthreads The number of threads that will use the structure
 *
 * @return The setup aggregate structure 
 */
aggregate_t * setup_aggregate(
    vtx_t const * nvtxs, 
    vtx_t * cnvtxs, 
    vtx_t ** match, 
    vtx_t ** cmap, 
    vtx_t ** fmap, 
    tid_t const nthreads);


#define free_aggregate __nerstrand_free_aggregate
/**
 * @brief Free an aggregate structure and its associated memory
 *
 * @param aggregate The structure to free
 *
 * @return !0 on success
 */
int free_aggregate(
    aggregate_t * aggregate);


#define aggregate_graph __nerstrand_aggregate_graph
/**
 * @brief Create an aggregation from the supplied graph using the parameters
 * specified by the supplied objective.
 *
 * @param objective The objective containing aggregation parameters
 * @param mgraph The mgraph to aggregate
 *
 * @return The new aggregation
 */
aggregate_t * aggregate_graph(
    objective_t * objective, 
    mgraph_t const * mgraph);


#define check_aggregate __nerstrand_check_aggregate
/**
 * @brief Verify the sanity of an aggregate
 *
 * @param aggregate The aggregate to verify
 * @param graph The aggregated graph
 *
 * @return 1 on success, 0 on invalid aggregation
 */
int check_aggregate(
    aggregate_t const * aggregate, 
    graph_t const * graph);




/******************************************************************************
* PARALLEL FUNCTION PROTOTYPES ************************************************
******************************************************************************/


#define par_setup_aggregate __nerstrand_par_setup_aggregate
/**
 * @brief Setup an aggregate structure -- to be called by each thread in a
 * parallel region. It should be noted that if the parameter mycnvtxs is passed
 * as NULL_VTX, then aggregate->nvtxs must be set later.
 *
 * @param mynvtxs The number of vertices owned by the current thread (can be
 *   NULL_VTX if arrays are supplied)
 * @param mycnvtxs The number of coarse vertices owned by the current thread
 *   (NULL_VTX if it is unknown and mynvtxs will be used)
 * @param match The match vector to use for aggregation (can be NULL)
 * @param cmap The cmap vector to use for aggregation (can be NULL)
 * @param fmap The fmap vector to use for aggregation (can be NULL, will be
 *   allocated at size of mycnvtxs or nvtxs)
 * @param comm The thread communicator.
 *
 * @return The newly setup aggregate structure.
 */
aggregate_t * par_setup_aggregate(
    vtx_t mynvtxs,
    vtx_t mycnvtxs,
    vtx_t * match,
    vtx_t * cmap,
    vtx_t * fmap,
    dlthread_comm_t comm);


#define par_free_aggregate __nerstrand_par_free_aggregate
/**
 * @brief Free and aggregate structure -- to be called by each thread in a
 * parallel region
 *
 * @param aggregate The aggregate structure to free
 *
 * @return !0 on success
 */
int par_free_aggregate(
    aggregate_t * aggregate);


#define par_aggregate_graph __nerstrand_par_aggregate_graph
/**
 * @brief Create an aggregation from the supplied graph using the parameters
 * specified by the supplied objective.
 *
 * @param objective The objective containing aggregation parameters
 * @param mgraph The mgraph to aggregate
 *
 * @return The new aggregation
 */
aggregate_t * par_aggregate_graph(
    objective_t * objective,
    mgraph_t const * mgraph);




#endif
