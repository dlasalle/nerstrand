/**
 * @file stats.h
 * @brief Types and function prototypes for keeping track of stats during
 * clustering
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-10-08
 */




#ifndef NERSTRAND_STATS_H
#define NERSTRAND_STATS_H




#include "base.h"




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef struct loadbalance_stat_t {
  real_t adj_imbalance;
  real_t vtx_imbalance;
  real_t bnd_imbalance;
  real_t nbr_imbalance;
} loadbalance_stat_t;


typedef struct aggregation_stat_t {
  vtx_t nsupernodes;
  vtx_t nlonenodes;
  real_t coarsenrate;
  vtx_t maxnodesize;
  real_t avgnodesize;
  vtx_t nbrokenmatches;
} aggregation_stat_t;


typedef struct refinement_stat_t {
  vtx_t nmoves;
  vtx_t nunmoves;
  vtx_t nbnd;
  vtx_t actgain;
  vtx_t pergain;
  size_t npasses;
} refinement_stat_t;


typedef struct initialclustering_stats_t {
  size_t nclusterings;
  size_t nruns;
  size_t * lvls;
  cid_t * nclusters;
  cid_t * nbclusters;
  mod_t * mods;
  vtx_t * nvtxs;
  adj_t * nedges;
  wgt_t * tadjwgt;
} initialclustering_stats_t;




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


/* loadbalance_stat_t */
#define DLMEM_PREFIX loadbalance_stat
#define DLMEM_TYPE_T loadbalance_stat_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLLIST_PREFIX loadbalance_stat
#define DLLIST_TYPE_T loadbalance_stat_t
#define DLLIST_LINKED
#include "dllist_headers.h"
#undef DLLIST_LINKED
#undef DLLIST_TYPE_T
#undef DLLIST_PREFIX


/* aggregation_stat_t */
#define DLMEM_PREFIX aggregation_stat
#define DLMEM_TYPE_T aggregation_stat_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLLIST_PREFIX aggregation_stat
#define DLLIST_TYPE_T aggregation_stat_t
#define DLLIST_LINKED
#include "dllist_headers.h"
#undef DLLIST_LINKED
#undef DLLIST_TYPE_T
#undef DLLIST_PREFIX


/* refinement_stat_t */
#define DLMEM_PREFIX refinement_stat
#define DLMEM_TYPE_T refinement_stat_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLLIST_PREFIX refinement_stat
#define DLLIST_TYPE_T refinement_stat_t
#define DLLIST_LINKED
#include "dllist_headers.h"
#undef DLLIST_LINKED
#undef DLLIST_TYPE_T
#undef DLLIST_PREFIX


/* initialclustering_stat_t */
#define DLMEM_PREFIX initialclustering_stats
#define DLMEM_TYPE_T initialclustering_stats_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


/* rename the lists */
#define loadbalance_stats_t loadbalance_stat_list_t
#define loadbalance_stats_create loadbalance_stat_list_create
#define loadbalance_stats_push loadbalance_stat_list_push
#define loadbalance_stats_pop loadbalance_stat_list_pop
#define loadbalance_stats_get loadbalance_stat_list_get
#define loadbalance_stats_replace loadbalance_stat_list_replace
#define loadbalance_stats_free loadbalance_stat_list_free
#define aggregation_stats_t aggregation_stat_list_t
#define aggregation_stats_create aggregation_stat_list_create
#define aggregation_stats_push aggregation_stat_list_push
#define aggregation_stats_pop aggregation_stat_list_pop
#define aggregation_stats_free aggregation_stat_list_free
#define refinement_stats_t refinement_stat_list_t
#define refinement_stats_create refinement_stat_list_create
#define refinement_stats_push refinement_stat_list_push
#define refinement_stats_pop refinement_stat_list_pop
#define refinement_stats_free refinement_stat_list_free




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define init_loadbalance_stat __nerstrand_init_loadbalance_stat
/**
 * @brief Initialize a load balancing statistics structure.
 *
 * @param stat The load balancing statistics structure.
 */
void init_loadbalance_stat(
    loadbalance_stat_t * stat);


#define init_aggregation_stat __nerstrand_init_aggregation_stat
/**
 * @brief Initialize a aggregation statistics structure.
 *
 * @param stat The aggregation statistics structure.
 */
void init_aggregation_stat(
    aggregation_stat_t * stat);


#define init_refinement_stat __nerstrand_init_refinement_stat
/**
 * @brief Initialize a refinement statistics structure.
 *
 * @param stat The refinement statistics structure.
 */
void init_refinement_stat(
    refinement_stat_t * stat);


#define init_initialclustering_stats __nerstrand_init_initialclustering_stats
/**
 * @brief Initialize an initial clustering statistics structure.
 *
 * @param stats The initial clustering statistics structure.
 */
void init_initialclustering_stats(
    initialclustering_stats_t * stats);


#define initialclustering_stats_create \
    __nerstrand_initialclustering_stats_create
/**
 * @brief Allocate and initialize an initial clustering statistics structure.
 *
 * @param nruns The number of runs (clusterings to be generated) to be 
 *   performed.
 * @param ncuts The numbe rof cuts (initial clusterings) to be made per run.
 *
 * @return The initial clustering statistics structure.
 */
initialclustering_stats_t * initialclustering_stats_create(
    size_t nruns, 
    size_t ncuts);


#define initialclustering_stats_free __nerstrand_initialclustering_stats_free
/**
 * @brief Free an initial clustering statistics structure and associated
 *   memory.
 *
 * @param stats The structure to free.
 */
void initialclustering_stats_free(
    initialclustering_stats_t * stats);




#endif
