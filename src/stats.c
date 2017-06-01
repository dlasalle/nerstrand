/**
 * @file stats.c
 * @brief Functions for clustering statistic tracking
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-10-08
 */




#ifndef NERSTRAND_STATS_C
#define NERSTRAND_STATS_C




#include "stats.h"




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


/* loadbalance_stat_t */
#define DLMEM_PREFIX loadbalance_stat
#define DLMEM_TYPE_T loadbalance_stat_t
#define DLMEM_DLTYPE DLTYPE_STRUCT
#define DLMEM_INITFUNCTION init_loadbalance_stat
#include "dlmem_funcs.h"
#undef DLMEM_DLTYPE
#undef DLMEM_INITFUNCTION
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLLIST_PREFIX loadbalance_stat
#define DLLIST_TYPE_T loadbalance_stat_t
#define DLLIST_DLTYPE DLTYPE_STRUCT
#define DLLIST_LINKED
#include "dllist_funcs.h"
#undef DLLIST_DLTYPE
#undef DLLIST_LINKED
#undef DLLIST_TYPE_T
#undef DLLIST_PREFIX


/* aggregation_stat_t */
#define DLMEM_PREFIX aggregation_stat
#define DLMEM_TYPE_T aggregation_stat_t
#define DLMEM_DLTYPE DLTYPE_STRUCT
#define DLMEM_INITFUNCTION init_aggregation_stat
#include "dlmem_funcs.h"
#undef DLMEM_DLTYPE
#undef DLMEM_INITFUNCTION
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLLIST_PREFIX aggregation_stat
#define DLLIST_TYPE_T aggregation_stat_t
#define DLLIST_DLTYPE DLTYPE_STRUCT
#define DLLIST_LINKED
#include "dllist_funcs.h"
#undef DLMEM_DLTYPE
#undef DLLIST_LINKED
#undef DLLIST_TYPE_T
#undef DLLIST_PREFIX


/* refinement_stat_t */
#define DLMEM_PREFIX refinement_stat
#define DLMEM_TYPE_T refinement_stat_t
#define DLMEM_DLTYPE DLTYPE_STRUCT
#define DLMEM_INITFUNCTION init_refinement_stat
#include "dlmem_funcs.h"
#undef DLMEM_DLTYPE
#undef DLMEM_INITFUNCTION
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLLIST_PREFIX refinement_stat
#define DLLIST_TYPE_T refinement_stat_t
#define DLLIST_LINKED
#define DLLIST_DLTYPE DLTYPE_STRUCT
#include "dllist_funcs.h"
#undef DLMEM_DLTYPE
#undef DLLIST_LINKED
#undef DLLIST_TYPE_T
#undef DLLIST_PREFIX




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void init_loadbalance_stat(
    loadbalance_stat_t * const stat)
{
  stat->adj_imbalance = 0.0;
  stat->vtx_imbalance = 0.0;
  stat->bnd_imbalance = 0.0;
  stat->nbr_imbalance = 0.0;
}


void init_aggregation_stat(
    aggregation_stat_t * const stat)
{
  stat->nsupernodes = 0;
  stat->nlonenodes = 0;
  stat->coarsenrate = 0.0;
  stat->maxnodesize = 0;
  stat->avgnodesize = 0;
  stat->nbrokenmatches = 0;
}


void init_refinement_stat(
    refinement_stat_t * const stat)
{
  stat->nmoves = 0;
  stat->nunmoves = 0;
  stat->nbnd = 0;
  stat->actgain = 0.0;
  stat->pergain = 0.0;
  stat->npasses = 0;
}


void init_initialclustering_stats(
    initialclustering_stats_t * const stats)
{
  stats->nclusterings = 0;
  stats->nruns = 0;
  stats->lvls = NULL;
  stats->mods = NULL;
  stats->nvtxs = NULL;
  stats->nedges = NULL;
  stats->nclusters = NULL;
  stats->nbclusters = NULL;
}


initialclustering_stats_t * initialclustering_stats_create(
    size_t const nruns, 
    size_t const ncuts)
{
  initialclustering_stats_t * stats;
  
  stats = (initialclustering_stats_t*) \
          malloc(sizeof(initialclustering_stats_t));
  stats->nclusterings = nruns*ncuts;
  stats->nruns = nruns;
  stats->lvls = size_alloc(stats->nclusterings);
  stats->mods = mod_alloc(stats->nclusterings);
  stats->nvtxs = vtx_alloc(stats->nclusterings);
  stats->nedges = adj_alloc(stats->nclusterings);
  stats->tadjwgt = wgt_alloc(stats->nclusterings);
  stats->nclusters = cid_alloc(stats->nclusterings);
  stats->nbclusters = cid_alloc(stats->nruns);
  return stats;
}


void initialclustering_stats_free(
    initialclustering_stats_t * stats)
{
  dl_free(stats->lvls);
  dl_free(stats->mods);
  dl_free(stats->nvtxs);
  dl_free(stats->nedges);
  dl_free(stats->tadjwgt);
  dl_free(stats->nclusters);
  dl_free(stats->nbclusters);
  dl_free(stats);
}




#endif
