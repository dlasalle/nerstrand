/**
 * @file nerstrand.c
 * @brief Top level functions for nerstrand
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-07-23
 */




#ifndef NERSTRAND_C
#define NERSTRAND_C




#include "base.h"
#include "graph.h"
#include "cluster.h"
#include "objective.h"
#include "multilevel.h"
#include "io.h"





/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef struct arg_t {
  objective_t * objective;
  graph_t const * graph;
  clustering_t * cluster;
} arg_t;




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static void __thread_cluster(
    void * const ptr)
{
  size_t i;
  unsigned int seed;
  mod_t maxmod, curmod;
  tid_t myid;
  arg_t * arg;
  clustering_t * maxcluster, * curcluster;
  const graph_t * graph;
  objective_t * objective;

  arg = ptr;

  objective = arg->objective;
  graph = arg->graph;
  
  myid = dlthread_get_id(objective->comm);

  maxmod = 0;
  maxcluster = NULL;
  curcluster = NULL;

  seed = objective->seed;

  if (myid == 0) {
    if (objective->verbosity >= NERSTRAND_VERBOSITY_LOW) {
      print_objective(objective);
    }
  }

  for (i=0;i<objective->nruns;++i) {
    objective->runnum = i;

    curcluster = par_cluster_graph(objective,graph); 

    curmod = calc_total_modularity(curcluster,graph);

    if (myid == 0) {
      if (objective->modstats) {
        objective->modstats[i] = curmod;
      }
      objective->seed = ++seed;
    }
    if (maxcluster == NULL) {
      maxcluster = curcluster;
      maxmod = curmod;
    } else if (curmod > maxmod) {
      par_free_clustering(maxcluster,objective->comm);
      maxcluster = curcluster;
      maxmod = curmod;
    } else {
      par_free_clustering(curcluster,objective->comm);
    }
  }

  if (myid == 0) {
    /* copy the real clustering out */
    arg->cluster = maxcluster;
  }
}


/**
 * @brief Internal main clustering function for nerstrand.
 *
 * @param objective The objective containing parameters for the clustering.
 * @param graph The graph to cluster.
 *
 * @return The generated clustering.
 */
static clustering_t * __nerstrand_cluster(
    objective_t * const objective, 
    graph_t const * const graph)
{
  arg_t arg;
  tid_t const nthreads = objective->nthreads;
  clustering_t * finalcluster = NULL;

  objective->comm = DLTHREAD_COMM_ROOT;

  if (nthreads > 1) {
    arg.objective = objective;
    arg.graph = graph;

    dlthread_launch(nthreads,&__thread_cluster,&arg);

    finalcluster = arg.cluster;
  } else {
    size_t i;
    mod_t maxmod = 0;
    mod_t curmod;
    clustering_t * curcluster = NULL;

    if (objective->verbosity >= NERSTRAND_VERBOSITY_LOW) {
      print_objective(objective);
    }

    for (i=0;i<objective->nruns;++i) {
      objective->runnum = i;

      curcluster = cluster_graph(objective,graph); 

      curmod = calc_total_modularity(curcluster,graph);

      if (objective->modstats) {
        objective->modstats[i] = curmod;
      }

      if (finalcluster == NULL) {
        finalcluster = curcluster;
        maxmod = curmod;
      } else if (curmod > maxmod) {
        par_free_clustering(finalcluster,objective->comm);
        finalcluster = curcluster;
        maxmod = curmod;
      } else {
        par_free_clustering(curcluster,objective->comm);
      }
    }
  }

  /* output results */
  if (objective->modstats) {
    print_modularity_stats(objective->modstats,objective->nruns);
  }
  if (objective->icstats) {
    print_initialclustering_stats(objective->icstats);
  }
  if (objective->lbstats) {
    print_loadbalance_stats(objective->lbstats);
  }
  if (objective->aggstats) {
    print_aggregation_stats(objective->aggstats);
  }
  if (objective->refstats) {
    print_refinement_stats(objective->refstats);
  }

  return finalcluster;
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


double * nerstrand_init_options(void)
{
  double * options = double_init_alloc(NERSTRAND_VAL_OFF,NERSTRAND_NOPTIONS);
  return options; 
}


clustering_t * nerstrand_internal_cluster(
    objective_t * objective, 
    const graph_t * graph);


clustering_t * nerstrand_internal_cluster(
    objective_t * const objective, 
    const graph_t * const graph)
{
  return __nerstrand_cluster(objective,graph);
}


int nerstrand_cluster_explicit(
    const vtx_t * const r_nvtxs, 
    const adj_t * const xadj, 
    const vtx_t * const adjncy, 
    const wgt_t * const adjwgt, 
    const double * const options, 
    cid_t * const r_nclusters, 
    cid_t * const cid, 
    double * const r_mod)
{
  graph_t * graph;
  objective_t * objective;
  clustering_t * clustering;

  if (parse_objective(options,&objective) != NERSTRAND_SUCCESS) {
    return NERSTRAND_ERROR_INVALIDOPTIONS;
  }

  graph = distribute_graph(*r_nvtxs,xadj,adjncy,adjwgt,NULL,NULL,
      objective->nthreads,objective->distribution,objective->blocksize);
   
  clustering = __nerstrand_cluster(objective,graph);

  if (cid != NULL) {
    unify_where((const cid_t * const *)clustering->where,graph,cid);
  }
  if (r_mod != NULL) {
    *r_mod = (double)calc_total_modularity(clustering,graph);
  }
  if (r_nclusters != NULL) {
    *r_nclusters = clustering->nclusters;
  }

  free_clustering(clustering);
  free_graph(graph);
  free_objective(objective);

  return 1;
}


int nerstrand_cluster_kway(
    const vtx_t * const r_nvtxs, 
    const adj_t * const xadj, 
    const vtx_t * const adjncy, 
    const wgt_t * const adjwgt, 
    const cid_t * const r_nclusters, 
    cid_t * const cid, 
    double * const r_mod)
{
  graph_t * graph;
  objective_t * objective;
  clustering_t * clustering;

  objective = create_objective();
  objective->nclusters = *r_nclusters;
  setup_objective(objective);

  graph = distribute_graph(*r_nvtxs,xadj,adjncy,adjwgt,NULL,NULL,
      objective->nthreads,objective->distribution,objective->blocksize);
   
  clustering = __nerstrand_cluster(objective,graph);

  if (cid != NULL) {
    unify_where((const cid_t * const *)clustering->where,graph,cid);
  }
  if (r_mod != NULL) {
    *r_mod = (double)calc_total_modularity(clustering,graph);
  }

  free_clustering(clustering);
  free_graph(graph);
  free_objective(objective);

  return 1;
}


int nerstrand_cluster_anyway(
    const vtx_t * const r_nvtxs, 
    const adj_t * const xadj, 
    const vtx_t * const adjncy, 
    const wgt_t * const adjwgt, 
    cid_t * const r_nclusters, 
    cid_t * const cid, 
    double * const r_mod)
{
  graph_t * graph;
  objective_t * objective;
  clustering_t * clustering;

  objective = create_objective();
  objective->nclusters = 0;
  objective->parttype = NERSTRAND_PARTITION_ANYWAY;
  setup_objective(objective);

  graph = distribute_graph(*r_nvtxs,xadj,adjncy,adjwgt,NULL,NULL,
      objective->nthreads,objective->distribution,objective->blocksize);
   
  clustering = __nerstrand_cluster(objective,graph);

  if (cid != NULL) {
    unify_where((const cid_t * const *)clustering->where,graph,cid);
  }
  if (r_nclusters != NULL) {
    *r_nclusters = clustering->nclusters;
  }
  if (r_mod != NULL) {
    *r_mod = (double)calc_total_modularity(clustering,graph);
  }

  free_clustering(clustering);
  free_graph(graph);
  free_objective(objective);

  return 1;
}




#endif
