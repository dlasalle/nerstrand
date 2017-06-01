/**
 * @file coarsen.c
 * @brief Functions for coarsening a graph
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2014-03-03
 */




#ifndef NERSTRAND_COARSEN_C
#define NERSTRAND_COARSEN_C




#include "coarsen.h"
#include "aggregate.h"
#include "contract.h"
#include "sparsen.h"




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


mgraph_t * coarsen_graph(
    objective_t * const objective,
    mgraph_t * const mgraph)
{
  graph_t * cgraph;
  aggregate_t * agg;

  dl_start_timer(&(objective->timers.coarsening));

  /* aggregate */
  agg = aggregate_graph(objective,mgraph);  

  mgraph->cmap = agg->cmap;
  mgraph->free_cmap = 1;

  vprintf(objective->verbosity,NERSTRAND_VERBOSITY_MAXIMUM,
      "Aggregation generated "PF_VTX_T" coarse vertices\n",
      vtx_sum(agg->nvtxs,agg->npar));

  /* contract */
  cgraph = contract_graph(objective,mgraph,agg);

  vprintf(objective->verbosity,NERSTRAND_VERBOSITY_MAXIMUM,
      "Contracted graph from "PF_VTX_T" vertices to "PF_VTX_T" vertices\n",
      mgraph->graph->nvtxs,cgraph->nvtxs);

  free_aggregate(agg);

  sparsen_graph(objective,cgraph,mgraph->graph->nedges/2);

  mgraph->coarser = setup_mgraph(mgraph->level+1,cgraph,NULL,NULL,mgraph);
  mgraph->coarser->free_graph = 1;

  dl_stop_timer(&(objective->timers.coarsening));

  return mgraph->coarser;
}




/******************************************************************************
* PUBLIC PARALLEL FUNCTIONS ***************************************************
******************************************************************************/


mgraph_t * par_coarsen_graph(
    objective_t * const objective, 
    mgraph_t * const mgraph)
{
  graph_t * cgraph;
  aggregate_t * agg;

  tid_t const myid = dlthread_get_id(objective->comm);

  if (myid == 0) {
    dl_start_timer(&(objective->timers.coarsening));
  }

  /* aggregate */
  agg = par_aggregate_graph(objective,mgraph);  

  if (myid == 0) {
    mgraph->cmap = agg->cmap;
    mgraph->free_cmap = 1;
  }
  dlthread_barrier(objective->comm);

  if (myid == 0) {
    vprintf(objective->verbosity,NERSTRAND_VERBOSITY_MAXIMUM,
        "Aggregation generated "PF_VTX_T" coarse vertices\n",
        vtx_sum(agg->nvtxs,agg->npar));
  }

  /* contract */
  cgraph = par_contract_graph(objective,mgraph,agg);

  if (myid == 0) {
    vprintf(objective->verbosity,NERSTRAND_VERBOSITY_MAXIMUM,
        "Contracted graph from "PF_VTX_T" vertices to "PF_VTX_T" vertices\n",
        mgraph->graph->nvtxs,cgraph->nvtxs);
  }

  par_free_aggregate(agg);

  DL_ASSERT(mgraph->cmap != NULL,"Free'd cmap before too early\n");

  mgraph->coarser = par_setup_mgraph(mgraph->level+1,cgraph,NULL,NULL,
      mgraph,objective->comm);

  if (myid == 0) {
    mgraph->coarser->free_graph = 1;
    dl_stop_timer(&(objective->timers.coarsening));
  }

  return mgraph->coarser;
}


#endif
