/**
 * @file mgraph.c
 * @brief Functions for manipulating and allocating multi-level graphs
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-10-08
 */




#ifndef NERSTRAND_MGRAPH_C
#define NERSTRAND_MGRAPH_C




#include "mgraph.h"




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLMEM_PREFIX mgraph
#define DLMEM_TYPE_T mgraph_t
#define DLMEM_STATIC
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T




/******************************************************************************
* PUBLIC SERIAL FUNCTIONS *****************************************************
******************************************************************************/


mgraph_t * init_mgraph(
    mgraph_t * const mgraph) 
{
  mgraph->level = 0;
  mgraph->graph = NULL;
  mgraph->cmap = NULL;
  mgraph->coarser = NULL;
  mgraph->finer = NULL;
  mgraph->free_cmap = 0;
  mgraph->free_graph = 0;
  mgraph->free_ucinfo = 0;

  return mgraph;
}


mgraph_t * setup_mgraph(
    size_t level, 
    graph_t const * const graph, 
    vtx_t ** cmap,
    mgraph_t * coarser,
    mgraph_t * finer)
{
  mgraph_t * mgraph;
  
  mgraph = mgraph_calloc(1);

  mgraph->level = level;
  mgraph->graph = graph;
  mgraph->cmap = cmap;
  mgraph->coarser = coarser;
  mgraph->finer = finer;

  return mgraph;
}


void free_mgraph(
    mgraph_t * mgraph)
{
  tid_t const nthreads = mgraph->graph->npar;

  if (mgraph->free_graph) {
    free_graph((graph_t*)mgraph->graph);
  }
  if (mgraph->free_cmap) {
    r_vtx_free(mgraph->cmap,nthreads);
  }
  if (mgraph->free_ucinfo) {
    free_ucinfos(mgraph->ucinfo,nthreads);
  }
  dl_free(mgraph);
}


void adjust_cmap(
    vtx_t * const * const cmap,
    vtx_t const * const nvtxs, 
    tid_t const nthreads,
    graphdist_t const olddist,
    graphdist_t const newdist)
{
  tid_t myid;
  vtx_t i,k;
  tid_t o;

  for (myid=0;myid<nthreads;++myid) {
    for (i=0;i<nvtxs[myid];++i) {
      if (cmap[myid][i] >= olddist.offset) { /* remote vertex */
        k = gvtx_to_lvtx(cmap[myid][i],olddist);
        o = gvtx_to_tid(cmap[myid][i],olddist);
        cmap[myid][i] = lvtx_to_gvtx(k,o,newdist);
      }
    }
  }
}




/******************************************************************************
* PUBLIC PARALLEL FUNCTIONS ***************************************************
******************************************************************************/


mgraph_t * par_setup_mgraph(
    size_t level, 
    graph_t const * const graph, 
    vtx_t * const cmap, 
    mgraph_t * const coarser, 
    mgraph_t * const finer,
    dlthread_comm_t const comm)
{
  mgraph_t * mgraph;
  mgraph_t ** r_mgraph;

  tid_t const myid = dlthread_get_id(comm);
  tid_t const nthreads = dlthread_get_nthreads(comm);

  DL_ASSERT_EQUALS(nthreads,graph->npar,PF_TID_T);

  r_mgraph = dlthread_get_shmem(sizeof(mgraph),comm);

  if (myid == 0) {
    mgraph = mgraph_calloc(1);
    mgraph->level = level;
    mgraph->graph = graph;
    mgraph->coarser = coarser;
    mgraph->finer = finer;
    if (cmap) {
      mgraph->cmap = r_vtx_alloc(nthreads);
    } else {
      mgraph->cmap = NULL;
    }
    *r_mgraph = mgraph;
  }
  dlthread_barrier(comm);

  mgraph = *r_mgraph;
  dlthread_free_shmem(r_mgraph,comm);

  if (cmap) {
    mgraph->cmap[myid] = cmap;
    dlthread_barrier(comm);
  }

  return mgraph;
}


void par_free_mgraph(
    mgraph_t * mgraph,
    dlthread_comm_t const comm)
{
  vtx_t const myid = dlthread_get_id(comm);

  if (mgraph->free_graph) {
    par_free_graph((graph_t*)mgraph->graph,comm);
  }

  if (mgraph->free_cmap) {
    dl_free(mgraph->cmap[myid]);

    dlthread_barrier(comm);
    if (myid == 0) {
      dl_free(mgraph->cmap);
    }
  }
  if (mgraph->free_ucinfo) {
    free_ucinfo(mgraph->ucinfo[myid]);

    dlthread_barrier(comm);
    if (myid == 0) {
      dl_free(mgraph->ucinfo);
    }
  }

  dlthread_barrier(comm);
  if (myid == 0) {
    dl_free(mgraph);
  }
}


void par_adjust_cmap(
    vtx_t * const cmap, 
    vtx_t const mynvtxs, 
    graphdist_t const olddist, 
    graphdist_t const newdist)
{
  vtx_t i,k;
  tid_t o;
  for (i=0;i<mynvtxs;++i) {
    if (cmap[i] >= olddist.offset) { /* remote vertex */
      k = gvtx_to_lvtx(cmap[i],olddist);
      o = gvtx_to_tid(cmap[i],olddist);
      cmap[i] = lvtx_to_gvtx(k,o,newdist);
    }
  }
}




#endif
