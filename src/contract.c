/**
 * @file contract.c
 * @brief Functions for contracting a graph
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2014-03-03
 */




#ifndef NERSTRAND_CONTRACT_C
#define NERSTRAND_CONTRACT_C




#include "contract.h"




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


#define DEF_MASK_SIZE (0x1000)
static vtx_t const MASK_SIZE = DEF_MASK_SIZE;




/******************************************************************************
* PRIVATE UTILITY FUNCTIONS ***************************************************
******************************************************************************/


/**
 * @brief Reverse the bits of a key for a given mask. Use for hashing into a
 * secondary hash table to prevent collisions.
 *
 * @param n The key to reverse.
 * @param mask The mask of the bits to be reversed.
 *
 * @return The reversed key.
 */
static inline vtx_t __reverse(
    vtx_t const n, 
    vtx_t const mask)
{
  vtx_t r = vtx_reversebits(n);
  size_t const mb = vtx_downlog2(mask);
  size_t const vs = sizeof(vtx_t)*8;
  if (vs >= 2*mb) {
    return (r >> (vs - (2*mb))) & mask;
  } else {
    return r >> (vs - mb);
  }
}


/**
 * @brief Internal function for collapsing edges together using a double-hash
 * table.
 *
 * @param mask The bit mask to use for hashing.
 * @param c The local coarse vertex edges are being collapsed for.
 * @param cg The global coarse vertex edges are being collapsed for.
 * @param k The coarse neighbor of the edge.
 * @param ewgt The weight of the edge.
 * @param cadjncy The coarse adjacency list.
 * @param cadjwgt The coarse edge weight.
 * @param ciadjwgt The coarse vertex-interal edge weight.
 * @param htable1 The primary hash table.
 * @param htable2 The secondary hash table (for reverse hashing).
 * @param cnedges The number of coarse edges genrated so far.
 *
 * @return The new number of coarse edges.
 */
static inline adj_t __handle_edge_masked(
    vtx_t const mask,
    vtx_t const c, 
    vtx_t const cg,
    vtx_t const k,
    wgt_t const ewgt,
    vtx_t * cadjncy, 
    wgt_t * cadjwgt,
    wgt_t * ciadjwgt,
    wgt_t * ceadjwgt,
    adj_t * htable1,
    adj_t * htable2, 
    adj_t cnedges)
{
  vtx_t l;
  adj_t i;

  if (k == c || k == cg) {
    /* internal edge */
    ciadjwgt[c] += ewgt;
    ceadjwgt[c] -= ewgt;
  } else {
    /* external edge */
    l = k&mask;
    i = htable1[l];
    if (i == NULL_ADJ) {
      /* new edge */
      cadjncy[cnedges] = k;
      cadjwgt[cnedges] = ewgt;
      htable1[l] = cnedges++; 
    } else if (k == cadjncy[i]) {
      /* duplicate edge */
      cadjwgt[i] += ewgt;
    } else {
      l = __reverse(l,mask);
      for (;;l=((l+1)&mask)) {
        i = htable2[l];
        if (i == NULL_ADJ) {
          /* new edge */
          cadjncy[cnedges] = k;
          cadjwgt[cnedges] = ewgt;
          htable2[l] = cnedges++; 
          break;
        } else if (k == cadjncy[i]) {
          /* duplicate edge */
          cadjwgt[i] += ewgt;
          break;
        } 
      }
    }
  }
  return cnedges;
}




/******************************************************************************
* PRIVATE SERIAL FUNCTIONS ****************************************************
******************************************************************************/


/**
 * @brief Contract a graph using summation of weights
 *
 * @param mgraph The graph to contract
 * @param aggregate The aggregation to contract based on
 *
 * @return The contracted (coarser) graph 
 */
static graph_t * __contract_SUM(
    mgraph_t * const mgraph, 
    aggregate_t const * const aggregate)
{
  adj_t cnedges,l,j;
  tid_t myid,o,t;
  vtx_t v,c,cg,mask,k,omask;
  adj_t * htable1, * htable2;

  adj_t * mycxadj;
  vtx_t * mycadjncy; 
  wgt_t * mycadjwgt;
  wgt_t * myciadjwgt;
  wgt_t * myceadjwgt;

  graphdist_t dist;

  /* make accessing my old graph easy */
  graph_t const * const graph = mgraph->graph;
  tid_t const nthreads = graph->npar;
  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy;
  wgt_t const * const * const adjwgt = (wgt_t const * const *)graph->adjwgt;
  wgt_t const * const * const iadjwgt = (wgt_t const * const *)graph->iadjwgt;
  wgt_t const * const * const eadjwgt = (wgt_t const * const *)graph->eadjwgt;

  calc_graph_dist(vtx_max_value(aggregate->nvtxs,nthreads),&dist);

  adjust_cmap(aggregate->cmap,graph->mynvtxs,nthreads,graph->dist,
      dist);

  /* easy access to aggregation */
  vtx_t * const cnvtxs = aggregate->nvtxs;
  vtx_t const * const * const cmap = (vtx_t const * const *)aggregate->cmap;
  vtx_t const * const * const match = (vtx_t const * const *)aggregate->match;
  vtx_t const * const * const fmap = (vtx_t const * const *)aggregate->fmap;

  /* create me new graph and give me easy access */
  adj_t ** const cxadj = r_adj_alloc(nthreads);
  vtx_t ** const cadjncy = r_vtx_alloc(nthreads);
  wgt_t ** const cadjwgt = r_wgt_alloc(nthreads);

  adj_t * const maxdeg = adj_alloc(nthreads);

  for (myid=0;myid<nthreads;++myid) {
    cxadj[myid] = adj_alloc(cnvtxs[myid]+1);
    /* count possible edges */
    cnedges = 0;
    maxdeg[myid] = 0;
    for (c=0;c<cnvtxs[myid];++c) {
      v = fmap[myid][c];
      o = myid;
      l= 0;
      do {
        l += xadj[o][v+1] - xadj[o][v];
        v = match[o][v];
        if (v >= graph->mynvtxs[o]) {
          o = gvtx_to_tid(v,graph->dist);
          v = gvtx_to_lvtx(v,graph->dist);
        }
      } while (!(o == myid && v == fmap[myid][c]));
      dl_storemax(maxdeg[myid],l);
      cnedges += l;
    }

    cadjncy[myid] = vtx_alloc(cnedges);
    cadjwgt[myid] = wgt_alloc(cnedges);
  }

  wgt_t ** const ciadjwgt = r_wgt_dalloc(cnvtxs,sizeof(vtx_t),nthreads);
  wgt_t ** const ceadjwgt = r_wgt_dalloc(cnvtxs,sizeof(vtx_t),nthreads);

  htable1 = NULL;
  htable2 = NULL;
  omask = 0;

  for (myid=0;myid<nthreads;++myid) {
    mycxadj = cxadj[myid];
    mycadjncy = cadjncy[myid];
    mycadjwgt = cadjwgt[myid];
    myciadjwgt = ciadjwgt[myid];
    myceadjwgt = ceadjwgt[myid];

    /* set up edge hash table */
    mask = MASK_SIZE;

    /* this could be problematic regarding max degree */
    while (maxdeg[myid] > (mask >> 3) && max_gvtx(graph) > mask) {
      dprintf("Bumping up mask size from "PF_VTX_T" to "PF_VTX_T"\n",mask,
          (mask << 1));
      mask <<=1;
    }
    if (mask > omask) {
      if (htable1) {
        dl_free(htable1);
        dl_free(htable2);
      }
      htable1 = adj_init_alloc(NULL_ADJ,mask);
      htable2 = adj_init_alloc(NULL_ADJ,mask);
    }
    omask = mask;
    --mask;

    dprintf("Using a mask of "PF_VTX_T" (%lX) for "PF_VTX_T" coarse "
        "vertices\n",mask,(long unsigned int)mask,cnvtxs[myid]);

    cnedges = 0;
    mycxadj[0] = 0;

    /* set max degree for the coarse graph */
    maxdeg[myid] = 0;
    for (c=0;c<cnvtxs[myid];++c) {
      cg = lvtx_to_gvtx(c,myid,dist);
      /* initialize the coarse vertex */
      myciadjwgt[c] = 0;

      v = fmap[myid][c];
      o = myid;
      DL_ASSERT_EQUALS(myid,gvtx_to_tid(lvtx_to_gvtx(v,myid,graph->dist),
          graph->dist), PF_TID_T);
      myceadjwgt[c] = 0;
      do { 
        DL_ASSERT_EQUALS(c,gvtx_to_lvtx(cmap[o][v],dist),PF_VTX_T);

        /* transfer over vertex stuff from v and u */
        if (iadjwgt) {
          myciadjwgt[c] += iadjwgt[o][v];
        }
        myceadjwgt[c] += eadjwgt[o][v];

        if (adjwgt) {
          for (j=xadj[o][v];j<xadj[o][v+1];++j) {
            k = adjncy[o][j];
            if (k < graph->mynvtxs[o]) {
              t = o;
            } else {
              t = gvtx_to_tid(k,graph->dist);
              k = gvtx_to_lvtx(k,graph->dist);
            }
            k = cmap[t][k];
            if (gvtx_to_tid(k,dist) == myid) {
              k = gvtx_to_lvtx(k,dist);
            }
            cnedges = __handle_edge_masked(mask,c,cg,k,adjwgt[o][j], \
                mycadjncy,mycadjwgt,myciadjwgt,myceadjwgt,htable1,htable2, \
                cnedges);
          }
        } else {
          for (j=xadj[o][v];j<xadj[o][v+1];++j) {
            k = adjncy[o][j];
            if (k < graph->mynvtxs[o]) {
              t = o;
            } else {
              t = gvtx_to_tid(k,graph->dist);
              k = gvtx_to_lvtx(k,graph->dist);
            }
            k = cmap[t][k];
            if (gvtx_to_tid(k,dist) == myid) {
              k = gvtx_to_lvtx(k,dist);
            }
            cnedges = __handle_edge_masked(mask,c,cg,k,1.0,mycadjncy, \
                mycadjwgt,myciadjwgt,myceadjwgt,htable1,htable2,cnedges);
          }
        }

        v = match[o][v];
        if (v >= graph->mynvtxs[o]) {
          o = gvtx_to_tid(v,graph->dist);
          v = gvtx_to_lvtx(v,graph->dist);
        }
      } while (!(myid == o && v == fmap[myid][c]));

      /* clear the htable */
      for (j = cnedges;j > mycxadj[c];) {
        --j;
        l = (mycadjncy[j]&mask);
        htable1[l] = NULL_ADJ;
        l = __reverse(l,mask);
        for (;htable2[l] != NULL_ADJ;l=((l+1)&mask)) {
          htable2[l] = NULL_ADJ;
        }
      }

      mycxadj[c+1] = cnedges;
      dl_storemax(maxdeg[myid],mycxadj[c+1] - mycxadj[c]);
    }
  }
  dl_free(htable1);
  dl_free(htable2);


  graph_t * cgraph = setup_graph(cnvtxs,cxadj,cadjncy,cadjwgt,ciadjwgt,
      ceadjwgt,maxdeg,nthreads);

  cgraph->dist = dist;

  cgraph->free_xadj = 1;
  cgraph->free_adjncy = 1;
  cgraph->free_adjwgt = 1;
  cgraph->free_iadjwgt = 1;

  /* phantom edges are contagious */
  cgraph->phantom_edges = graph->phantom_edges;

  DL_ASSERT(check_graph(cgraph) == 1, "Bad graph generated in coarsening\n");

  return cgraph;
}




/******************************************************************************
* PRIVATE PARALLEL FUNCTIONS **************************************************
******************************************************************************/


/**
 * @brief Contract a graph using summation of weights
 *
 * @param mgraph The graph to contract
 * @param aggregate The aggegration to contract based on
 *
 * @return The contracted (coarser) graph
 */
static graph_t * __par_contract_SUM(
    mgraph_t * const mgraph, 
    aggregate_t const * const aggregate,
    dlthread_comm_t const comm)
{
  adj_t cnedges,l,maxdeg,j;
  tid_t o,t;
  vtx_t v,c,cg,mask,k;
  adj_t * htable1, * htable2;
  graphdist_t dist;

  tid_t const myid = dlthread_get_id(comm);

  /* make accessing my old graph easy */
  graph_t const * const graph = mgraph->graph;
  tid_t const nthreads = graph->npar;
  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy;
  wgt_t const * const * const adjwgt = (wgt_t const * const *)graph->adjwgt;
  wgt_t const * const * const iadjwgt = (wgt_t const * const *)graph->iadjwgt;
  wgt_t const * const * const eadjwgt = (wgt_t const * const *)graph->eadjwgt;

  calc_graph_dist(vtx_max_value(aggregate->nvtxs,nthreads),&dist);

  par_adjust_cmap(aggregate->cmap[myid],mynvtxs,graph->dist,dist);

  /* easy access to aggregation */
  vtx_t const mycnvtxs = aggregate->nvtxs[myid];
  vtx_t const * const * const cmap = (vtx_t const * const *)aggregate->cmap;
  vtx_t const * const * const match = (vtx_t const * const *)aggregate->match;
  vtx_t const * const fmap = aggregate->fmap[myid];

  /* count possible edges */
  cnedges = 0;
  maxdeg = 0;
  for (c=0;c<mycnvtxs;++c) {
    v = fmap[c];
    o = myid;
    l = 0;
    do {
      l += xadj[o][v+1] - xadj[o][v];
      v = match[o][v];
      if (v >= graph->mynvtxs[o]) {
        o = gvtx_to_tid(v,graph->dist);
        v = gvtx_to_lvtx(v,graph->dist);
      }
    } while (!(o == myid && v == fmap[c]));
    dl_storemax(maxdeg,l);
    cnedges += l;
  }

  adj_t * const mycxadj = adj_alloc(mycnvtxs+1);
  vtx_t * const mycadjncy = vtx_alloc(cnedges);
  wgt_t * const mycadjwgt = wgt_alloc(cnedges);
  wgt_t * const myciadjwgt = wgt_alloc(mycnvtxs);
  wgt_t * const myceadjwgt = wgt_alloc(mycnvtxs);

  htable1 = NULL;
  htable2 = NULL;

  /* set up edge hash table */
  mask = MASK_SIZE;

  /* this could be problematic regarding max degree */
  while (maxdeg > (mask >> 3) && max_gvtx(graph) > mask) {
    dprintf("["PF_TID_T"] Bumping up mask size from "PF_VTX_T" to "PF_VTX_T
        "\n",myid,mask,(mask << 1));
    mask <<=1;
  }
  htable1 = adj_init_alloc(NULL_ADJ,mask);
  htable2 = adj_init_alloc(NULL_ADJ,mask);

  --mask;

  dprintf("["PF_TID_T"] Using a mask of "PF_VTX_T" (%lX) for "PF_VTX_T
      " coarse vertices\n",myid,mask,(long unsigned int)mask,mycnvtxs);

  cnedges = 0;
  mycxadj[0] = 0;

  dlthread_barrier(comm);

  /* set max degree for the coarse graph */
  maxdeg = 0;
  for (c=0;c<mycnvtxs;++c) {
    cg = lvtx_to_gvtx(c,myid,dist);
    /* initialize the coarse vertex */
    myciadjwgt[c] = 0;

    v = fmap[c];
    o = myid;
    DL_ASSERT_EQUALS(myid,gvtx_to_tid(lvtx_to_gvtx(v,myid,graph->dist),
        graph->dist), PF_TID_T);
    myceadjwgt[c] = 0;
    do { 
      DL_ASSERT_EQUALS(c,gvtx_to_lvtx(cmap[o][v],dist),PF_VTX_T);

      /* transfer over vertex stuff from v and u */
      if (iadjwgt) {
        myciadjwgt[c] += iadjwgt[o][v];
      }
      myceadjwgt[c] += eadjwgt[o][v];

      if (adjwgt) {
        for (j=xadj[o][v];j<xadj[o][v+1];++j) {
          k = adjncy[o][j];
          if (k < graph->mynvtxs[o]) {
            t = o;
          } else {
            t = gvtx_to_tid(k,graph->dist);
            k = gvtx_to_lvtx(k,graph->dist);
          }
          k = cmap[t][k];
          if (gvtx_to_tid(k,dist) == myid) {
            k = gvtx_to_lvtx(k,dist);
          }
          cnedges = __handle_edge_masked(mask,c,cg,k,adjwgt[o][j], \
              mycadjncy,mycadjwgt,myciadjwgt,myceadjwgt,htable1,htable2, \
              cnedges);
        }
      } else {
        for (j=xadj[o][v];j<xadj[o][v+1];++j) {
          k = adjncy[o][j];
          if (k < graph->mynvtxs[o]) {
            t = o;
          } else {
            t = gvtx_to_tid(k,graph->dist);
            k = gvtx_to_lvtx(k,graph->dist);
          }
          k = cmap[t][k];
          if (gvtx_to_tid(k,dist) == myid) {
            k = gvtx_to_lvtx(k,dist);
          }
          cnedges = __handle_edge_masked(mask,c,cg,k,1.0,mycadjncy, \
              mycadjwgt,myciadjwgt,myceadjwgt,htable1,htable2,cnedges);
        }
      }

      v = match[o][v];
      if (v >= graph->mynvtxs[o]) {
        o = gvtx_to_tid(v,graph->dist);
        v = gvtx_to_lvtx(v,graph->dist);
      }
    } while (!(myid == o && v == fmap[c]));

    /* clear the htable */
    for (j = cnedges;j > mycxadj[c];) {
      --j;
      l = (mycadjncy[j]&mask);
      htable1[l] = NULL_ADJ;
      l = __reverse(l,mask);
      for (;htable2[l] != NULL_ADJ;l=((l+1)&mask)) {
        htable2[l] = NULL_ADJ;
      }
    }

    mycxadj[c+1] = cnedges;
    dl_storemax(maxdeg,mycxadj[c+1] - mycxadj[c]);
  }

  dl_free(htable1);
  dl_free(htable2);

  graph_t * cgraph = par_setup_graph(mycnvtxs,mycxadj,mycadjncy,mycadjwgt,
      myciadjwgt,myceadjwgt,maxdeg,comm);

  if (myid == 0) {
    cgraph->dist = dist;
    cgraph->free_xadj = 1;
    cgraph->free_adjncy = 1;
    cgraph->free_adjwgt = 1;
    cgraph->free_iadjwgt = 1;

    /* phantom edges are contagious */
    cgraph->phantom_edges = graph->phantom_edges;
  }
  dlthread_barrier(comm);

  DL_ASSERT(check_graph(cgraph) == 1, "Bad graph generated in coarsening\n");

  return cgraph;
}




/******************************************************************************
* PUBLIC SERIAL FUNCTIONS *****************************************************
******************************************************************************/


graph_t * contract_graph(
    objective_t * const objective, 
    mgraph_t * const mgraph, 
    aggregate_t const * const agg)
{
  dl_start_timer(&(objective->timers.contraction));

  graph_t * cgraph;

  switch (objective->contype) {
    case NERSTRAND_CONTRACT_SUM:
      cgraph = __contract_SUM(mgraph,agg);
      break;
    default:
      dl_error("Unknown contraction type %d\n",objective->contype);
  }

  dl_stop_timer(&(objective->timers.contraction));

  return cgraph;
}




/******************************************************************************
* PUBLIC PARALLEL FUNCTIONS ***************************************************
******************************************************************************/


graph_t * par_contract_graph(
    objective_t * const objective, 
    mgraph_t * const mgraph, 
    aggregate_t const * const agg)
{
  graph_t * cgraph = NULL;

  tid_t const myid = dlthread_get_id(objective->comm);

  if (myid == 0) {
    dl_start_timer(&(objective->timers.contraction));
  }

  switch (objective->contype) {
    case NERSTRAND_CONTRACT_SUM:
      cgraph = __par_contract_SUM(mgraph,agg,objective->comm);
      break;
    default:
      dl_error("Unknown contraction type %d\n",objective->contype);
  }

  if (myid == 0) {
    dl_stop_timer(&(objective->timers.contraction));
  }

  return cgraph;
}




#endif
