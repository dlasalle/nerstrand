/**
 * @file graph.c
 * @brief Functions for graph operations
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-10-08
 */




#ifndef NERSTRAND_GRAPH_C
#define NERSTRAND_GRAPH_C




#include "graph.h"




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLMEM_PREFIX graph
#define DLMEM_TYPE_T graph_t
#define DLMEM_DLTYPE DLTYPE_STRUCT
#define DLMEM_INITFUNCTION init_graph
#define DLMEM_STATIC
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_INITFUNCTION
#undef DLMEM_DLTYPE
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX




/******************************************************************************
* INLINE FUNCTIONS ************************************************************
******************************************************************************/


/**
 * @brief Determine threshold for a vertex to be a "supernode".
 *
 * @param graph The graph for calibrate the threshold for.
 *
 * @return The threshold.
 */
static inline wgt_t __calc_snthresh(
    graph_t const * const graph)
{
  return graph->gadjwgt / sqrt(graph->nvtxs);
}




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


/**
 * @brief Distribute the vertices of a graph linearly (in continigous blocks). 
 *
 * @param nvtxs The number of vertices.
 * @param xadj The adjacecny list pointer.
 * @param nthreads The number of threads to distribute the vertices over.
 * @param mynvtxs The number of vertices owned by each thread.
 * @param mynedges The number of edges owned by each thread.
 * @param lvtx The local vertex numbers.
 * @param owner The array assigning vertices to threads.
 */
static void __distribute_block(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    tid_t const nthreads,
    vtx_t * const mynvtxs,
    adj_t * const mynedges,
    vtx_t * const lvtx,
    tid_t * const owner)
{
  vtx_t v;
  tid_t myid;
  adj_t j, nedgesleft;

  adj_t const nedges = xadj[nvtxs];
  adj_t const avgnedges = nedges / nthreads;

  /* distribute vertices based on edge count */
  nedgesleft = nedges;
  myid = 0;
  for (v=0;v<nvtxs;++v) {
    if ((nthreads-myid-1)*avgnedges > nedgesleft) {
      ++myid;
    }
    owner[v] = myid;
    j = xadj[v+1] - xadj[v]; 
    mynedges[myid] += j;
    nedgesleft -= j;
  }
  for (v=0;v<nvtxs;++v) {
    myid = owner[v];
    lvtx[v] = mynvtxs[myid]++;
  }
}


/**
 * @brief Distribute the vertices of a graph cyclicly. 
 *
 * @param nvtxs The number of vertices.
 * @param xadj The adjacecny list pointer.
 * @param nthreads The number of threads to distribute the vertices over.
 * @param mynvtxs The number of vertices owned by each thread.
 * @param mynedges The number of edges owned by each thread.
 * @param lvtx The local vertex numbers.
 * @param owner The array assigning vertices to threads.
 */
static void __distribute_cyclic(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    tid_t const nthreads,
    vtx_t * const mynvtxs,
    adj_t * const mynedges,
    vtx_t * const lvtx,
    tid_t * const owner)
{
  vtx_t i, v;
  tid_t myid;
  adj_t j, nedgesleft;
  vtx_t * perm;

  adj_t const nedges = xadj[nvtxs];
  adj_t const avgnedges = nedges / nthreads;

  perm = vtx_alloc(nvtxs);

  vtx_cyclicperm(perm,nthreads,nvtxs);

  /* distribute vertices based on edge count */
  nedgesleft = nedges;
  myid = 0;
  for (i=0;i<nvtxs;++i) {
    if ((nthreads-myid-1)*avgnedges > nedgesleft) {
      ++myid;
    }
    v = perm[i];
    owner[v] = myid;
    j = xadj[v+1] - xadj[v]; 
    mynedges[myid] += j;
    nedgesleft -= j;
  }
  for (i=0;i<nvtxs;++i) {
    myid = owner[i];
    lvtx[i] = mynvtxs[myid]++;
  }

  dl_free(perm);
}


/**
 * @brief Distribute the vertices of a graph block-cyclicly. 
 *
 * @param nvtxs The number of vertices.
 * @param xadj The adjacecny list pointer.
 * @param nthreads The number of threads to distribute the vertices over.
 * @param mynvtxs The number of vertices owned by each thread.
 * @param mynedges The number of edges owned by each thread.
 * @param lvtx The local vertex numbers.
 * @param owner The array assigning vertices to threads.
 * @param block This size of a block.
 */
static void __distribute_blockcyclic(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    tid_t const nthreads,
    vtx_t * const mynvtxs,
    adj_t * const mynedges,
    vtx_t * const lvtx,
    tid_t * const owner,
    vtx_t block)
{
  vtx_t i, v;
  tid_t myid;
  adj_t j, nedgesleft;
  vtx_t * perm;

  adj_t const nedges = xadj[nvtxs];
  adj_t const avgnedges = nedges / nthreads;

  perm = vtx_alloc(nvtxs);

  /* create cyclic permutation */
  if (nthreads * block > nvtxs) {
    /* adjust the block if it will be imbalanced */
    block = dl_max(1,nvtxs / nthreads);
  }
  vtx_blockcyclicperm(perm,nthreads,block,nvtxs);

  /* distribute vertices based on edge count */
  nedgesleft = nedges;
  myid = 0;
  for (i=0;i<nvtxs;++i) {
    if ((nthreads-myid-1)*avgnedges > nedgesleft) {
      ++myid;
    }
    v = perm[i];
    owner[v] = myid;
    j = xadj[v+1] - xadj[v]; 
    mynedges[myid] += j;
    nedgesleft -= j;
  }
  for (i=0;i<nvtxs;++i) {
    myid = owner[i];
    lvtx[i] = mynvtxs[myid]++;
  }

  dl_free(perm);
}




/******************************************************************************
* PUBLIC SERIAL FUNCTIONS *****************************************************
******************************************************************************/


graph_t * init_graph(
    graph_t * graph)
{
  /* should be a valid graph with 0 vertices */
  graph->nvtxs = 0;
  graph->nedges = 0;
  graph->tadjwgt = 0.0;
  graph->gadjwgt = 0.0;
  graph->snthresh = (wgt_t)FLT_MAX;
  /* null out pointers */
  graph->xadj = NULL;
  graph->adjncy = NULL;
  graph->adjwgt = NULL;
  graph->iadjwgt = NULL;
  graph->eadjwgt = NULL;
  graph->alias = NULL;
  /* setup distribution for one thread */
  graph->dist.mask = 0x0;
  graph->dist.shift = 0;
  graph->dist.offset = 0;
  /* no need to free nulls */
  graph->free_xadj = 0;
  graph->free_adjncy = 0;
  graph->free_adjwgt = 0;
  graph->free_iadjwgt = 0;
  /* phantom edge weights */
  graph->phantom_edges = 0;
  return graph;
}


graph_t * setup_graph(
    vtx_t * const nvtxs, 
    adj_t ** const xadj, 
    vtx_t ** const adjncy, 
    wgt_t ** const adjwgt,
    wgt_t ** const iadjwgt, 
    wgt_t ** const eadjwgt, 
    adj_t * const maxdeg, 
    tid_t const nthreads)
{
  vtx_t mymd,i;
  tid_t myid;

  graph_t * const graph = graph_calloc(1); 

  DL_ASSERT(xadj != NULL,"Null xadj passed into setup_graph()");
  DL_ASSERT(adjncy != NULL,"Null xadj passed into setup_graph()");

  /* set nvtxs */
  graph->nvtxs = vtx_sum(nvtxs,nthreads);
  graph->mynvtxs = nvtxs;

  /* set nedges */
  graph->mynedges = adj_alloc(nthreads);
  graph->nedges = 0;
  for (myid=0;myid<nthreads;++myid) {
    graph->nedges += (graph->mynedges[myid] = xadj[myid][graph->mynvtxs[myid]]);
  }

  /* set thread stuff */
  graph->npar = nthreads;
  calc_graph_dist(vtx_max_value(graph->mynvtxs,nthreads),&(graph->dist));
  
  /* these have been asserted to not be NULL */
  graph->xadj = xadj;
  graph->adjncy = adjncy;

  /* these are either NULL or their properly initialized -- good either way */
  graph->adjwgt = adjwgt;
  graph->iadjwgt = iadjwgt;

  /* if it's NULL, need to derive it */
  if (eadjwgt) {
    graph->eadjwgt = eadjwgt;
  } else {
    /* calculate eadjwgt */
    graph->eadjwgt = r_wgt_dalloc(graph->mynvtxs,sizeof(vtx_t),nthreads);
    if (adjwgt) {
      /* sum adjwgts */
      for (myid=0;myid<nthreads;++myid) {
        for (i=0;i<graph->mynvtxs[myid];++i) {
          graph->eadjwgt[myid][i] = wgt_sum(graph->adjwgt[myid]+xadj[myid][i],
              xadj[myid][i+1]-xadj[myid][i]);
        }
      }
    } else {
      /* use vertex degree */
      for (myid=0;myid<nthreads;++myid) {
        for (i=0;i<graph->mynvtxs[myid];++i) {
          graph->eadjwgt[myid][i] = xadj[myid][i+1]-xadj[myid][i];
        }
      }
    }
  }

  /* handle max degree */
  if (maxdeg) {
    graph->maxdeg = maxdeg;
  } else {
    graph->maxdeg = adj_alloc(nthreads);
    for (myid=0;myid<nthreads;++myid) {
      mymd = 0;
      for (i=0;i<graph->mynvtxs[myid];++i) {
        dl_storemax(mymd,xadj[myid][i+1]-xadj[myid][i]);
      }
      graph->maxdeg[myid] = mymd;
    }
  }

  /* calculate tadjwgt */
  graph->tadjwgt = 0.0;
  for (myid=0;myid<nthreads;++myid) {
    graph->tadjwgt += wgt_fa_sum(graph->eadjwgt[myid],graph->mynvtxs[myid]);
  }

  /* calculate gadjwgt */
  graph->gadjwgt = 0;
  if (iadjwgt) {
    for (myid=0;myid<nthreads;++myid) {
      graph->gadjwgt += wgt_fa_sum(graph->iadjwgt[myid],graph->mynvtxs[myid]);
    }
  }
  graph->gadjwgt += graph->tadjwgt;

  /* set supernode info */
  graph->snthresh = __calc_snthresh(graph); 

  dprintf("Setup a graph with "PF_VTX_T" vertices, "PF_ADJ_T" edges, split "
      "among "PF_TID_T" threads\n",graph->nvtxs,graph->nedges,nthreads);

  return graph;
}


graph_t * distribute_graph(
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    wgt_t const * const adjwgt, 
    wgt_t const * const iadjwgt, 
    wgt_t const * const eadjwgt, 
    tid_t const nthreads, 
    int const distribution,
    vtx_t const block)
{
  vtx_t i,k,v,deg,mynvtxs;
  adj_t j,l;
  tid_t myid;
  vtx_t * dmynvtxs;
  adj_t * dmynedges;
  adj_t ** dxadj;
  vtx_t ** dadjncy, ** dalias;
  wgt_t ** dadjwgt = NULL, ** diadjwgt = NULL, ** deadjwgt = NULL;
  graphdist_t dist;
  graph_t * graph;

  DL_ASSERT(nvtxs>0, "setup_graph() called with nvtxs="PF_VTX_T"\n",
      nvtxs);

  vtx_t * lvtx = vtx_alloc(nvtxs);
  tid_t * owner = tid_alloc(nvtxs);

  /* set arrays */
  dmynvtxs = vtx_calloc(nthreads);
  dmynedges = adj_calloc(nthreads);
  dalias = r_vtx_alloc(nthreads);
  dxadj = r_adj_alloc(nthreads);
  dadjncy = r_vtx_alloc(nthreads);
  if (adjwgt) {
    dadjwgt = r_wgt_alloc(nthreads);
  }
  if (iadjwgt) {
    diadjwgt = r_wgt_alloc(nthreads);
  }

  /* I always want an eadjwgt */
  deadjwgt = r_wgt_alloc(nthreads);

  switch(distribution) {
    case NERSTRAND_DISTRIBUTION_BLOCK:
      __distribute_block(nvtxs,xadj,nthreads,dmynvtxs,dmynedges, \
          lvtx,owner);
      break;
    case NERSTRAND_DISTRIBUTION_CYCLIC:
      __distribute_cyclic(nvtxs,xadj,nthreads,dmynvtxs,dmynedges, \
          lvtx,owner);
      break;
    case NERSTRAND_DISTRIBUTION_BLOCKCYCLIC:
      __distribute_blockcyclic(nvtxs,xadj,nthreads,dmynvtxs, \
          dmynedges,lvtx,owner,block);
      break;
    default:
      dl_error("Unknown distribution '%d'\n",distribution);
  }

  calc_graph_dist(vtx_max_value(dmynvtxs,nthreads),&dist);

  /* allocate arrays */
  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = dmynvtxs[myid];
    dxadj[myid] = adj_alloc(mynvtxs+1);
    dxadj[myid][0] = 0;
    dalias[myid] = vtx_alloc(mynvtxs);
    dadjncy[myid] = vtx_alloc(dmynedges[myid]);
    deadjwgt[myid] = wgt_alloc(mynvtxs);
    if (iadjwgt) {
      diadjwgt[myid] = wgt_alloc(mynvtxs);
    }
    if (adjwgt) {
      dadjwgt[myid] = wgt_alloc(dmynedges[myid]);
    }
    /* zero counts for insertion later */
    dmynvtxs[myid] = 0;
    dmynedges[myid] = 0;
  }

  /* set xadj and iadjwgt */
  for (v=0;v<nvtxs;++v) { 
    myid = owner[v];
    i = dmynvtxs[myid]++; 
    dalias[myid][i] = v;
    DL_ASSERT_EQUALS(i,lvtx[v],PF_VTX_T);
    deg = xadj[v+1] - xadj[v];
    dxadj[myid][i+1] = dxadj[myid][i] + deg;
    l = dmynedges[myid];
    if (adjwgt) {
      deadjwgt[myid][i] = 0;
      for (j=xadj[v];j<xadj[v+1];++j) {
        k = adjncy[j];
        deadjwgt[myid][i] += adjwgt[j];
        if (owner[k] == myid) {
          dadjncy[myid][l] = lvtx[k];
        } else {
          dadjncy[myid][l] = lvtx_to_gvtx(lvtx[k],owner[k],dist);
        }
        dadjwgt[myid][l++] = adjwgt[j]; 
      }
    } else {
      deadjwgt[myid][i] = (wgt_t)deg;
      for (j=xadj[v];j<xadj[v+1];++j) {
        k = adjncy[j];
        if (owner[k] == myid) {
          dadjncy[myid][l] = lvtx[k];
        } else {
          dadjncy[myid][l] = lvtx_to_gvtx(lvtx[k],owner[k],dist);
        }
        ++l;
      }
    }
    DL_ASSERT_EQUALS(dxadj[myid][i+1],l,PF_ADJ_T);
    dmynedges[myid] = l;
    if (iadjwgt) {
      diadjwgt[myid][i] = iadjwgt[v];
    }
  }

  dl_free(dmynedges);
  dl_free(owner);
  dl_free(lvtx);

  graph = setup_graph(dmynvtxs,dxadj,dadjncy,dadjwgt,diadjwgt,deadjwgt,NULL, \
      nthreads);
  graph->free_xadj = 1;
  graph->free_adjncy = 1;
  if (dadjwgt) {
    graph->free_adjwgt = 1;
  }
  graph->alias = dalias;

  return graph;
}


int free_graph(
    graph_t * graph)
{
  tid_t const nthreads = graph->npar;
  tid_t t;

  /* free what should be free'd */
  if (graph->free_xadj) {
    for (t=0;t<nthreads;++t) {
      dl_free(graph->xadj[t]);
    }
    dl_free(graph->xadj);
  }
  if (graph->free_adjncy) {
    for(t=0;t<nthreads;++t){
      dl_free(graph->adjncy[t]);
    }
    dl_free(graph->adjncy);
  }
  if (graph->free_adjwgt) {
    for(t=0;t<nthreads;++t){
      dl_free(graph->adjwgt[t]);
    }
    dl_free(graph->adjwgt);
  }
  if (graph->free_iadjwgt) {
    for(t=0;t<nthreads;++t){
      dl_free(graph->iadjwgt[t]);
    }
    dl_free(graph->iadjwgt);
  }

  /* free what we've allocated */
  if (graph->eadjwgt) {
    for(t=0;t<nthreads;++t){
      dl_free(graph->eadjwgt[t]);
    }
    dl_free(graph->eadjwgt);
  }
  if (graph->alias) {
    for(t=0;t<nthreads;++t){
      dl_free(graph->alias[t]);
    }
    dl_free(graph->alias);
  }

  if (graph->maxdeg) {
    dl_free(graph->maxdeg);
  }

  /* free thread arrays */
  dl_free(graph->mynvtxs);
  dl_free(graph->mynedges);

  /* free the memory */
  dl_free(graph);

  /* we'll always assume we're successful */
  return 1;
}


adj_t ** build_reverse_adjidx(
    graph_t const * const graph)
{
  vtx_t i, v, k, l, mynvtxs;
  adj_t j;
  tid_t t, o;
  vtx_t ** trans;
  adj_t ** tadj, ** padj, ** txadj;
  adj_t ** radj;

  tid_t const nthreads = graph->npar;
  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy;

  /* build the radj array */
  radj = r_adj_alloc(nthreads);
  txadj = r_adj_alloc(nthreads);
  tadj = r_adj_alloc(nthreads);
  for (t=0;t<nthreads;++t) {
    mynvtxs = graph->mynvtxs[t];
    txadj[t] = adj_duplicate(xadj[t],mynvtxs+1);
    tadj[t] = adj_alloc(graph->mynedges[t]);
    radj[t] = adj_alloc(graph->mynedges[t]);
  }

  /* initial pass for tadj and radj */
  for (t=0;t<nthreads;++t) {
    mynvtxs=graph->mynvtxs[t];
    for (i=0;i<mynvtxs;++i) {
      /* for each vertex */
      for (j=xadj[t][i];j<xadj[t][i+1];++j) {
        k = adjncy[t][j];
        if (k < mynvtxs) {
          o = t;
          l = k;
          v = i;
        } else {
          o = gvtx_to_tid(k,graph->dist);
          l = gvtx_to_lvtx(k,graph->dist);
          v = lvtx_to_gvtx(i,t,graph->dist);
        }
        if (i < k) {
          radj[t][txadj[t][i]] = k;
          tadj[t][txadj[t][i]] = txadj[o][l];
          radj[o][txadj[o][l]] = v;
          tadj[o][txadj[o][l]] = txadj[t][i];
          ++txadj[o][l];
          ++txadj[t][i];
        }
      }
    }
  }
  r_adj_free(txadj,nthreads);

  trans = r_vtx_alloc(nthreads);
  padj = r_adj_alloc(nthreads);
  for (t=0;t<nthreads;++t) {
    trans[t] = vtx_alloc(graph->mynvtxs[t]);
    padj[t] = adj_alloc(graph->mynedges[t]);
  }

  /* turn tadj into a translate array */
  for (t=0;t<nthreads;++t) {
    mynvtxs=graph->mynvtxs[t];
    for (i=0;i<mynvtxs;++i) {
      /* write new indices */
      for (j=xadj[t][i];j<xadj[t][i+1];++j) {
        k = adjncy[t][j];
        if (k < mynvtxs) {
          o = t;
          l = k;
        } else {
          o = gvtx_to_tid(k,graph->dist);
          l = gvtx_to_lvtx(k,graph->dist);
        }
        trans[o][l] = j;
      }
      /* read new indices from old indices */
      for (j=xadj[t][i];j<xadj[t][i+1];++j) {
        k = radj[t][j];
        if (k < mynvtxs) {
          o = t;
          l = k;
        } else {
          o = gvtx_to_tid(k,graph->dist);
          l = gvtx_to_lvtx(k,graph->dist);
        }
        padj[t][j] = trans[o][l];
      }
    }
  }
  r_vtx_free(trans,nthreads);

  /* correct radj -- save it to tadj */
  for (t=0;t<nthreads;++t) {
    mynvtxs=graph->mynvtxs[t];
    for (i=0;i<mynvtxs;++i) {
      for (j=xadj[t][i];j<xadj[t][i+1];++j) {
        radj[t][padj[t][j]] = padj[t][tadj[t][j]];
      }
    }
  }

  dl_free(padj);
  dl_free(tadj);

  return radj;
}


int calc_graph_dist(
    vtx_t const nvtxs, 
    graphdist_t * const dist) 
{
  dist->offset = vtx_uppow2(nvtxs+1);
  dist->mask = dist->offset - 1;
  dist->shift = vtx_downlog2(dist->offset);

  DL_ASSERT(nvtxs < (vtx_t)(1<<dist->shift),"Shift of %d for "PF_VTX_T
      " vertices\n",dist->shift,nvtxs);
  DL_ASSERT_EQUALS(dist->offset,(vtx_t)(1<<dist->shift),PF_VTX_T);

  return 1;
}


int check_graph(
    graph_t const * graph)
{
  vtx_t i,k,v;
  adj_t j,l,maxdeg;
  wgt_t t;
  twgt_t gadjwgt, tadjwgt;
  tid_t myid, o;
  vtx_t mynvtxs;
  adj_t mynedges;

  tid_t const nthreads = graph->npar;
  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy;
  wgt_t const * const * const adjwgt = (wgt_t const * const *)graph->adjwgt;
  wgt_t const * const * const iadjwgt = (wgt_t const * const *)graph->iadjwgt;
  wgt_t const * const * const eadjwgt = (wgt_t const * const *)graph->eadjwgt;

  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    mynedges = graph->mynedges[myid];

    /* xadj checks */
    if (xadj[myid][0] != 0) {
      eprintf("graph->xadj[0] = "PF_ADJ_T"\n",xadj[myid][0]);
      return 0;
    }
    if (mynedges != xadj[myid][mynvtxs]) {
      eprintf("mynedges = "PF_ADJ_T", xadj[nvtxs] = "PF_ADJ_T"\n",
          mynedges,xadj[myid][mynvtxs]);
      return 0;
    }

    /* check max degree */
    maxdeg = 0;
    for (i=0;i<mynvtxs;++i) {
      if (maxdeg < xadj[myid][i+1] - xadj[myid][i]) {
        maxdeg = xadj[myid][i+1] - xadj[myid][i];
      }
    }
    if (maxdeg != graph->maxdeg[myid]) {
      eprintf("Incorrect max degree of "PF_ADJ_T", should be "PF_ADJ_T"\n",
          graph->maxdeg[myid],maxdeg);
      return 0;
    }

    /* edge doubliness checkes */
    for (i=0;i<mynvtxs;++i) {
      t = 0;
      for (j=xadj[myid][i];j<xadj[myid][i+1];++j) {
        k = adjncy[myid][j];
        if (k < mynvtxs) {
          v = i;
          o = myid;
        } else {
          o = gvtx_to_tid(k,graph->dist);
          DL_ASSERT(o<nthreads,"Bad thread id ("PF_TID_T"/"PF_TID_T") for "
              "vertex "PF_VTX_T"/"PF_VTX_T":"PF_VTX_T"\n",o,nthreads,k,
              graph->nvtxs,graph->dist.offset);
          k = gvtx_to_lvtx(k,graph->dist);
          v = lvtx_to_gvtx(i,myid,graph->dist);
        }
        for (l=xadj[o][k];l<xadj[o][k+1];++l) {
          if (adjncy[o][l] == v) {
            break;
          }
        }
        if (l == xadj[o][k+1]) {
          eprintf("Edge from "PF_VTX_T":"PF_TID_T" to "PF_VTX_T":"PF_TID_T" "
              "missing in other direction\n",i,myid,k,o);
          return 0;
        }
        if (adjwgt) {
          if (!dl_near_equal(adjwgt[myid][j],adjwgt[o][l])) {
            eprintf("Edge weight of edge "PF_ADJ_T"={"PF_VTX_T","PF_VTX_T"} is "
                "not symmetric, ["PF_WGT_T","PF_WGT_T"]\n",j,i,k,
                adjwgt[myid][j], adjwgt[o][l]);
            return 0;
          }
          t += adjwgt[myid][j];
        } else {
          t += 1.0;
        }
      }
      if (graph->phantom_edges) {
        if (eadjwgt[myid][i] < t) {
          eprintf("Wrong total degree of "PF_WGT_T" (should be greater than "
              PF_WGT_T") of vertex "PF_VTX_T"\n",eadjwgt[myid][i],t,i);
          return 0;
        }
      } else {
        if (!dl_near_equal(eadjwgt[myid][i],t)) {
          eprintf("Wrong total degree of "PF_WGT_T" (should be "PF_WGT_T") of "
              "vertex "PF_VTX_T"\n",eadjwgt[myid][i],t,i);
          return 0;
        }
      }
    }
  }

  /* calculate tadjwgt */
  tadjwgt = 0.0;
  for (myid=0;myid<nthreads;++myid) {
    tadjwgt += wgt_fa_sum(eadjwgt[myid],graph->mynvtxs[myid]);
  }

  /* calculate gadjwgt */
  gadjwgt = 0.0;
  if (iadjwgt) {
    for (myid=0;myid<nthreads;++myid) {
      gadjwgt += wgt_fa_sum(iadjwgt[myid],graph->mynvtxs[myid]);
    }
  }
  gadjwgt += tadjwgt;

  /* check totals */
  if (!dl_near_equal(tadjwgt,graph->tadjwgt)) {
    eprintf("Incorrect total exposed edge weight of "PF_TWGT_T", actual "
        PF_TWGT_T"\n",graph->tadjwgt,tadjwgt);
    return 0;
  }
  if (!dl_near_equal(gadjwgt,graph->gadjwgt)) {
    eprintf("Incorrect total edge weight of "PF_TWGT_T", actual "
        PF_TWGT_T"\n",graph->gadjwgt,gadjwgt);
    return 0;
  }

  return 1;
}




/******************************************************************************
* PUBLIC PARALLEL FUNCTIONS ***************************************************
******************************************************************************/


graph_t * par_setup_graph(
    vtx_t const mynvtxs,
    adj_t * xadj,
    vtx_t * adjncy,
    wgt_t * adjwgt,
    wgt_t * iadjwgt,
    wgt_t * eadjwgt,
    adj_t const maxdeg,
    dlthread_comm_t const comm)
{
  vtx_t mymd,i;
  twgt_t wgtsum;
  wgt_t * myeadjwgt;
  graph_t * graph;
  graph_t ** r_graph;

  DL_ASSERT(xadj != NULL,"Null xadj passed into setup_graph()");
  DL_ASSERT(adjncy != NULL,"Null xadj passed into setup_graph()");

  tid_t const myid = dlthread_get_id(comm);
  tid_t const nthreads = dlthread_get_nthreads(comm);

  adj_t const mynedges = xadj[mynvtxs];

  r_graph = dlthread_get_shmem(sizeof(graph),comm);

  if (myid == 0) {
    graph = graph_calloc(1); 
    graph->npar = nthreads;
    graph->mynvtxs = vtx_alloc(nthreads); 
    graph->mynedges = adj_alloc(nthreads);
    graph->xadj = r_adj_alloc(nthreads);
    graph->adjncy = r_vtx_alloc(nthreads);
    if (adjwgt) {
      graph->adjwgt = r_wgt_alloc(nthreads);
    }
    if (iadjwgt) {
      graph->iadjwgt = r_wgt_alloc(nthreads);
    }
    graph->eadjwgt = r_wgt_alloc(nthreads);
    graph->maxdeg = adj_alloc(nthreads);

    *r_graph = graph;
  }
  dlthread_barrier(comm);

  graph = *r_graph;

  graph->mynvtxs[myid] = mynvtxs;
  graph->mynedges[myid] = mynedges;

  dlthread_free_shmem(r_graph,comm);

  if (myid == 0) {
    /* set nvtxs and nedges */
    graph->nvtxs = vtx_sum(graph->mynvtxs,nthreads);
    graph->nedges = adj_sum(graph->mynedges,nthreads);
    calc_graph_dist(vtx_max_value(graph->mynvtxs,nthreads),
        &(graph->dist));
  }

  /* these have been asserted to not be NULL */
  graph->xadj[myid] = xadj;
  graph->adjncy[myid] = adjncy;

  graph->xadj[myid][0] = 0;

  if (adjwgt) {
    graph->adjwgt[myid] = adjwgt;
  }
  if (iadjwgt) {
    graph->iadjwgt[myid] = iadjwgt;
  }

  /* if it's NULL or we adjust adjwgt, need to derive it */
  if (eadjwgt) {
    myeadjwgt = graph->eadjwgt[myid] = eadjwgt;
  } else {
    /* calculate eadjwgt */
    myeadjwgt = graph->eadjwgt[myid] = wgt_alloc(mynvtxs);
    if (adjwgt) {
      /* sum adjwgts */
      for (i=0;i<graph->mynvtxs[myid];++i) {
        myeadjwgt[i] = wgt_sum(adjwgt+xadj[i],xadj[i+1]-xadj[i]);
      }
    } else {
      /* use vertex degree */
      for (i=0;i<mynvtxs;++i) {
        myeadjwgt[i] = xadj[i+1]-xadj[i];
      }
    }
  }

  /* handle max degree */
  if (maxdeg != NULL_VTX) {
    graph->maxdeg[myid] = maxdeg;
  } else {
    mymd = 0;
    for (i=0;i<mynvtxs;++i) {
      dl_storemax(mymd,xadj[i+1]-xadj[i]);
    }
    graph->maxdeg[myid] = mymd;
  }

  /* calculate tadjwgt */
  if (adjwgt) {
    wgtsum = twgt_dlthread_sumreduce(wgt_fa_sum(myeadjwgt,mynvtxs),comm);
    if (myid == 0) {
      graph->tadjwgt = wgtsum;
    } 
  } else {
    if (myid == 0) {
      graph->tadjwgt = (wgt_t)graph->nedges; 
    }
  }

  /* calculate gadjwgt */
  if (iadjwgt) {
    wgtsum = twgt_dlthread_sumreduce(wgt_fa_sum(iadjwgt,mynvtxs),comm);
    if (myid == 0) {
      graph->gadjwgt = graph->tadjwgt + wgtsum;
    }
  } else {
    if (myid == 0) {
      graph->gadjwgt = graph->tadjwgt;
    }
  }

  if (myid == 0) {
    graph->snthresh = __calc_snthresh(graph);
  }

  dlthread_barrier(comm);

  if (myid == 0) {
    dprintf("Setup a graph with "PF_VTX_T" vertices, "PF_ADJ_T" edges, split "
        "among "PF_TID_T" threads\n",graph->nvtxs,graph->nedges,
        nthreads);
  }
     
  return graph;
}


void par_free_graph(
    graph_t * graph,
    dlthread_comm_t const comm)
{
  tid_t const myid = dlthread_get_id(comm);

  /* free what should be free'd */
  if (graph->free_xadj) {
    dl_free(graph->xadj[myid]);
  }
  if (graph->free_adjncy) {
    dl_free(graph->adjncy[myid]);
  }
  if (graph->free_adjwgt) {
    dl_free(graph->adjwgt[myid]);
  }
  if (graph->free_iadjwgt) {
    dl_free(graph->iadjwgt[myid]);
  }
  if (graph->eadjwgt) {
    dl_free(graph->eadjwgt[myid]);
  }
  if (graph->alias) {
    dl_free(graph->alias[myid]);
  }

  dlthread_barrier(comm);

  if (myid == 0) {
    if (graph->free_xadj) {
      dl_free(graph->xadj);
    }
    if (graph->free_adjncy) {
      dl_free(graph->adjncy);
    }
    if (graph->free_adjwgt) {
      dl_free(graph->adjwgt);
    }
    if (graph->free_iadjwgt) {
      dl_free(graph->iadjwgt);
    }
    if (graph->eadjwgt) {
      dl_free(graph->eadjwgt);
    }
    if (graph->alias) {
      dl_free(graph->alias);
    }
    if (graph->maxdeg) {
      dl_free(graph->maxdeg);
    }
    if (graph->mynvtxs) {
      dl_free(graph->mynvtxs);
    }
    if (graph->mynedges) {
      dl_free(graph->mynedges);
    }
    /* free the memory */
    dl_free(graph);
  }

  dlthread_barrier(comm);
}




#endif
