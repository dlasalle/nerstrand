/**
 * @file sparsen.c
 * @brief Edge spasrening functions
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-07-25
 */




#ifndef NERSTRAND_SPARSEN_C
#define NERSTRAND_SPARSEN_C




#include "sparsen.h"




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


static const adj_t NBUCKETS = 1024;
static const adj_t RRANGE = 128;




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static inline int __myedge(
    const vtx_t i, 
    const vtx_t j)
{
  if ((i + j) % 2 == 0) {
    return i > j;
  } else {
    return i < j;
  }
}


/**
 * @brief Sparsify the graphs edges randomly
 *
 * @param objective The objective providign parameters for sparsification
 * @param graph The graph to sparsen
 * @param targetnedges The target number of edges to leave
 * @param radj The array to place edges selected for removal in
 *
 * @return NERSTRAND_SUCCESS, or an error code
 */
static int __sparsify_RANDOM(
    objective_t * const objective, 
    const graph_t * const graph, 
    const adj_t targetnedges, 
    adj_t ** const radj)
{
  unsigned int seed;
  vtx_t i, g, l, mynvtxs;
  adj_t j, k, nadj, mv, mynradj, mynedges;
  tid_t myid, o;
  adj_t * myxadj, * rc, * myradj; 
  vtx_t * myadjncy;
  small_t * myrating, ** rating;
  
  const tid_t nthreads = graph->npar;
  const adj_t over = graph->nedges > targetnedges ? 
      graph->nedges - targetnedges : 0;

  rating = r_small_dalloc(graph->mynedges,sizeof(adj_t),nthreads);
  rc = adj_calloc(NBUCKETS);

  if (over > 0) {
    for (myid=0;myid<nthreads;++myid) {
      mynvtxs = graph->mynvtxs[myid];
      myxadj = graph->xadj[myid];
      myadjncy = graph->adjncy[myid];
      myrating = rating[myid];
      /* assign numbers to edges */
      for (i=0;i<mynvtxs;++i) {
        for (j=myxadj[i];j<myxadj[i+1];++j) {
          k = myadjncy[j];
          if (k < mynvtxs) {
            o = myid;
            l = k;
            k = lvtx_to_gvtx(l,o,graph->dist);
            g = i;
          } else {
            o = gvtx_to_tid(k,graph->dist);
            l = gvtx_to_lvtx(k,graph->dist);
            g = lvtx_to_gvtx(i,myid,graph->dist);
          }
          seed = ((unsigned int)((g^k)^(g+k)))%NBUCKETS;
          myrating[j] = seed;
          ++rc[seed];
        }
      }
    }
    /* decide what edges to prune */
    mv = nadj = 0;
    for (j=0;j<NBUCKETS;++j) {
      if ((nadj += rc[j]) > over) {
        mv = j;
        break;
      }
    }
    for (myid=0;myid<nthreads;++myid) {
      mynedges = graph->mynedges[myid];
      myrating = rating[myid];
      myradj = radj[myid];
      mynradj = 0;
      /* assign numbers to edges */
      for (j=0;j<mynedges;++j) {
        if (myrating[j] <= mv) {
          myradj[mynradj++] = j;
        }
      }
      myradj[mynradj] = mynedges; /* the cap */
    }
  } else {
    for (myid=0;myid<nthreads;++myid) {
      myradj = radj[myid];
      myradj[0] = graph->mynedges[myid];
    }
  }

  r_small_free(rating,nthreads);
  dl_free(rc);

  return NERSTRAND_SUCCESS;
}


/**
 * @brief Remove the lightest edges 
 *
 * @param objective The objective providing parameters for sparsification
 * @param graph The graph to sparsen
 * @param targetnedges The target number of edges to leave
 * @param radj The array to place edges selected for removal in
 *
 * @return The edges flagged to be removed 
 */
static int __sparsify_LIGHT(
    objective_t * const objective, 
    const graph_t * const graph, 
    const adj_t targetnedges, 
    adj_t ** const radj)
{
  vtx_t i, g, l, mynvtxs;
  adj_t j, k, nadj, mv, seed, mynradj, mynedges;
  tid_t myid, o;
  adj_t * myxadj, * rc, * myradj; 
  vtx_t * myadjncy;
  wgt_t * myadjwgt;
  small_t * myrating, ** rating;
  
  const tid_t nthreads = graph->npar;
  const adj_t over = graph->nedges > targetnedges ? 
      graph->nedges - targetnedges : 0;
  const wgt_t avgwgt = graph->tadjwgt / graph->nedges;
  const wgt_t bsize = (avgwgt*2.0) / NBUCKETS;

  if (graph->adjwgt == NULL) {
    return __sparsify_RANDOM(objective,graph,targetnedges,radj);
  }

  rating = r_small_dalloc(graph->mynedges,sizeof(adj_t),nthreads);
  rc = adj_calloc(NBUCKETS);

  if (over > 0) {
    for (myid=0;myid<nthreads;++myid) {
      mynvtxs = graph->mynvtxs[myid];
      myxadj = graph->xadj[myid];
      myadjncy = graph->adjncy[myid];
      myadjwgt = graph->adjwgt[myid];
      myrating = rating[myid];
      /* assign numbers to edges */
      for (i=0;i<mynvtxs;++i) {
        for (j=myxadj[i];j<myxadj[i+1];++j) {
          k = myadjncy[j];
          if (k < mynvtxs) {
            o = myid;
            l = k;
            k = lvtx_to_gvtx(l,o,graph->dist);
            g = i;
          } else {
            o = gvtx_to_tid(k,graph->dist);
            l = gvtx_to_lvtx(k,graph->dist);
            g = lvtx_to_gvtx(i,myid,graph->dist);
          }
          seed = (adj_t)((myadjwgt[j]/bsize)+(((g^k)^(g+k))%(RRANGE)));
          if (seed >= NBUCKETS) {
            seed = NBUCKETS-1;
          }
          myrating[j] = seed;
          ++rc[seed];
        }
      }
    }
    /* decide what edges to prune */
    mv = nadj = 0;
    for (j=0;j<NBUCKETS;++j) {
      if ((nadj += rc[j]) > over) {
        mv = j;
        break;
      }
    }
    for (myid=0;myid<nthreads;++myid) {
      mynedges = graph->mynedges[myid];
      myrating = rating[myid];
      myradj = radj[myid];
      mynradj = 0;
      /* assign numbers to edges */
      for (j=0;j<mynedges;++j) {
        if (myrating[j] <= mv) {
          myradj[mynradj++] = j;
        }
      }
      myradj[mynradj] = mynedges; /* the cap */
    }
  } else {
    for (myid=0;myid<nthreads;++myid) {
      myradj = radj[myid];
      myradj[0] = graph->mynedges[myid];
    }
  }

  r_small_free(rating,nthreads);
  dl_free(rc);

  return NERSTRAND_SUCCESS;
}


/**
 * @brief Remove the highest degree edges 
 *
 * @param objective The objective providing parameters for sparsification
 * @param graph The graph to sparsen
 * @param targetnedges The target number of edges to leave
 * @param radj The array to place edges selected for removal in
 *
 * @return The edges flagged to be removed 
 */
static int __sparsify_DEGREE(
    objective_t * const objective, 
    const graph_t * const graph, 
    const adj_t targetnedges, 
    adj_t * const * const radj)
{
  vtx_t i, g, l, mynvtxs;
  adj_t j, k, nadj, mv, seed, deg, mydeg, mynradj, mynedges;
  tid_t myid, o;
  adj_t * myxadj, * rc, * myradj;
  small_t * myrating, ** rating; 
  vtx_t * myadjncy;
  
  const tid_t nthreads = graph->npar;
  const adj_t over = graph->nedges > targetnedges ? 
      graph->nedges - targetnedges : 0;
  const adj_t maxdeg = adj_max_value(graph->maxdeg,nthreads); 
  const wgt_t bsize = (2.0*maxdeg) / NBUCKETS;

  rating = r_small_dalloc(graph->mynedges,sizeof(adj_t),nthreads);
  rc = adj_calloc(NBUCKETS);

  if (over > 0) {
    for (myid=0;myid<nthreads;++myid) {
      mynvtxs = graph->mynvtxs[myid];
      myxadj = graph->xadj[myid];
      myadjncy = graph->adjncy[myid];
      myrating = rating[myid];
      /* assign numbers to edges */
      for (i=0;i<mynvtxs;++i) {
        mydeg = myxadj[i+1] - myxadj[i];
        for (j=myxadj[i];j<myxadj[i+1];++j) {
          k = myadjncy[j];
          if (k < mynvtxs) {
            o = myid;
            l = k;
            k = lvtx_to_gvtx(l,o,graph->dist);
            g = i;
          } else {
            o = gvtx_to_tid(k,graph->dist);
            l = gvtx_to_lvtx(k,graph->dist);
            g = lvtx_to_gvtx(i,myid,graph->dist);
          }
          deg = maxdeg - (mydeg + (graph->xadj[o][l+1]-graph->xadj[o][l]));
          seed = (adj_t)((deg/bsize) + (((g^k)^(g+k))%(RRANGE)));
          if (seed >= NBUCKETS) {
            seed = NBUCKETS-1;
          }
          myrating[j] = seed;
          ++rc[seed];
        }
      }
    }
    /* decide what edges to prune */
    mv = nadj = 0;
    for (j=0;j<NBUCKETS;++j) {
      if ((nadj += rc[j]) > over) {
        mv = j;
        break;
      }
    }
    for (myid=0;myid<nthreads;++myid) {
      mynedges = graph->mynedges[myid];
      myrating = rating[myid];
      myradj = radj[myid];
      mynradj = 0;
      /* assign numbers to edges */
      for (j=0;j<mynedges;++j) {
        if (myrating[j] <= mv) {
          myradj[mynradj++] = j;
        }
      }
      myradj[mynradj] = mynedges; /* the cap */
    }
  } else {
    for (myid=0;myid<nthreads;++myid) {
      myradj = radj[myid];
      myradj[0] = graph->mynedges[myid];
    }
  }

  r_small_free(rating,nthreads);
  dl_free(rc);

  return NERSTRAND_SUCCESS;
}


/**
 * @brief Remove the selected edges from the graph and drop their
 * weight.
 *
 * @param radj The list of edges to remove
 * @param graph The graph to remove edges from
 *
 * @return NERSTRAND_SUCCESS, or an error code if unsuccessful
 */
static int __remove_DROP(
    const adj_t * const * const radj, 
    graph_t * const graph)
{
  vtx_t i, mynvtxs;
  adj_t j, k, maxdeg, deg, nedges;
  wgt_t dwgt;
  tid_t myid;
  adj_t * myxadj;
  vtx_t * myadjncy;
  wgt_t * myadjwgt, * myeadjwgt;
  const adj_t * myradj;

  const tid_t nthreads = graph->npar;

  for (myid=0;myid<nthreads;++myid) {
    /* alias things that are constant in this iteration */
    mynvtxs = graph->mynvtxs[myid];
    myxadj = graph->xadj[myid];
    myadjncy = graph->adjncy[myid];
    myadjwgt = graph->adjwgt[myid];
    myeadjwgt = graph->eadjwgt[myid];
    myradj = radj[myid];
    /* zero out variables for this iteration */
    dwgt = 0.0;
    nedges = 0;
    maxdeg = 0;
    j = 0;
    k = 0;
    /* loop through each vertex */
    for (i=0;i<mynvtxs;++i) {
      /* loop through the edges of this vertex */
      while (j<myxadj[i+1]) {
        /* check to see if this edge is removed */
        if (j == myradj[k]) {
          /* remove the edge */ 
          myeadjwgt[i] -= myadjwgt[j];
          dwgt += myadjwgt[j];
          ++k;
        } else {
          /* shift the edge to the next slot */
          myadjncy[nedges] = myadjncy[j];
          myadjwgt[nedges] = myadjwgt[j];
          ++nedges;
        }
        ++j;
      }
      myxadj[i+1] = nedges;
      deg = nedges - myxadj[i];
      if (deg > maxdeg) {
        maxdeg = deg; 
      }
    }
    DL_ASSERT_EQUALS(myradj[k],graph->mynedges[myid],PF_ADJ_T);
    graph->mynedges[myid] = nedges;
    graph->maxdeg[myid] = maxdeg;
    graph->gadjwgt -= dwgt;
    graph->tadjwgt -= dwgt;
  }
  graph->nedges = adj_sum(graph->mynedges,nthreads);

  return NERSTRAND_SUCCESS;
}


/**
 * @brief Remove the selected edges from the graph and add their
 * weight as internal.
 *
 * @param radj The list of edges to remove
 * @param graph The graph to remove edges from
 *
 * @return NERSTRAND_SUCCESS, or an error code if unsuccessful
 */
static int __remove_LOOP(
    const adj_t * const * const radj, 
    graph_t * const graph)
{
  vtx_t i, mynvtxs;
  adj_t j, k, nedges, maxdeg, deg;
  wgt_t dwgt, iwgt;
  tid_t myid;
  adj_t * myxadj;
  vtx_t * myadjncy;
  wgt_t * myadjwgt, * myeadjwgt;
  const adj_t * myradj;

  const tid_t nthreads = graph->npar;

  for (myid=0;myid<nthreads;++myid) {
    /* alias things that are constant in this iteration */
    mynvtxs = graph->mynvtxs[myid];
    myxadj = graph->xadj[myid];
    myadjncy = graph->adjncy[myid];
    myadjwgt = graph->adjwgt[myid];
    myeadjwgt = graph->eadjwgt[myid];
    myradj = radj[myid];
    /* zero out variables for this iteration */
    dwgt = 0.0;
    iwgt = 0.0;
    nedges = 0;
    maxdeg = 0;
    j = 0;
    k = 0;
    /* loop through each vertex */
    for (i=0;i<mynvtxs;++i) {
      /* loop through the edges of this vertex */
      while (j<myxadj[i+1]) {
        /* check to see if this edge is removed */
        if (j == myradj[k]) {
          /* remove the edge */ 
          if (graph->iadjwgt) {
            graph->iadjwgt[myid][i] += myadjwgt[j];
            iwgt += myadjwgt[j];
          }
          myeadjwgt[i] -= myadjwgt[j];
          dwgt += myadjwgt[j];
          ++k;
        } else {
          /* shift the edge to the next slot */
          myadjncy[nedges] = myadjncy[j];
          myadjwgt[nedges] = myadjwgt[j];
          ++nedges;
        }
        ++j;
      }
      myxadj[i+1] = nedges;
      deg = nedges - myxadj[i];
      if (deg > maxdeg) {
        maxdeg = deg; 
      }
    }
    DL_ASSERT_EQUALS(myradj[k],graph->mynedges[myid],PF_ADJ_T);
    graph->mynedges[myid] = nedges;
    graph->maxdeg[myid] = maxdeg;
    graph->tadjwgt -= dwgt;
    graph->gadjwgt -= dwgt - iwgt;
  }
  graph->nedges = adj_sum(graph->mynedges,nthreads);

  return NERSTRAND_SUCCESS;
}


/**
 * @brief Remove the selected edges from the graph and redistribute their
 * weight across remaining edges.
 *
 * @param radj The list of edges to remove
 * @param graph The graph to remove edges from
 *
 * @return NERSTRAND_SUCCESS, or an error code if unsuccessful
 */
static int __remove_DISTRIBUTE(
    const adj_t * const * const radj, 
    graph_t * const graph)
{
  vtx_t i, k, l, mynvtxs;
  adj_t j, nedges, maxdeg, deg;
  tid_t myid,o;
  adj_t * myxadj;
  vtx_t * myadjncy;
  wgt_t * myadjwgt, * myeadjwgt;
  wgt_t ** dist;

  const adj_t * myradj;
  
  const tid_t nthreads = graph->npar;

  dist = r_wgt_dcalloc(graph->mynvtxs,sizeof(vtx_t),nthreads);

  for (myid=0;myid<nthreads;++myid) {
    /* alias things that are constant in this iteration */
    mynvtxs = graph->mynvtxs[myid];
    myxadj = graph->xadj[myid];
    myadjncy = graph->adjncy[myid];
    myadjwgt = graph->adjwgt[myid];
    myeadjwgt = graph->eadjwgt[myid];
    myradj = radj[myid];
    /* zero out variables for this iteration */
    nedges = 0;
    maxdeg = 0;
    j = 0;
    k = 0;
    /* loop through each vertex */
    for (i=0;i<mynvtxs;++i) {
      /* loop through the edges of this vertex */
      while (j<myxadj[i+1]) {
        /* check to see if this edge is removed */
        if (j == myradj[k]) {
          /* remove the edge */ 
          myeadjwgt[i] -= myadjwgt[j];
          dist[myid][i] += (myadjwgt[j]/2.0);
          ++k;
        } else {
          /* shift the edge to the next slot */
          myadjncy[nedges] = myadjncy[j];
          myadjwgt[nedges] = myadjwgt[j];
          ++nedges;
        }
        ++j;
      }
      myxadj[i+1] = nedges;
      deg = nedges - myxadj[i];
      if (deg > maxdeg) {
        maxdeg = deg; 
      }
    }
    DL_ASSERT_EQUALS(myradj[k],graph->mynedges[myid],PF_ADJ_T);
    graph->mynedges[myid] = nedges;
    graph->maxdeg[myid] = maxdeg;
  }
  graph->nedges = adj_sum(graph->mynedges,nthreads);
  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    myxadj = graph->xadj[myid];
    for (i=0;i<mynvtxs;++i) {
      if (myxadj[i+1] > myxadj[i]) {
        dist[myid][i] = floor(dist[myid][i]/(myxadj[i+1] - myxadj[i]));
      } else if (graph->iadjwgt) {
        graph->iadjwgt[myid][i] += dist[myid][i]*2.0;
      }
    }
  }
  graph->gadjwgt = graph->tadjwgt = 0;
  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    myxadj = graph->xadj[myid];
    myadjncy = graph->adjncy[myid];
    myadjwgt = graph->adjwgt[myid];
    for (i=0;i<mynvtxs;++i) {
      /* re-sum to be numerically stable */
      graph->eadjwgt[myid][i] = 0;
      for (j=myxadj[i];j<myxadj[i+1];++j) {
        k = myadjncy[j];
        if (k < mynvtxs) {
          o = myid;
          l = k;
        } else {
          o = gvtx_to_tid(k,graph->dist);
          l = gvtx_to_lvtx(k,graph->dist);
        }
        myadjwgt[j] += dist[myid][i] + dist[o][l];
        graph->eadjwgt[myid][i] += myadjwgt[j];
      }
      if (graph->iadjwgt) {
        graph->gadjwgt += graph->iadjwgt[myid][i];
      }
      graph->gadjwgt += graph->eadjwgt[myid][i];
      graph->tadjwgt += graph->eadjwgt[myid][i];
    }
  }

  r_wgt_free(dist,nthreads);

  return NERSTRAND_SUCCESS;
}


/**
 * @brief Remove the selected edges from the graph but keep track of the
 * vertex's extern edge weight.
 *
 * @param radj The list of edges to remove
 * @param graph The graph to remove edges from
 *
 * @return NERSTRAND_SUCCESS, or an error code if unsuccessful
 */
static int __remove_PHANTOM(
    const adj_t * const * const radj, 
    graph_t * const graph)
{
  vtx_t i, k, mynvtxs;
  adj_t j, nedges, maxdeg, deg;
  tid_t myid;
  adj_t * myxadj;
  vtx_t * myadjncy;
  wgt_t * myadjwgt;
  const adj_t * myradj;
  
  const tid_t nthreads = graph->npar;

  for (myid=0;myid<nthreads;++myid) {
    /* alias things that are constant in this iteration */
    mynvtxs = graph->mynvtxs[myid];
    myxadj = graph->xadj[myid];
    myadjncy = graph->adjncy[myid];
    myadjwgt = graph->adjwgt[myid];
    myradj = radj[myid];
    /* zero out variables for this iteration */
    nedges = 0;
    maxdeg = 0;
    j = 0;
    k = 0;
    /* loop through each vertex */
    for (i=0;i<mynvtxs;++i) {
      /* loop through the edges of this vertex */
      while (j<myxadj[i+1]) {
        /* check to see if this edge is removed */
        if (j == myradj[k]) {
          /* remove the edge */ 
          ++k;
        } else {
          /* shift the edge to the next slot */
          myadjncy[nedges] = myadjncy[j];
          myadjwgt[nedges] = myadjwgt[j];
          ++nedges;
        }
        ++j;
      }
      myxadj[i+1] = nedges;
      deg = nedges - myxadj[i];
      if (deg > maxdeg) {
        maxdeg = deg; 
      }
    }
    DL_ASSERT_EQUALS(myradj[k],graph->mynedges[myid],PF_ADJ_T);
    graph->mynedges[myid] = nedges;
    graph->maxdeg[myid] = maxdeg;
  }
  graph->nedges = adj_sum(graph->mynedges,nthreads);
  graph->phantom_edges = 1;

  return NERSTRAND_SUCCESS;
}


/**
 * @brief Check that an edge removal array is symmetric
 *
 * @param graph The graph with removed edges
 * @param radj The edge removal array
 *
 * @return !0 on success
 */
static int __check_radj(
    const graph_t * const graph, 
    const adj_t * const * const radj)
{
  vtx_t i,k,mynvtxs,l;
  adj_t j,jj,x,xx;
  tid_t myid,o;

  const tid_t nthreads = graph->npar;
  const adj_t * const * const xadj = (const adj_t * const *)graph->xadj;
  const vtx_t * const * const adjncy = (const vtx_t * const *)graph->adjncy;

  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    j = 0;
    while (1) {
      if (radj[myid][j] < xadj[myid][mynvtxs]) {
        if (radj[myid][j] >= radj[myid][j+1]) {
          eprintf("Out of order edge index in radj "PF_ADJ_T" follwed by "
              PF_ADJ_T"\n",radj[myid][j],radj[myid][j+1]);
          return 0;
        }
      } else if (radj[myid][j] == xadj[myid][mynvtxs]) {
        break;
      } else {
        eprintf("Invalid edge index in radj "PF_ADJ_T"/"PF_ADJ_T"\n",
            radj[myid][j],xadj[myid][mynvtxs]);
        return 0;
      }
    }
    j = 0;
    x = 0;
    for (i=0;i<mynvtxs;++i) {
      while (j<xadj[myid][i+1]) {
        if (j == radj[myid][x]) {
          k = adjncy[myid][j];
          if (k < mynvtxs) {
            o = myid;
            l = k;
          } else {
            o = gvtx_to_tid(k,graph->dist);
            l = gvtx_to_lvtx(k,graph->dist);
          }
          for (jj=xadj[o][l];jj<xadj[o][l+1];++jj) {
            if (adjncy[o][jj] == i) {
              xx = 0;
              while (radj[o][xx] < xadj[o][graph->mynvtxs[o]]) {
                if (radj[o][xx] == jj) {
                  goto FOUND;
                }
                ++xx;
              }
            }
          }
          eprintf("Edge "PF_VTX_T":"PF_TID_T"-"PF_TID_T":"PF_VTX_T" ["PF_ADJ_T
              "] is not removed in other direction\n",i,myid,o,l,j);
          return 0;
          FOUND:
          ++x;
        }
        ++j;
      }
    }
  }
  return 1;
}




/******************************************************************************
* PUBLIC SERIAL FUNCTIONS *****************************************************
******************************************************************************/


int sparsen_graph(
    objective_t * objective, 
    graph_t * graph, 
    const adj_t targetnedges)
{
  adj_t ** radj = NULL;

  const tid_t nthreads = graph->npar;
  const adj_t nedges = graph->nedges;

  dl_start_timer(&(objective->timers.sparsification));

  if (objective->sparsifytype != NERSTRAND_SPARSIFY_NONE) {
    radj = r_adj_dalloc(graph->mynedges,sizeof(adj_t),nthreads);

    dl_start_timer(&(objective->timers.edgeselection));
    switch (objective->sparsifytype) {
      case NERSTRAND_SPARSIFY_RANDOM:
        __sparsify_RANDOM(objective,graph,targetnedges,radj);
        break;
      case NERSTRAND_SPARSIFY_LIGHT:
        __sparsify_LIGHT(objective,graph,targetnedges,radj);
        break;
      case NERSTRAND_SPARSIFY_DEGREE:
        __sparsify_DEGREE(objective,graph,targetnedges,radj);
        break;
      default:
        eprintf("Unknown sparsify type '%d'\n",objective->sparsifytype);
    }
    dl_stop_timer(&(objective->timers.edgeselection));

    DL_ASSERT(__check_radj(graph,(const adj_t * const *)radj) == 1,
        "Unsymmetric radj vector\n");

    dl_start_timer(&(objective->timers.edgeremoval));
    switch (objective->edgeremtype) {
      case NERSTRAND_EDGEREMOVAL_DROP:
        __remove_DROP((const adj_t * const *)radj,graph);
        break;
      case NERSTRAND_EDGEREMOVAL_LOOP:
        __remove_LOOP((const adj_t * const *)radj,graph);
        break;
      case NERSTRAND_EDGEREMOVAL_DISTRIBUTE:
        __remove_DISTRIBUTE((const adj_t * const *)radj,graph);
        break;
      case NERSTRAND_EDGEREMOVAL_PHANTOM:
        __remove_PHANTOM((const adj_t * const *)radj,graph);
        break;
      default:
        eprintf("Unknown edge removal type '%d'\n",objective->edgeremtype);
    }
    dl_stop_timer(&(objective->timers.edgeremoval));

    r_adj_free(radj,nthreads);

    vprintf(objective->verbosity,NERSTRAND_VERBOSITY_HIGH,
        "Sparsening left "PF_ADJ_T"/"PF_ADJ_T" of "PF_ADJ_T" edges\n",
        graph->nedges,targetnedges,nedges);
  }

  DL_ASSERT(check_graph(graph) == 1, "Bad graph after sparsening\n");

  dl_stop_timer(&(objective->timers.sparsification));

  return NERSTRAND_SUCCESS;
}





#endif
