/**
 * @file initialclustering.c
 * @brief Function for initial clustering
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-10-08
 */




#ifndef NERSTRAND_INITIALCLUSTERING_C
#define NERSTRAND_INITIALCLUSTERING_C




#include "initialclustering.h"
#include "uncoarsen.h"
#include "refine.h"
#include "newman.h"




/******************************************************************************
* PRIVATE CONSTANTS ***********************************************************
******************************************************************************/


static cid_t const ENQUEUED = -2;
static vtx_t const MAX_FIXES = 8;
static vtx_t const MIN_SNVTXS = 32;
static vtx_t const HILLSIZE = 0; //32;
static vtx_t const MAX_VTX_EDGES = 1048576;




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLPQ_PREFIX vv
#define DLPQ_KEY_T vtx_t
#define DLPQ_VAL_T vtx_t
#define DLPQ_STATIC
#include "dlpq_headers.h"
#undef DLPQ_STATIC
#undef DLPQ_VAL_T
#undef DLPQ_KEY_T
#undef DLPQ_PREFIX


#define DLPQ_PREFIX mv 
#define DLPQ_KEY_T mod_t
#define DLPQ_VAL_T vtx_t
#define DLPQ_STATIC
#include "dlpq_headers.h"
#undef DLPQ_STATIC
#undef DLPQ_VAL_T
#undef DLPQ_KEY_T
#undef DLPQ_PREFIX


#define DLDJSET_TYPE_T vtx_t
#define DLDJSET_PREFIX vtx
#define DLDJSET_STATIC
#include "dldjset_headers.h"
#undef DLDJSET_STATIC
#undef DLDJSET_PREFIX
#undef DLDJSET_TYPE_T




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


/**
 * @brief Determine the priority for merging two clusters.
 *
 * @param ewgt
 * @param dega
 * @param degb
 * @param invgadjwgt
 *
 * @return 
 */
static inline mod_t __merge_priority(
    wgt_t const ewgt, 
    wgt_t const dega, 
    wgt_t const degb, 
    wgt_t const invgadjwgt)
{
  /* Q_{a} = d_int{a} - (d_{a}^2)/m
     Q_{b} = d_int{b} - (d_{b}^2)/m
     Q_{a,b} = d_int{a} + d_int{b} + 2*ewgt - (d_{a}^2+d_{a}d_{b}+d_{b}^2)/m
     Q_merge = Q_{a,b} - Q_{a} - Q_{b}
     Q_merge = 2*ewgt - (d_{a}d_{b})/m */
  return (2*ewgt) - (dega*degb*invgadjwgt);
}


/**
 * @brief Determine the priority for moving a vertex from one cluster to
 * another.
 *
 * @param dega
 * @param degb
 * @param degv
 * @param awgt
 * @param bwgt
 * @param invgadjwgt
 *
 * @return 
 */
static inline mod_t __move_priority(
    wgt_t const dega, 
    wgt_t const degb, 
    wgt_t const degv, 
    wgt_t const awgt, 
    wgt_t const bwgt, 
    wgt_t const invgadjwgt)
{
  /* Q_{a} = d_int{a} - (d_{a}^2)/m
     Q_{b} = d_int{b} - (d_{b}^2)/m
     Q_{a-v} = d_int{a}-d_int{v}-2*d(a){v} - (d_{a}^2-d_{a}d_{v}+d_{v}^2)/m
     Q_{b+v} = d_int{b}+d_int{v}+2*d(b){v} - (d_{b}^2+d_{b}d_{v}+d_{v}^2)/m
     Q_move = Q_{a-v} + Q_{b+v} - (Q_{a} + Q_{b})
     Q_move = 2*d(b){v} - 2*d(a){v} +  (d_{v}(d_{a}- d_{b}-d_{v}))/m 
   */
  return bwgt - awgt + ((degv*(/*dega*/-degb-degv))*invgadjwgt);
}


/**
 * @brief Add a vertex to moved vertex list. 
 *
 * @param a
 * @param m
 * @param ind
 * @param ptr
 *
 * @return 
 */
static inline vtx_t __move_vertex(
    vtx_t const a, 
    vtx_t const m, 
    vtx_t * const ind, 
    vtx_t * const ptr)
{
  DL_ASSERT(ptr[a] >= m, "Bad parameters passed to move vertex (a="PF_VTX_T \
      " m="PF_VTX_T")\n",ptr[a],m);
  ind[ptr[a]] = ind[m];
  ptr[ind[m]] = ptr[a];
  ind[m] = a;
  ptr[a] = m;
  return m+1;
}


/**
 * @brief Remove a vertex from the moved vertex list.
 *
 * @param a
 * @param m
 * @param ind
 * @param ptr
 *
 * @return 
 */
static inline vtx_t __unmove_vertex(
    vtx_t const a, 
    vtx_t const m, 
    vtx_t * const ind, 
    vtx_t * const ptr)
{
  DL_ASSERT(ptr[a] >= m, "Bad parameters passed to move vertex (a="PF_VTX_T \
      " m="PF_VTX_T")\n",ptr[a],m);

  ind[ptr[a]] = ind[m];
  ptr[ind[m]] = ptr[a];
  ind[m] = a;
  ptr[a] = m;
  return m+1;
}


/**
 * @brief Assign each cluster to be a singleton vertex.
 *
 * @param mgraph
 * @param seed
 *
 * @return 
 */
static clustering_t * __cluster_VTX(
    mgraph_t const * const mgraph, 
    unsigned int seed)
{
  tid_t t;
  vtx_t mynvtxs,start;
  clustering_t * clustering;

  graph_t const * const graph = mgraph->graph;
  tid_t const nthreads = graph->npar;

  cid_t ** const where = r_cid_alloc(nthreads);

  start=0;
  for (t=0;t<nthreads;++t) {
    mynvtxs = graph->mynvtxs[t];
    where[t] = cid_alloc(mynvtxs);
    cid_incset(where[t],start,1,mynvtxs);
    start += mynvtxs;
  }

  clustering = setup_clustering(graph->nvtxs,where,graph);
  clustering->free_clusters = 1;
  clustering->free_where = 1;

  return clustering; 
}


/**
 * @brief Randomly assign a cluster to each vertex.
 *
 * @param nclusters The number of clusters to generate.
 * @param mgraph The graph.
 * @param seed The seed to use.
 *
 * @return 
 */
static clustering_t * __cluster_RAND(
    cid_t const nclusters,
    mgraph_t const * const mgraph, 
    unsigned int seed)
{
  tid_t myid;
  clustering_t * clustering;

  graph_t const * const graph = mgraph->graph;
  vtx_t const * const mynvtxs = graph->mynvtxs;
  tid_t const nthreads = graph->npar;

  cid_t ** const where = r_cid_dalloc(mynvtxs,sizeof(vtx_t),nthreads);

  for (myid=0;myid<nthreads;++myid) {
    cid_fill_rand_r(0,nclusters,where[myid],mynvtxs[myid],&seed);
  }

  clustering = setup_clustering(nclusters,where,graph);
  clustering->free_clusters = 1;
  clustering->free_where = 1;

  return clustering; 
}


/**
 * @brief Perform a seed based initial cluster, where by high degree vertices
 * are select to server as "seeds" for growing clusters.
 *
 * @param nclusters The number of clusters to generate.
 * @param mgraph The graph.
 * @param seed The random seed to use.
 *
 * @return The new cluster.
 */
static clustering_t * __cluster_SEED(
    cid_t const nclusters,
    mgraph_t const * const mgraph, 
    unsigned int seed)
{
  vtx_t k,l,v, mynvtxs, deg, g, nclustered, nprocessed;
  adj_t j;
  tid_t myid,o;
  cid_t c, other;
  clustering_t * clustering;

  graph_t const * const graph = mgraph->graph;
  tid_t const nthreads = graph->npar;

  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy; 

  cid_t ** const where = r_cid_alloc(nthreads);
  vv_pq_t * q = vv_pq_create(graph->dist.offset, \
      max_gvtx(graph)+1);
  vtx_t * horizon = vtx_alloc(graph->nvtxs);
  vtx_t * counts = vtx_calloc(nclusters);

  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid]; 
    where[myid] = cid_init_alloc(NULL_CID,mynvtxs);
    for (v=0;v<mynvtxs;++v) {
      deg = graph->xadj[myid][v+1]-graph->xadj[myid][v];
      deg += vtx_rand_r(0,deg,&seed); /* and randomness */
      g = lvtx_to_gvtx(v,myid,graph->dist);
      vv_pq_push(deg,g,q);
    }
  }

  nclustered = 0;
  for (c=0;c<nclusters;++c) {
    g = vv_pq_pop(q);
    v = gvtx_to_lvtx(g,graph->dist);
    myid = gvtx_to_tid(g,graph->dist);
    where[myid][v] = c;
    ++counts[c];
    horizon[nclustered++] = g;
  }

  nprocessed = 0;
  while (nclustered < graph->nvtxs) {
    g = horizon[nprocessed++];
    v = gvtx_to_lvtx(g,graph->dist);
    myid = gvtx_to_tid(g,graph->dist);

    mynvtxs = graph->mynvtxs[myid]; 
    c = where[myid][v];

    for (j=xadj[myid][v];j<xadj[myid][v+1];++j) {
      k = adjncy[myid][j];
      if (k < mynvtxs) {
        o = myid; 
        l = k;
        k = lvtx_to_gvtx(l,myid,graph->dist);
      } else {
        o = gvtx_to_tid(k,graph->dist);
        l = gvtx_to_lvtx(k,graph->dist);
      }
      other = where[o][l];
      if (other == NULL_CID) {
        where[o][l] = c;
        ++counts[c];
        horizon[nclustered++] = k;
        vv_pq_remove(k,q);
      } else if (counts[c] < counts[other]) {
        where[o][l] = c;
        ++counts[c];
        --counts[other];
      }
    }

    /* if theres disconnected components, be prepared to jump */
    if (nclustered < graph->nvtxs && nprocessed == nclustered) {
      g = vv_pq_pop(q);
      v = gvtx_to_lvtx(g,graph->dist);
      myid = gvtx_to_tid(g,graph->dist);
      horizon[nclustered++] = g;
      c = cid_rand_r(0,nclusters,&seed);
      ++counts[c];
      where[myid][v] = c;
    }
  }

  vv_pq_free(q);
  dl_free(counts);
  dl_free(horizon);

  clustering = setup_clustering(nclusters,where,graph);
  clustering->free_clusters = 1;
  clustering->free_where = 1;

  return clustering; 
}


/**
 * @brief Grow clusters one at a time to generate an initial clustering.
 *
 * @param nclusters The number of clusters to generate.
 * @param mgraph The graph.
 * @param fixed_clusters Use the number of clusters as a rough target (0), or
 * generate exactly that number of clusters (1).
 * @param seed The random seed to use.
 *
 * @return Return the clustering.
 */
static clustering_t * __cluster_GROW(
    cid_t nclusters, 
    mgraph_t const * const mgraph, 
    const int fixed_clusters, 
    unsigned int seed)
{
  vtx_t v,u,i,k,l,m,maxv,f;
  tid_t o,myid;
  cid_t c;
  adj_t j;
  mod_t priority;
  wgt_t ewgt;
  vtx_t * vtxs, *ptr;
  wgt_t * conn;
  mv_pq_t * q;
  clustering_t * clustering;

  /* expose the graph parts */
  graph_t const * const graph = mgraph->graph;
  vtx_t const nvtxs = graph->nvtxs;
  tid_t const nthreads = graph->npar;
  vtx_t const * const mynvtxs = graph->mynvtxs;
  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy;
  wgt_t const * const * const adjwgt = (wgt_t const * const *)graph->adjwgt;
  wgt_t const * const * const eadjwgt = (wgt_t const * const *)graph->eadjwgt;
  wgt_t const * const * const iadjwgt = (wgt_t const * const *)graph->iadjwgt;

  /* initialize where */
  cid_t ** const where = r_cid_alloc(nthreads);
  for (myid=0;myid<nthreads;++myid) {
    where[myid] = cid_init_alloc(NULL_CID,mynvtxs[myid]);
  }

  q = mv_pq_create(graph->dist.offset,max_gvtx(graph)+1);
  conn = wgt_alloc(max_gvtx(graph)+1);

  /* cluster stuff */
  wgt_t const invgadjwgt = 1.0 / graph->gadjwgt;
  wgt_t const tpwgt = graph->gadjwgt / nclusters;
  wgt_t * pwgts = wgt_calloc(nclusters);

  /* set vtxs to have gvtx numbers */
  vtxs = vtx_alloc(nvtxs);
  maxv = 0;
  for (myid=0;myid<nthreads;++myid) {
    if (mynvtxs[myid] > 0) {
      v = lvtx_to_gvtx(mynvtxs[myid]-1,myid,graph->dist);
      dl_storemax(maxv,v);
    }
  }
  ptr = vtx_alloc(maxv+1);
  i = 0;
  for (myid=0;myid<nthreads;++myid) {
    for (v=0;v<mynvtxs[myid];++v) {
      u = lvtx_to_gvtx(v,myid,graph->dist);
      ptr[u] = i;
      vtxs[i++] = u;
    }
  }

  c = 0;
  m = 0;
  /* grow growths */
  while (m<nvtxs) {
    v = vtxs[vtx_rand_r(m,nvtxs,&seed)];
    m = __move_vertex(v,m,vtxs,ptr);
    i = gvtx_to_lvtx(v,graph->dist);
    myid = gvtx_to_tid(v,graph->dist);
    where[myid][i] = c;
    pwgts[c] += eadjwgt[myid][i] + iadjwgt[myid][i];

    /* add v's neighbors to the queue */
    for (j=xadj[myid][i];j<xadj[myid][i+1];++j) {
      k = adjncy[myid][j];
      if (k < mynvtxs[myid]) {
        l = k;
        o = myid;
      } else {
        l = gvtx_to_lvtx(k,graph->dist);
        o = gvtx_to_tid(k,graph->dist);
      }
      if (where[o][l] == NULL_CID) {
        conn[lvtx_to_gvtx(l,o,graph->dist)] = adjwgt[myid][j];
        priority = __merge_priority(adjwgt[myid][j],pwgts[c], \
            eadjwgt[o][l]+iadjwgt[o][l],invgadjwgt);
        mv_pq_push(priority,lvtx_to_gvtx(l,o,graph->dist),q);
        where[o][l] = ENQUEUED;
      }
    }

    /* add vertices to the cluster unconditionally */
    while (pwgts[c] < 0.5*tpwgt && nvtxs - m > nclusters - c - 1) {
      if (q->size == 0) {
        break;
      }

      /* see if it has an accurate priority */
      f=0;
      while (1) {
        v = mv_pq_peek(q);
        i = gvtx_to_lvtx(v,graph->dist);
        myid = gvtx_to_tid(v,graph->dist);
        priority = __merge_priority(conn[v],pwgts[c], \
            eadjwgt[myid][i]+iadjwgt[myid][i],invgadjwgt);
        if (priority != mv_pq_top(q) && f < MAX_FIXES) {
          mv_pq_update(priority,v,q);
          ++f;
        } else {
          break;
        }
      }
      mv_pq_pop(q);

      m = __move_vertex(v,m,vtxs,ptr);
      where[myid][i] = c;
      pwgts[c] += eadjwgt[myid][i] + iadjwgt[myid][i];

      /* add v's neighbors to the queue */
      for (j=xadj[myid][i];j<xadj[myid][i+1];++j) {
        k = adjncy[myid][j];
        if (k < mynvtxs[myid]) {
          l = k;
          o = myid;
        } else {
          l = gvtx_to_lvtx(k,graph->dist);
          o = gvtx_to_tid(k,graph->dist);
        }
        if (where[o][l] == NULL_CID) {
          conn[lvtx_to_gvtx(l,o,graph->dist)] = adjwgt[myid][j];
          priority = __merge_priority(adjwgt[myid][j],pwgts[c], \
              eadjwgt[o][l]+iadjwgt[o][l],invgadjwgt);
          mv_pq_push(priority,lvtx_to_gvtx(l,o,graph->dist),q);
          where[o][l] = ENQUEUED;
        } else if (where[o][l] == ENQUEUED) {
          ewgt = (conn[lvtx_to_gvtx(l,o,graph->dist)] += adjwgt[myid][j]);
          priority = __merge_priority(ewgt,pwgts[c], \
              eadjwgt[o][l]+iadjwgt[o][l],invgadjwgt);
          mv_pq_update(priority,lvtx_to_gvtx(l,o,graph->dist),q);
        }
      }
    }
    /* add only positive gain vertices */
    while (!fixed_clusters || \
        (pwgts[c] < 2.0*tpwgt && nvtxs - m > nclusters - c - 1)) {
      if (q->size == 0) {
        break;
      }

      /* see if it has an accurate priority */
      f=0;
      while (1) {
        v = mv_pq_peek(q);
        i = gvtx_to_lvtx(v,graph->dist);
        myid = gvtx_to_tid(v,graph->dist);
        priority = __merge_priority(conn[v],pwgts[c], \
            eadjwgt[myid][i]+iadjwgt[myid][i],invgadjwgt);
        if (priority != mv_pq_top(q) && f < MAX_FIXES) {
          mv_pq_update(priority,v,q);
          ++f;
        } else {
          break;
        }
      }
      v = mv_pq_pop(q);

      /* see if it's actually a positive gain move */
      if (priority <= 0) {
        where[myid][i] = NULL_CID; /* reset popped value */
        break;
      }
      m = __move_vertex(v,m,vtxs,ptr);
      where[myid][i] = c;
      pwgts[c] += eadjwgt[myid][i] + iadjwgt[myid][i];

      /* add v's neighbors to the queue */
      for (j=xadj[myid][i];j<xadj[myid][i+1];++j) {
        k = adjncy[myid][j];
        if (k < mynvtxs[myid]) {
          l = k;
          o = myid;
        } else {
          l = gvtx_to_lvtx(k,graph->dist);
          o = gvtx_to_tid(k,graph->dist);
        }
        if (where[o][l] == NULL_CID) {
          conn[lvtx_to_gvtx(l,o,graph->dist)] = adjwgt[myid][j];
          priority = __merge_priority(adjwgt[myid][j],pwgts[c], \
              eadjwgt[o][l]+iadjwgt[o][l],invgadjwgt);
          mv_pq_push(priority,lvtx_to_gvtx(l,o,graph->dist),q);
          where[o][l] = ENQUEUED;
        } else if (where[o][l] == ENQUEUED) {
          ewgt = (conn[lvtx_to_gvtx(l,o,graph->dist)] += adjwgt[myid][j]);
          priority = __merge_priority(ewgt,pwgts[c], \
              eadjwgt[o][l]+iadjwgt[o][l],invgadjwgt);
          mv_pq_update(priority,lvtx_to_gvtx(l,o,graph->dist),q);
        }
      }
    }
    /* empty the priority queue and clear ENQUEUED */
    if (q->size > 0.125 * nvtxs) {
      mv_pq_clear(q);
      for (myid=0;myid<nthreads;++myid) {
        for (i=0;i<mynvtxs[myid];++i) {
          if (where[myid][i] == ENQUEUED) {
            where[myid][i] = NULL_CID;
          }
        }
      }
    } else { 
      while (q->size > 0) {
        v = q->elements[--q->size].val;
        q->index[v-graph->dist.offset] = (size_t)-1;
        i = gvtx_to_lvtx(v,graph->dist);
        myid = gvtx_to_tid(v,graph->dist);
        where[myid][i] = NULL_CID;
      }
    }
    ++c;
    if (c >= nclusters) {
      if (fixed_clusters) {
        c = 0;
      } else {
        nclusters *= 2;
        pwgts = wgt_realloc(pwgts,nclusters);
        wgt_set(pwgts+nclusters/2,0,nclusters/2); /* zero out new range */
      }
    }
  }

  if (!fixed_clusters) {
    nclusters = c;
  }

  dl_free(conn);
  dl_free(vtxs);
  dl_free(ptr);
  dl_free(pwgts);
  mv_pq_free(q);

  clustering = setup_clustering(nclusters,where,graph);
  clustering->free_clusters = 1;
  clustering->free_where = 1;

  return clustering;
}


/**
 * @brief Generate a clustering using label propagation.
 *
 * @param mgraph The graph.
 * @param seed The random seed to use.
 *
 * @return The generate clustering.
 */
static clustering_t * __cluster_LP(
    mgraph_t const * const mgraph, 
    unsigned int seed)
{
  vtx_t i,k,l,p,nrelabels,maxlabel,initlabel,mynvtxs,nlabels;
  adj_t j;
  tid_t t, o;
  clustering_t * clustering;
  cid_t nclusters;
  vtx_t ** label, ** perm;
  cid_t * clbl, ** where;

  graph_t const * const graph = mgraph->graph;
  tid_t const nthreads = graph->npar;
  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy;
  wgt_t const * const * const adjwgt = (wgt_t const * const *)graph->adjwgt;

  nlabels = 0;
  label = r_vtx_alloc(nthreads);
  perm = r_vtx_alloc(nthreads);
  for (t=0;t<nthreads;++t) {
    mynvtxs = graph->mynvtxs[t];
    label[t] = vtx_alloc(mynvtxs);
    vtx_incset(label[t],t,nthreads,mynvtxs);
    if (label[t][mynvtxs-1] > nlabels) {
      nlabels = label[t][mynvtxs-1];
    }
    vtx_pseudo_shuffle_r(label[t],mynvtxs/8,mynvtxs,&seed);
    perm[t] = vtx_alloc(mynvtxs);
    vtx_incset(perm[t],t,nthreads,mynvtxs);
  }
  ++nlabels; /* change from max to number */

  wgt_t * wgts = wgt_calloc(nlabels);

  nrelabels = 1;
  while (nrelabels > 0) {
    nrelabels = 0;
    for (t=0;t<nthreads;++t) {
      mynvtxs = graph->mynvtxs[t];
      vtx_pseudo_shuffle_r(perm[t],mynvtxs/8,mynvtxs,&seed);
      for (p=0;p<mynvtxs;++p) {
        i = perm[t][p];
        initlabel = maxlabel = label[t][i];
        #ifdef LP_IWGT
        wgts[maxlabel] = graph->iadjwgt[t][i] / 2.0;
        #endif
        for (j=xadj[t][i];j<xadj[t][i+1];++j) {
          k = adjncy[t][j];
          if (k < mynvtxs) {
            l = k;
            o = t;
          } else {
            l = gvtx_to_lvtx(k,graph->dist);
            o = gvtx_to_tid(k,graph->dist);
          }
          wgts[label[o][l]] += adjwgt[t][j];
          if (wgts[label[o][l]] > wgts[maxlabel]) {
            maxlabel = label[o][l];
          }
        }
        if (wgts[maxlabel] > wgts[initlabel]) {
          label[t][i] = maxlabel;
          ++nrelabels;
        }
        for (j=xadj[t][i];j<xadj[t][i+1];++j) {
          k = adjncy[t][j];
          if (k < mynvtxs) {
            l = k;
            o = t;
          } else {
            l = gvtx_to_lvtx(k,graph->dist);
            o = gvtx_to_tid(k,graph->dist);
          }
          wgts[label[o][l]] = 0.0;
        }
        wgts[initlabel] = 0.0;
      }
    }
  }
  dl_free(wgts);
  r_vtx_free(perm,nthreads);

  nclusters = 0;
  clbl = cid_init_alloc(NULL_CID,nlabels);
  where = r_cid_alloc(nthreads);
  for (t=0;t<nthreads;++t) {
    mynvtxs = graph->mynvtxs[t];
    where[t] = cid_alloc(mynvtxs);
    for (i=0;i<mynvtxs;++i) {
      if (clbl[label[t][i]] == NULL_CID) {
        clbl[label[t][i]] = nclusters++;
      }
      where[t][i] = clbl[label[t][i]];
    }
    dl_free(label[t]);
  }
  dl_free(label);
  dl_free(clbl);

  clustering = setup_clustering(nclusters,where,graph);
  clustering->free_clusters = 1;
  clustering->free_where = 1;

  return clustering; 
}


/**
 * @brief Generate a clustering using the CNM method.
 *
 * @param mgraph The graph.
 * @param seed The random seed to use.
 *
 * @return The generated clustering.
 */
static clustering_t * __cluster_NEWMAN(
    mgraph_t const * const mgraph, 
    unsigned int seed)
{
  return generate_newman_clustering(mgraph->graph);
}


/**
 * @brief Generate a clustering by growing clusters while using hill-climbing
 * to determine their size.
 *
 * @param mgraph The graph.
 * @param seed The random seed to use.
 *
 * @return The generated clustering.
 */
static clustering_t * __cluster_GROWKL(
    mgraph_t const * const mgraph, 
    unsigned int seed)
{
  vtx_t v,u,i,k,l,maxv,f,bs;
  tid_t o,myid;
  cid_t nclusters;
  adj_t j;
  mod_t priority, maxmod, mod;
  wgt_t pwgt, gwgt, deg;
  wgt_t * conn, * gconn;
  mv_pq_t * q;
  clustering_t * clustering;
  vtx_t movebuffer[HILLSIZE];
  vtx_iset_t * vtxs;

  /* expose the graph parts */
  graph_t const * const graph = mgraph->graph;
  vtx_t const nvtxs = graph->nvtxs;
  tid_t const nthreads = graph->npar;
  vtx_t const * const mynvtxs = graph->mynvtxs;
  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy;
  wgt_t const * const * const adjwgt = (wgt_t const * const *)graph->adjwgt;
  wgt_t const * const * const eadjwgt = (wgt_t const * const *)graph->eadjwgt;
  wgt_t const * const * const iadjwgt = (wgt_t const * const *)graph->iadjwgt;

  /* initialize where */
  cid_t ** const where = r_cid_alloc(nthreads);
  for (myid=0;myid<nthreads;++myid) {
    where[myid] = cid_init_alloc(NULL_CID,mynvtxs[myid]);
  }

  maxv = max_gvtx(graph)+1;

  q = mv_pq_create(graph->dist.offset,maxv);
  conn = wgt_alloc(maxv);
  gconn = wgt_alloc(maxv);


  /* cluster stuff */
  wgt_t const invgadjwgt = 1.0 / graph->gadjwgt;

  vtxs = vtx_iset_create(graph->dist.offset,maxv);

  i = 0;
  for (myid=0;myid<nthreads;++myid) {
    for (v=0;v<mynvtxs[myid];++v) {
      u = lvtx_to_gvtx(v,myid,graph->dist);
      vtx_iset_add(u,vtxs);
      gconn[u] = eadjwgt[myid][v];
    }
  }
  gwgt = graph->gadjwgt;

  dl_timer_t tmr;
  dl_init_timer(&tmr);
  dl_start_timer(&tmr);
  /* grow growths */
  nclusters = 0;
  while (vtxs->size > 0) {
    v = vtx_iset_remove_index(vtx_rand_r(0,vtxs->size,&seed),vtxs);
    i = gvtx_to_lvtx(v,graph->dist);
    myid = gvtx_to_tid(v,graph->dist);
    where[myid][i] = nclusters;
    pwgt = eadjwgt[myid][i] + iadjwgt[myid][i];
    gwgt -= pwgt;
    /* its just relative anyway */
    maxmod = mod = 0; 
    /* add v's neighbors to the queue */
    for (j=xadj[myid][i];j<xadj[myid][i+1];++j) {
      k = adjncy[myid][j];
      if (k < mynvtxs[myid]) {
        l = k;
        o = myid;
        k = lvtx_to_gvtx(l,o,graph->dist);
      } else {
        l = gvtx_to_lvtx(k,graph->dist);
        o = gvtx_to_tid(k,graph->dist);
      }
      if (where[o][l] == NULL_CID) {
        if (iadjwgt) {
          deg = eadjwgt[o][l]+iadjwgt[o][l];
        } else {
          deg = eadjwgt[o][l];
        }
        conn[k] = adjwgt[myid][j];
        gconn[k] -= adjwgt[myid][j];
        priority = __move_priority(gwgt,pwgt,deg,gconn[k],conn[k],invgadjwgt);
        mv_pq_push(priority,k,q);
        where[o][l] = ENQUEUED;
      }
    }

    bs = 0;
    
    /* add vertices to the cluster unconditionally */
    while (vtxs->size > 0) {
      if (q->size == 0) {
        break;
      }

      /* see if it has an accurate priority */
      f=0;
      while (1) {
        v = mv_pq_peek(q);
        i = gvtx_to_lvtx(v,graph->dist);
        myid = gvtx_to_tid(v,graph->dist);
        if (iadjwgt) {
          deg = eadjwgt[myid][i]+iadjwgt[myid][i];
        } else {
          deg = eadjwgt[myid][i];
        }
        priority = __move_priority(gwgt,pwgt,deg,gconn[v],conn[v],invgadjwgt);
        if (priority != mv_pq_top(q) && f < MAX_FIXES) {
          mv_pq_update(priority,v,q);
          ++f;
        } else {
          break;
        }
      }
      mv_pq_pop(q);

      mod += priority;
      if (mod > maxmod) {
        bs = 0;
        maxmod = mod;
      } else if (bs < HILLSIZE) {
        movebuffer[bs++] = v;
      } else {
        where[myid][i] = NULL_CID;
        break;
      }

      vtx_iset_remove(v,vtxs);
      where[myid][i] = nclusters;
      if (iadjwgt) {
        pwgt += eadjwgt[myid][i] + iadjwgt[myid][i];
        gwgt -= eadjwgt[myid][i] + iadjwgt[myid][i];
      } else {
        pwgt += eadjwgt[myid][i];
        gwgt -= eadjwgt[myid][i];
      }

      /* add v's neighbors to the queue */
      for (j=xadj[myid][i];j<xadj[myid][i+1];++j) {
        k = adjncy[myid][j];
        if (k < mynvtxs[myid]) {
          l = k;
          o = myid;
          k = lvtx_to_gvtx(l,o,graph->dist);
        } else {
          l = gvtx_to_lvtx(k,graph->dist);
          o = gvtx_to_tid(k,graph->dist);
        }
        if (iadjwgt) {
          deg = eadjwgt[myid][i]+iadjwgt[myid][i];
        } else {
          deg = eadjwgt[myid][i];
        }
        if (where[o][l] == NULL_CID) {
          gconn[k] -= adjwgt[myid][j];
          conn[k] = adjwgt[myid][j];
          priority = __move_priority(gwgt,pwgt,deg,gconn[k],conn[k], \
              invgadjwgt);
          mv_pq_push(priority,k,q);
          where[o][l] = ENQUEUED;
        } else if (where[o][l] == ENQUEUED) {
          gconn[k] -= adjwgt[myid][j];
          conn[k] += adjwgt[myid][j];
          priority = __move_priority(gwgt,pwgt,deg,gconn[k],conn[k], \
              invgadjwgt);
          mv_pq_update(priority,k,q);
        }
      }
    }

    /* revert to best state */
    while (bs>0) {
      v = movebuffer[--bs];
      i = gvtx_to_lvtx(v,graph->dist);
      myid = gvtx_to_tid(v,graph->dist);
      where[myid][i] = NULL_CID;
      vtx_iset_add(v,vtxs);
      pwgt -= eadjwgt[myid][i] + iadjwgt[myid][i];
      gwgt += eadjwgt[myid][i] + iadjwgt[myid][i];
      /* restore gconn */
      for (j=xadj[myid][i];j<xadj[myid][i+1];++j) {
        k = adjncy[myid][j];
        if (k < mynvtxs[myid]) {
          k = lvtx_to_gvtx(k,myid,graph->dist);
        }
        gconn[k] += adjwgt[myid][j];
      }
    }

    /* empty the priority queue and clear ENQUEUED */
    if (q->size > 0.125 * nvtxs) {
      mv_pq_clear(q);
      for (myid=0;myid<nthreads;++myid) {
        for (i=0;i<mynvtxs[myid];++i) {
          if (where[myid][i] == ENQUEUED) {
            where[myid][i] = NULL_CID;
          }
        }
      }
    } else { 
      while (q->size > 0) {
        v = q->elements[--q->size].val;
        q->index[v-graph->dist.offset] = (size_t)-1;
        i = gvtx_to_lvtx(v,graph->dist);
        myid = gvtx_to_tid(v,graph->dist);
        where[myid][i] = NULL_CID;
      }
    }
    ++nclusters;
  }
  dl_stop_timer(&tmr);

  dl_free(conn);
  vtx_iset_free(vtxs);
  mv_pq_free(q);

  clustering = setup_clustering(nclusters,where,graph);
  clustering->free_clusters = 1;
  clustering->free_where = 1;

  return clustering;
}


/**
 * @brief Generate singleton based clusterings with random boundary
 * refinement.
 *
 * @param objective The objective to supply parameters.
 * @param mgraph The graph.
 *
 * @return The generated clustering.
 */
static clustering_t * __rvtx_clusterings(
    objective_t * const objective,
    mgraph_t * const mgraph)
{
  size_t i,icnum;
  cid_t nclusters = 0;
  mod_t maxmod = -1.0, curmod;
  clustering_t * curclustering, * baseclustering, * clustering = NULL;
  ucinfo_t ** bestucinfo = NULL, ** baseinfo = NULL;

  graph_t const * const graph = mgraph->graph;
  vtx_t const nvtxs = graph->nvtxs;
  adj_t const nedges = graph->nedges;
  tid_t const nthreads = graph->npar;

  /* generate the initial clustering */
  switch (objective->ictype) {
    case NERSTRAND_INITIAL_CLUSTERING_RVTX :
      baseclustering = __cluster_VTX(mgraph,(objective->seed)++);
      break;
    default :
      dl_error("Unknown initial clustering type\n");
  }

  DL_ASSERT(check_clustering(baseclustering,graph) == 1,
      "Bad clustering generated\n");

  /* create the initial ucinfo */
  baseinfo = derive_ucinfo(baseclustering,mgraph);

  for (i=0;i<objective->initsolutions;++i) {
    curclustering = clone_clustering(graph->mynvtxs,baseclustering);
    mgraph->ucinfo = clone_ucinfos(graph->mynvtxs, \
        (ucinfo_t const * const*)baseinfo,nthreads);

    explicit_refine_graph(curclustering,mgraph, \
        NERSTRAND_REFINEMENT_RANDOM,(size_t)-1,1,(objective->seed)++, \
        objective->refstats);

    curmod = calc_total_modularity(curclustering,graph);

    if (objective->icstats) {
      icnum = i+(objective->runnum*objective->initsolutions);
      objective->icstats->lvls[icnum] = mgraph->level;
      objective->icstats->mods[icnum] = curmod;
      objective->icstats->nvtxs[icnum] = nvtxs;
      objective->icstats->nedges[icnum] = nedges;
      objective->icstats->tadjwgt[icnum] = graph->tadjwgt;
      objective->icstats->nclusters[icnum] = curclustering->nclusters;
    }

    vprintf(objective->verbosity,NERSTRAND_VERBOSITY_MAXIMUM,"Created an " \
        "initial solution of "PF_CID_T"/"PF_CID_T" clusters with a " \
        "modularity of "PF_MOD_T" at level "PF_SIZE_T"\n", \
        curclustering->nclusters,nclusters, \
        calc_total_modularity(curclustering,graph),mgraph->level); 

    if (clustering == NULL) {
      clustering = curclustering;
      maxmod = curmod;
      bestucinfo = mgraph->ucinfo;
    } else if (curmod > maxmod) {
      free_clustering(clustering);
      free_ucinfos(bestucinfo,nthreads);
      bestucinfo = mgraph->ucinfo;
      clustering = curclustering;
      maxmod = curmod;
    } else {
      free_clustering(curclustering);
      free_ucinfos(mgraph->ucinfo,nthreads);
      mgraph->ucinfo = bestucinfo;
    }
  }

  free_clustering(baseclustering);
  free_ucinfos(baseinfo,nthreads);

  return clustering;
}


/**
 * @brief Generate non-deterministic clusterings with an unknown number of
 * clusters.
 *
 * @param objective The parameters to use.
 * @param mgraph The graph to cluster.
 *
 * @return The clustering.
 */
static clustering_t * __random_clusterings(
    objective_t * const objective,
    mgraph_t * const mgraph)
{
  size_t i,icnum;
  cid_t nclusters = 0;
  mod_t maxmod = -1.0, curmod;
  clustering_t * curclustering, * clustering = NULL;
  ucinfo_t ** bestucinfo = NULL;

  graph_t const * const graph = mgraph->graph;
  vtx_t const nvtxs = graph->nvtxs;
  adj_t const nedges = graph->nedges;
  tid_t const nthreads = graph->npar;

  for (i=0;i<objective->initsolutions;++i) {
    /* generate the initial clustering */
    switch (objective->ictype) {
      case NERSTRAND_INITIAL_CLUSTERING_LP :
        curclustering = __cluster_LP(mgraph,(objective->seed)++);
        break;
      case NERSTRAND_INITIAL_CLUSTERING_GROWKL :
        curclustering = __cluster_GROWKL(mgraph,(objective->seed)++);
        break;
      default :
        dl_error("Unknown initial clustering type\n");
    }
    
    DL_ASSERT(check_clustering(curclustering,mgraph->graph) == 1, \
        "Bad clustering generated\n");

    /* create the initial ucinfo */
    mgraph->ucinfo = derive_ucinfo(curclustering,mgraph);

    DL_ASSERT(check_ucinfo(mgraph->graph, \
          (ucinfo_t const * const *)mgraph->ucinfo, \
          (cid_t const * const *)curclustering->where, \
          curclustering->nclusters) == 1,  \
        "Bad ucinfo derived\n");
    DL_ASSERT(check_bnd(mgraph, \
         (cid_t const * const *)curclustering->where) == 1, \
        "Bad boundary derived\n");

    curclustering = refine_graph(objective,curclustering,mgraph);

    curmod = calc_total_modularity(curclustering,mgraph->graph);

    if (objective->icstats) {
      icnum = i+(objective->runnum*objective->initsolutions);
      objective->icstats->lvls[icnum] = mgraph->level;
      objective->icstats->mods[icnum] = curmod;
      objective->icstats->nvtxs[icnum] = nvtxs;
      objective->icstats->nedges[icnum] = nedges;
      objective->icstats->tadjwgt[icnum] = mgraph->graph->tadjwgt;
      objective->icstats->nclusters[icnum] = curclustering->nclusters;
    }

    vprintf(objective->verbosity,NERSTRAND_VERBOSITY_MAXIMUM,"Created an " \
        "initial solution of "PF_CID_T"/"PF_CID_T" clusters with a " \
        "modularity of "PF_MOD_T" at level "PF_SIZE_T"\n", \
        curclustering->nclusters,nclusters, \
        calc_total_modularity(curclustering,mgraph->graph),mgraph->level); 

    if (clustering == NULL) {
      clustering = curclustering;
      maxmod = curmod;
      bestucinfo = mgraph->ucinfo;
    } else if (curmod > maxmod) {
      free_clustering(clustering);
      free_ucinfos(bestucinfo,nthreads);
      bestucinfo = mgraph->ucinfo;
      clustering = curclustering;
      maxmod = curmod;
    } else {
      free_clustering(curclustering);
      free_ucinfos(mgraph->ucinfo,nthreads);
      mgraph->ucinfo = bestucinfo;
    }
  }

  return clustering;
}


/**
 * @brief Generate an initial clustering by useing a k-way clustering method
 * makeing several clustering with different values of k.
 *
 * @param objective The parameters to use.
 * @param mgraph The graph to cluster.
 *
 * @return The generated clustering.
 */
static clustering_t * __search_clusterings(
    objective_t * const objective,
    mgraph_t * const mgraph)
{
  size_t i,icnum;
  cid_t nclusters = 0;
  mod_t maxmod = -1.0, curmod;
  clustering_t * curclustering, * clustering = NULL;
  ucinfo_t ** bestucinfo = NULL;

  graph_t const * const graph = mgraph->graph;
  vtx_t const nvtxs = graph->nvtxs;
  adj_t const nedges = graph->nedges;
  tid_t const nthreads = graph->npar;
  real_t const scale = graph->gadjwgt > 0 ? 
      (graph->gadjwgt - graph->tadjwgt) / graph->gadjwgt :
      1.0;
  vtx_t const snvtxs = dl_max(
      dl_min(nvtxs*scale*scale/objective->cnvtxs_per_cluster,nvtxs),
      MIN_SNVTXS);

  for (i=0;i<objective->initsolutions;++i) {
    if (snvtxs > 4) {
      nclusters = cid_rand_r(log2(snvtxs),snvtxs,&(objective->seed));
    } else if (nvtxs > 1) {
      nclusters = cid_rand_r(2,nvtxs,&(objective->seed));
    } else {
      nclusters = 1;
    }

    /* generate the initial clustering */
    switch (objective->ictype) {
      case NERSTRAND_INITIAL_CLUSTERING_RANDOM :
        curclustering = __cluster_RAND(nclusters,mgraph,(objective->seed)++);
        break;
      case NERSTRAND_INITIAL_CLUSTERING_SEED :
        curclustering = __cluster_SEED(nclusters,mgraph,(objective->seed)++);
        break;
      case NERSTRAND_INITIAL_CLUSTERING_GROW :
        curclustering = __cluster_GROW(nclusters,mgraph,0,(objective->seed)++);
        break;
      default :
        dl_error("Unknown initial clustering type\n");
    }
    
    DL_ASSERT(check_clustering(curclustering,mgraph->graph) == 1,
        "Bad clustering generated\n");

    /* create the initial ucinfo */
    mgraph->ucinfo = derive_ucinfo(curclustering,mgraph);

    DL_ASSERT(check_ucinfo(mgraph->graph,
          (ucinfo_t const * const *)mgraph->ucinfo,
          (cid_t const * const *)curclustering->where,
          curclustering->nclusters) == 1, 
        "Bad ucinfo derived\n");
    DL_ASSERT(check_bnd(mgraph,
         (cid_t const * const *)curclustering->where) == 1,
        "Bad boundary derived\n");

    /* initial refinement */
    curclustering = refine_graph(objective,curclustering,mgraph);

    curmod = calc_total_modularity(curclustering,mgraph->graph);

    if (objective->icstats) {
      icnum = i+(objective->runnum*objective->initsolutions);
      objective->icstats->lvls[icnum] = mgraph->level;
      objective->icstats->mods[icnum] = curmod;
      objective->icstats->nvtxs[icnum] = nvtxs;
      objective->icstats->nedges[icnum] = nedges;
      objective->icstats->tadjwgt[icnum] = mgraph->graph->tadjwgt;
      objective->icstats->nclusters[icnum] = curclustering->nclusters;
    }

    vprintf(objective->verbosity,NERSTRAND_VERBOSITY_MAXIMUM,"Created an " \
        "initial solution of "PF_CID_T"/"PF_CID_T" clusters with a " \
        "modularity of "PF_MOD_T" at level "PF_SIZE_T"\n", \
        curclustering->nclusters,nclusters, \
        calc_total_modularity(curclustering,mgraph->graph),mgraph->level); 

    if (clustering == NULL) {
      clustering = curclustering;
      maxmod = curmod;
      bestucinfo = mgraph->ucinfo;
    } else if (curmod > maxmod) {
      free_clustering(clustering);
      free_ucinfos(bestucinfo,nthreads);
      bestucinfo = mgraph->ucinfo;
      clustering = curclustering;
      maxmod = curmod;
    } else {
      free_clustering(curclustering);
      free_ucinfos(mgraph->ucinfo,nthreads);
      mgraph->ucinfo = bestucinfo;
    }
  }

  return clustering;
}




/******************************************************************************
* PUBLIC SERIAL FUNCTIONS *****************************************************
******************************************************************************/


clustering_t * generate_kclustering(
    objective_t * const objective,
    mgraph_t * const mgraph)
{
  dl_start_timer(&(objective->timers.initcluster));

  size_t i, icnum;
  mgraph_t * mymgraph;
  mod_t maxmod = -1.0, curmod;
  clustering_t * clustering = NULL;
  clustering_t * curclustering;
  ucinfo_t ** bestucinfo = NULL;
  dl_timer_t refine_tmr_bkup;

  tid_t const nthreads = mgraph->graph->npar;

  mymgraph = setup_mgraph(mgraph->level,mgraph->graph,NULL,NULL,NULL);

  /* save the state of the refinement timer */
  refine_tmr_bkup = objective->timers.refinement;

  for (i=0;i<objective->initsolutions;++i) {
    /* generate the initial clustering */
    switch (objective->ictype) {
      case NERSTRAND_INITIAL_CLUSTERING_RANDOM :
        curclustering = __cluster_RAND(objective->nclusters,mymgraph, \
            objective->seed++);
        break;
      case NERSTRAND_INITIAL_CLUSTERING_SEED :
        curclustering = __cluster_SEED(objective->nclusters,mymgraph, \
            objective->seed++);
        break;
      case NERSTRAND_INITIAL_CLUSTERING_GROW :
        curclustering = __cluster_GROW(objective->nclusters,mymgraph,1, \
            objective->seed++);
        break;
      default :
        dl_error("Unknown initial clustering type\n");
    }

    DL_ASSERT(check_clustering(curclustering,mymgraph->graph) == 1, \
        "Bad clustering generated\n");

    /* create the initial ucinfo */
    mymgraph->ucinfo = derive_ucinfo(curclustering,mymgraph);

    DL_ASSERT(check_ucinfo(mymgraph->graph, \
          (ucinfo_t const * const *)mymgraph->ucinfo, \
          (cid_t const * const *)curclustering->where, \
          curclustering->nclusters) == 1, \
        "Bad ucinfo derived\n");
    DL_ASSERT(check_bnd(mymgraph, \
         (cid_t const * const *)curclustering->where) == 1, \
        "Bad boundary derived\n");

    /* initial refinement */
    curclustering = refine_graph(objective,curclustering,mymgraph);

    curmod = calc_total_modularity(curclustering,mymgraph->graph);

    if (objective->icstats) {
      icnum = i+(objective->runnum*objective->initsolutions);
      objective->icstats->lvls[icnum] = mgraph->level;
      objective->icstats->mods[icnum] = curmod;
      objective->icstats->nvtxs[icnum] = mymgraph->graph->nvtxs;
      objective->icstats->nedges[icnum] = mymgraph->graph->nedges;
      objective->icstats->tadjwgt[icnum] = mymgraph->graph->tadjwgt;
      objective->icstats->nclusters[icnum] = curclustering->nclusters;
    }

    vprintf(objective->verbosity,NERSTRAND_VERBOSITY_MAXIMUM,"Created an " \
        "initial solution of "PF_CID_T" clusters with a modularity of " \
        PF_MOD_T" at level "PF_SIZE_T"\n",curclustering->nclusters,curmod, \
        mymgraph->level); 

    if (clustering == NULL) {
      clustering = curclustering;
      maxmod = curmod;
      bestucinfo = mymgraph->ucinfo;
    } else if (curmod > maxmod) {
      free_clustering(clustering);
      free_ucinfos(bestucinfo,nthreads);
      bestucinfo = mymgraph->ucinfo;
      clustering = curclustering;
      maxmod = curmod;
    } else {
      free_clustering(curclustering);
      free_ucinfos(mymgraph->ucinfo,nthreads);
      mymgraph->ucinfo = NULL;
    }
  }
  mgraph->ucinfo = bestucinfo;
  mymgraph->free_ucinfo = 0;
  free_mgraph(mymgraph);
  mgraph->free_ucinfo = 1;

  /* restore the refinement timer */
  objective->timers.refinement = refine_tmr_bkup;

  if (objective->icstats) {
    objective->icstats->nbclusters[objective->runnum] = \
        clustering->nclusters;
  }

  vprintf(objective->verbosity,NERSTRAND_VERBOSITY_MEDIUM,"Initial " \
      "clustering with modularity of "PF_MOD_T" selected.\n",maxmod);

  dl_stop_timer(&(objective->timers.initcluster));

  return clustering;
}


clustering_t * generate_aclustering(
    objective_t * const objective, 
    mgraph_t * const mgraph)
{
  mgraph_t * mymgraph;
  clustering_t * clustering = NULL;
  dl_timer_t refine_tmr_bkup;

  graph_t const * const graph = mgraph->graph;

  dl_start_timer(&(objective->timers.initcluster));

  /* save the state of the refinement timer */
  refine_tmr_bkup = objective->timers.refinement;

  if (graph->nedges > MAX_VTX_EDGES) {
    objective->ictype = NERSTRAND_INITIAL_CLUSTERING_GROW;
  }

  mymgraph = setup_mgraph(mgraph->level,graph,NULL,NULL,NULL);
  DL_ASSERT_EQUALS(mymgraph->free_graph,0,"%d");

  switch (objective->ictype) {
    case NERSTRAND_INITIAL_CLUSTERING_RANDOM :
    case NERSTRAND_INITIAL_CLUSTERING_SEED :
    case NERSTRAND_INITIAL_CLUSTERING_GROW :
      clustering = __search_clusterings(objective,mymgraph);
      break;
    case NERSTRAND_INITIAL_CLUSTERING_GROWKL :
    case NERSTRAND_INITIAL_CLUSTERING_LP :
      clustering = __random_clusterings(objective,mymgraph);
      break;
    case NERSTRAND_INITIAL_CLUSTERING_RVTX:
      clustering = __rvtx_clusterings(objective,mymgraph);
      break;
    case NERSTRAND_INITIAL_CLUSTERING_VTX :
      clustering = __cluster_VTX(mymgraph,(objective->seed)++);
      mymgraph->ucinfo = derive_ucinfo(clustering,mymgraph);
      explicit_refine_graph(clustering,mymgraph, \
          NERSTRAND_REFINEMENT_GREEDY,(size_t)-1,1,(objective->seed)++, \
          objective->refstats);
      if (objective->icstats) {
        objective->icstats->lvls[objective->runnum] = mymgraph->level;
        objective->icstats->mods[objective->runnum] = \
            calc_total_modularity(clustering,mgraph->graph);
        objective->icstats->nvtxs[objective->runnum] = graph->nvtxs;
        objective->icstats->nedges[objective->runnum] = graph->nedges;
        objective->icstats->tadjwgt[objective->runnum] = \
            mgraph->graph->tadjwgt;
        objective->icstats->nclusters[objective->runnum] = \
            clustering->nclusters;
      }
      break;
    case NERSTRAND_INITIAL_CLUSTERING_NEWMAN :
      clustering = __cluster_NEWMAN(mymgraph,(objective->seed)++);
      mymgraph->ucinfo = derive_ucinfo(clustering,mymgraph);
      refine_graph(objective,clustering,mymgraph);
      if (objective->icstats) {
        objective->icstats->lvls[objective->runnum] = mymgraph->level;
        objective->icstats->mods[objective->runnum] = \
            calc_total_modularity(clustering,mgraph->graph);
        objective->icstats->nvtxs[objective->runnum] = graph->nvtxs;
        objective->icstats->nedges[objective->runnum] = graph->nedges;
        objective->icstats->tadjwgt[objective->runnum] = \
            mgraph->graph->tadjwgt;
        objective->icstats->nclusters[objective->runnum] = \
            clustering->nclusters;
      }
      break;
    default :
      dl_error("Unknown initial clustering type\n");
  }

  mgraph->ucinfo = mymgraph->ucinfo;
  mymgraph->free_ucinfo = 0;
  free_mgraph(mymgraph);
  mgraph->free_ucinfo = 1;

  vprintf(objective->verbosity,NERSTRAND_VERBOSITY_MEDIUM,"Initial "
      "clustering with modularity of "PF_MOD_T" selected.\n", \
      calc_total_modularity(clustering,mgraph->graph));

  objective->nclusters = clustering->nclusters;
  if (objective->icstats) {
    objective->icstats->nbclusters[objective->runnum] = \
      clustering->nclusters;
  }

  /* restore the refinement timer */
  objective->timers.refinement = refine_tmr_bkup;

  dl_stop_timer(&(objective->timers.initcluster));

  return clustering;
}




/******************************************************************************
* PUBLIC PARALLEL FUNCTIONS ***************************************************
******************************************************************************/


clustering_t * __pgk_clustering = NULL;
clustering_t * par_generate_kclustering(
    objective_t * const objective, 
    mgraph_t * const mgraph)
{
  
  size_t i,icnum,myics,myicstart;
  tid_t bestidx;
  mgraph_t * mymgraph;
  mod_t maxmod = -1.0, curmod;
  clustering_t * myclustering = NULL;
  clustering_t * curclustering;
  ucinfo_t ** bestucinfo = NULL;
  dl_timer_t refine_tmr_bkup;

  tid_t const myid = dlthread_get_id(objective->comm);
  tid_t const nthreads = dlthread_get_nthreads(objective->comm);

  dl_init_timer(&refine_tmr_bkup);

  if (myid == 0) {
    dl_start_timer(&(objective->timers.initcluster));

    /* save the state of the refinement timer */
    refine_tmr_bkup = objective->timers.refinement;
  }
  dlthread_barrier(objective->comm);

  unsigned int myseed = objective->seed+(myid*objective->initsolutions);

  mymgraph = setup_mgraph(mgraph->level,mgraph->graph,NULL,NULL,NULL);
  DL_ASSERT_EQUALS(mymgraph->free_graph,0,"%d");


  myics = (objective->initsolutions/nthreads) + \
      ((objective->initsolutions % nthreads) > (size_t)myid ? 1 : 0); 
  myicstart = (objective->initsolutions/nthreads)*myid + \
      dl_min((size_t)myid,objective->initsolutions%nthreads); 

  for (i=0;i<myics;++i) {
    /* generate the initial clustering */
    switch (objective->ictype) {
      case NERSTRAND_INITIAL_CLUSTERING_RANDOM :
        curclustering = __cluster_RAND(objective->nclusters,mymgraph,myseed++);
        break;
      case NERSTRAND_INITIAL_CLUSTERING_SEED :
        curclustering = __cluster_SEED(objective->nclusters,mymgraph, \
            objective->seed++);
        break;
      case NERSTRAND_INITIAL_CLUSTERING_GROW :
        curclustering = __cluster_GROW(objective->nclusters,mymgraph,1, \
            objective->seed++);
        break;
      default :
        dl_error("Unknown initial clustering type\n");
    }
    
    DL_ASSERT(check_clustering(curclustering,mymgraph->graph) == 1, \
        "Bad clustering generated\n");

    /* create the initial ucinfo */
    mymgraph->ucinfo = derive_ucinfo(curclustering,mymgraph);

    DL_ASSERT(check_ucinfo(mymgraph->graph, \
          (ucinfo_t const * const *)mymgraph->ucinfo, \
          (cid_t const * const *)curclustering->where, \
          curclustering->nclusters) == 1, \
        "Bad ucinfo derived\n");
    DL_ASSERT(check_bnd(mymgraph, \
         (cid_t const * const *)curclustering->where) == 1, \
        "Bad boundary derived\n");

    /* initial refinement */
    curclustering = refine_graph(objective,curclustering,mymgraph);

    curmod = calc_total_modularity(curclustering,mymgraph->graph);

    if (objective->icstats) {
      icnum = i+myicstart+(objective->runnum*objective->initsolutions);
      objective->icstats->lvls[icnum] = mgraph->level;
      objective->icstats->mods[icnum] = curmod;
      objective->icstats->nvtxs[icnum] = mymgraph->graph->nvtxs;
      objective->icstats->nedges[icnum] = mymgraph->graph->nedges;
      objective->icstats->tadjwgt[icnum] = mymgraph->graph->tadjwgt;
      objective->icstats->nclusters[icnum] = curclustering->nclusters;
    }

    vprintf(objective->verbosity,NERSTRAND_VERBOSITY_MAXIMUM,"Created an " \
        "initial solution of "PF_CID_T" clusters with a modularity of " \
        PF_MOD_T" at level "PF_SIZE_T"\n",curclustering->nclusters, \
        calc_total_modularity(curclustering,mymgraph->graph), \
        mymgraph->level); 

    if (myclustering == NULL) {
      myclustering = curclustering;
      maxmod = curmod;
      bestucinfo = mymgraph->ucinfo;
    } else if (curmod > maxmod) {
      free_clustering(myclustering);
      free_ucinfos(bestucinfo,nthreads);
      bestucinfo = mymgraph->ucinfo;
      myclustering = curclustering;
      maxmod = curmod;
    } else {
      free_clustering(curclustering);
      free_ucinfos(mymgraph->ucinfo,nthreads);
      mymgraph->ucinfo = NULL;
    }
  }
  mymgraph->free_ucinfo = 0;
  free_mgraph(mymgraph);

  bestidx = mod_dlthread_maxreduce_index(maxmod,objective->comm);
 
  /* find the best clustering */
  if (bestidx == myid) {
    __pgk_clustering = myclustering;
    mgraph->ucinfo = bestucinfo;
    mgraph->free_ucinfo = 1;
  } else if (myics > 0) {
    free_clustering(myclustering);
    free_ucinfos(bestucinfo,nthreads);
  }

  /* make sure everyone has copied out "seed" */
  if (myid == 0) {
    vprintf(objective->verbosity,NERSTRAND_VERBOSITY_MEDIUM,"Initial " \
        "clustering with modularity of "PF_MOD_T" selected.\n", \
        calc_total_modularity(__pgk_clustering,mgraph->graph));

    if (objective->icstats) {
      objective->icstats->nbclusters[objective->runnum] = \
          __pgk_clustering->nclusters;
    }

    /* restore the state of the refinement timer */
    objective->timers.refinement = refine_tmr_bkup;

    objective->seed = myseed;
    dl_stop_timer(&(objective->timers.initcluster));
  }

  return __pgk_clustering;
}


clustering_t * __pga_clustering = NULL;
clustering_t * par_generate_aclustering(
    objective_t * const objective, 
    mgraph_t * const mgraph)
{
  size_t i,icnum,myics,myicstart;
  tid_t bestidx;
  cid_t nclusters;
  mgraph_t * mymgraph;
  mod_t maxmod = -1.0, curmod;
  clustering_t * myclustering = NULL;
  clustering_t * curclustering;
  ucinfo_t ** bestucinfo = NULL;
  dl_timer_t refine_tmr_bkup;

  tid_t const myid = dlthread_get_id(objective->comm);
  tid_t const nthreads = dlthread_get_nthreads(objective->comm);

  graph_t const * const graph = mgraph->graph;
  vtx_t const nvtxs = graph->nvtxs;
  real_t const scale = graph->gadjwgt > 0 ? \
                       graph->tadjwgt / graph->gadjwgt : 1.0;
  vtx_t const snvtxs = dl_max( \
      dl_min(nvtxs*scale*scale/objective->cnvtxs_per_cluster,nvtxs), \
      MIN_SNVTXS);

  unsigned int myseed = objective->seed+(myid*objective->initsolutions);

  if (myid == 0) {
    dl_start_timer(&(objective->timers.initcluster));

    /* save the state of the refinement timer */
    refine_tmr_bkup = objective->timers.refinement;

    if (graph->nedges > MAX_VTX_EDGES) {
      objective->ictype = NERSTRAND_INITIAL_CLUSTERING_GROW;
    }
  }
  dlthread_barrier(objective->comm);

  mymgraph = setup_mgraph(mgraph->level,graph,NULL,NULL,NULL);
  DL_ASSERT_EQUALS(mymgraph->free_graph,0,"%d");

  myics = (objective->initsolutions/nthreads) + \
      ((objective->initsolutions % nthreads) > (size_t)myid ? 1 : 0); 
  myicstart = (objective->initsolutions/nthreads)*myid + \
      dl_min((size_t)myid,objective->initsolutions%nthreads); 
  for (i=0;i<myics;++i) {
    if (snvtxs > 4) {
      nclusters = cid_rand_r(log2(snvtxs),snvtxs,&myseed);
    } else if (nvtxs > 1) {
      nclusters = cid_rand_r(2,nvtxs,&myseed);
    } else {
      nclusters = 1;
    }

    /* generate the initial clustering */
    switch (objective->ictype) {
      case NERSTRAND_INITIAL_CLUSTERING_RANDOM :
        curclustering = __cluster_RAND(nclusters,mymgraph,myseed++);
        break;
      case NERSTRAND_INITIAL_CLUSTERING_SEED :
        curclustering = __cluster_SEED(nclusters,mymgraph,myseed++);
        break;
      case NERSTRAND_INITIAL_CLUSTERING_GROW :
        curclustering = __cluster_GROW(nclusters,mymgraph,0,myseed++);
        break;
      case NERSTRAND_INITIAL_CLUSTERING_GROWKL :
        curclustering = __cluster_GROWKL(mymgraph,myseed++);
        break;
      case NERSTRAND_INITIAL_CLUSTERING_VTX :
      case NERSTRAND_INITIAL_CLUSTERING_RVTX :
        curclustering = __cluster_VTX(mymgraph,myseed++);
        break;
      default :
        dl_error("Unknown initial clustering type\n");
    }
    
    DL_ASSERT(check_clustering(curclustering,mymgraph->graph) == 1,
        "Bad clustering generated\n");

    /* create the initial ucinfo */
    mymgraph->ucinfo = derive_ucinfo(curclustering,mymgraph);

    DL_ASSERT(check_ucinfo(mymgraph->graph, \
          (ucinfo_t const * const *)mymgraph->ucinfo, \
          (cid_t const * const *)curclustering->where, \
          curclustering->nclusters) == 1,  \
        "Bad ucinfo derived\n"); \
    DL_ASSERT(check_bnd(mymgraph, \
         (cid_t const * const *)curclustering->where) == 1, \
        "Bad boundary derived\n");

    /* initial refinement */
    switch (objective->ictype) {
      case NERSTRAND_INITIAL_CLUSTERING_RVTX:
        explicit_refine_graph(curclustering,mymgraph, \
            NERSTRAND_REFINEMENT_RANDOM,(size_t)-1,1,myseed++, \
            objective->refstats);
        break;
      case NERSTRAND_INITIAL_CLUSTERING_VTX:
        explicit_refine_graph(curclustering,mymgraph, \
            NERSTRAND_REFINEMENT_GREEDY,(size_t)-1,1,myseed++, \
            objective->refstats);
        break;
      default:
        refine_graph(objective,curclustering,mymgraph);
        break;
    }

    curmod = calc_total_modularity(curclustering,mymgraph->graph);

    if (objective->icstats) {
      icnum = i+myicstart+(objective->runnum*objective->initsolutions);
      objective->icstats->lvls[icnum] = mgraph->level;
      objective->icstats->mods[icnum] = curmod;
      objective->icstats->nvtxs[icnum] = mymgraph->graph->nvtxs;
      objective->icstats->nedges[icnum] = mymgraph->graph->nedges;
      objective->icstats->tadjwgt[icnum] = mymgraph->graph->tadjwgt;
      objective->icstats->nclusters[icnum] = curclustering->nclusters;
    }

    vprintf(objective->verbosity,NERSTRAND_VERBOSITY_MAXIMUM,"Created an " \
        "initial solution of "PF_CID_T"/"PF_CID_T" clusters with a " \
        "modularity of "PF_MOD_T" at level "PF_SIZE_T"\n", \
        curclustering->nclusters,nclusters, \
        calc_total_modularity(curclustering,mymgraph->graph),mymgraph->level); 

    if (myclustering == NULL) {
      myclustering = curclustering;
      maxmod = curmod;
      bestucinfo = mymgraph->ucinfo;
    } else if (curmod > maxmod) {
      free_clustering(myclustering);
      free_ucinfos(bestucinfo,nthreads);
      bestucinfo = mymgraph->ucinfo;
      myclustering = curclustering;
      maxmod = curmod;
    } else {
      free_clustering(curclustering);
      free_ucinfos(mymgraph->ucinfo,nthreads);
      mymgraph->ucinfo = NULL;
    }
  }
  mymgraph->free_ucinfo = 0;
  free_mgraph(mymgraph);

  bestidx = mod_dlthread_maxreduce_index(maxmod,objective->comm);
 
  /* find the best clustering */
  if (bestidx == myid) {
    __pga_clustering = myclustering;
    mgraph->ucinfo = bestucinfo;
    mgraph->free_ucinfo = 1;
  } else if (myics > 0) {
    free_clustering(myclustering);
    free_ucinfos(bestucinfo,nthreads);
  }

  /* make sure everyone has copied out "seed" */
  dlthread_barrier(objective->comm);

  if (myid == 0) {
    vprintf(objective->verbosity,NERSTRAND_VERBOSITY_MEDIUM,"Initial " \
        "clustering with modularity of "PF_MOD_T" selected.\n", \
        calc_total_modularity(__pga_clustering,mgraph->graph));

    /* restore the state of the refinement timer */
    objective->timers.refinement = refine_tmr_bkup;

    objective->seed = myseed;
    objective->nclusters = __pga_clustering->nclusters;
    if (objective->icstats) {
      objective->icstats->nbclusters[objective->runnum] = \
          __pga_clustering->nclusters;
    }
    dl_stop_timer(&(objective->timers.initcluster));
  }

  return __pga_clustering;
}




#endif
