/**
 * @file aggregate.c
 * @brief Functions for performing aggregation in serial and parallel
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2014-08-06
 */




#ifndef NERSTRAND_AGGREGATE_C
#define NERSTRAND_AGGREGATE_C




#include "aggregate.h"




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


#define DEF_MASK_SIZE (0x1000)
#define CYCLE_MASK_SIZE (0x100000)
static vtx_t const CYCLE_MASK = CYCLE_MASK_SIZE-1;
static vtx_t const MAX_CYCLE_SIZE = CYCLE_MASK_SIZE >> 1;
static adj_t const MAXDEG2HOP = 3;
static adj_t const MAXSEARCH2HOP = 32;
static vtx_t const DESIRED_VERTEX_SIZE = 4;
static vtx_t const MAX_VERTEX_SIZE = 32;
static wgt_t const REMOTE_PENALTY = 0.75;




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLMEM_PREFIX aggregate
#define DLMEM_TYPE_T aggregate_t
#define DLMEM_DLTYPE DLTYPE_STRUCT
#define DLMEM_INITFUNCTION init_aggregate
#define DLMEM_STATIC
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_INITFUNCTION
#undef DLMEM_DLTYPE
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX


#define DLHT_PREFIX vw
#define DLHT_KEY_T vtx_t
#define DLHT_VAL_T wgt_t
#define DLHT_STATIC
#include "dlht_headers.h"
#undef DLHT_STATIC
#undef DLHT_VAL_T
#undef DLHT_KEY_T
#undef DLHT_PREFIX




/******************************************************************************
* INLINE FUNCTIONS ************************************************************
******************************************************************************/


/**
 * @brief Determine if a fine vertex owns the coarse vertex it is a part of.
 *
 * @param v The global vertex number of mv vertex.
 * @param u The global vertex number of the other vertex.
 * @param dv The degree of my vertex.
 * @param du The degree of the other vertex.
 *
 * @return 1 if the vertex v owns the coarse vertex.
 */
static inline int is_my_cvtx(
    vtx_t const v, 
    vtx_t const u, 
    wgt_t const dv,
    wgt_t const du)
{
  return (v == u) /* If they're the same vertex */
    || dv > du || /* If my vertex has a higher degree */
    (dv == du &&  /* Or break ties evenly */
        ((((v+u)%2) == 0 && v > u) || (((v+u)%2) == 1 && v < u)));
}


/**
 * @brief Determine if enough vertices have been aggregated together.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param collapsed The number of vertices that have been collapsed.
 * @param max_agg_rate The minimum ratio of coarse vertices to fine vertices.
 *
 * @return 1 if aggregation should stop.
 */
static inline int __agg_limit_reached(
    vtx_t const nvtxs, 
    vtx_t const collapsed, 
    real_t const max_agg_rate)
{
  real_t const agg_rate = (real_t)(nvtxs - collapsed);
  return (agg_rate < (real_t)(max_agg_rate*nvtxs));
}


/**
 * @brief Add a vertex i to the same cluster as the vertex maxidx.
 *
 * @param i The vertex to add.
 * @param maxidx The vertex who's cluster to join.
 * @param myid The calling thread's id.
 * @param match The matching vector.
 * @param graph The graph.
 *
 * @return The number of vertices collapsed (1 if formed a cluster, 0 if
 * aggreagted with itself). 
 */
static inline vtx_t __cluster(
    vtx_t const i, 
    vtx_t const maxidx,
    tid_t const myid, 
    vtx_t ** const match, 
    graph_t const * const graph)
{
  vtx_t l, matched;
  tid_t o;
  vtx_t const mynvtxs = graph->mynvtxs[myid];
  if (maxidx == i) {
    /* if we didn't match */
    match[myid][i] = i;
    return 0;
  } else if (maxidx < mynvtxs) {
    /* if I'm matching with something local */
    matched = match[myid][maxidx];
    if (matched == NULL_VTX) {
      match[myid][i] = maxidx;
      match[myid][maxidx] = i;
      return 1;
    } else {
      match[myid][i] = matched;
      match[myid][maxidx] = i;
      return 1;
    }
  } else {
    /* matching with something remote */
    o = gvtx_to_tid(maxidx,graph->dist);
    l = gvtx_to_lvtx(maxidx,graph->dist);
    matched = match[o][l];
    if (matched == NULL_VTX) {
      /* setting up a single match */
      match[myid][i] = maxidx;
      match[o][l] = lvtx_to_gvtx(i,myid,graph->dist);
      return 1;
    } else {
      /* joining a cluster */
      if (matched < graph->mynvtxs[o]) { /* o owns it */
        match[myid][i] = lvtx_to_gvtx(matched,o,graph->dist);
        match[o][l] = lvtx_to_gvtx(i,myid,graph->dist);
      } else if (gvtx_to_tid(matched,graph->dist) == myid) { /* i own it */
        match[myid][i] = gvtx_to_lvtx(matched,graph->dist);
        match[o][l] = lvtx_to_gvtx(i,myid,graph->dist);
      } else { /* third party */
        match[myid][i] = matched;
        match[o][l] = lvtx_to_gvtx(i,myid,graph->dist);
      }
      return 1;
    }
  }
}


/**
 * @brief Match a vertex with itself.
 *
 * @param i The vertex to match with itself.
 * @param cnvtxs The number of coarse vertices.
 * @param myid The id of the calling thread.
 * @param match The calling thread's matching vector.
 * @param cmap The calling thread's coarse vertex map.
 * @param fmap The calling thread's fine vertex map.
 * @param dist The distribution structure of the graph.
 *
 * @return The new number of coarse vertices. 
 */
static inline vtx_t __match_self(
    vtx_t const i, 
    vtx_t const cnvtxs, 
    vtx_t const myid, 
    vtx_t * const match, 
    vtx_t * const cmap, 
    vtx_t * const fmap, 
    graphdist_t const dist)
{
  match[i] = i; 
  fmap[cnvtxs] = i;
  cmap[i] = lvtx_to_gvtx(cnvtxs,myid,dist);

  return cnvtxs+1;
}



/**
 * @brief Find an index in the hash table corresponding to the key.
 *
 * @param i The key.
 * @param htable The hashtable.
 *
 * @return The index corresponding to the key.
 */
static inline vtx_t __htable_idx(
    vtx_t const i, 
    vtx_t const * const htable)
{
  vtx_t idx,m,j;

  idx = i&CYCLE_MASK;
  m = htable[idx];

  if (m == NULL_VTX || m == i) {
    return idx;
  } else {
    for (j=(idx+1)&CYCLE_MASK;;j=(j+1)&CYCLE_MASK) {
      if (htable[j] == NULL_VTX || htable[j] == i) {
        break;
      }
    }
    return j;
  }
}


/**
 * @brief Match the vertex i with the vertex maxidx.
 *
 * @param i The primary vertex to match.
 * @param maxidx The secondary vertex to match.
 * @param myid The id of the calling thread.
 * @param match The thread's matching vector.
 * @param graph The graph.
 *
 * @return The number of vertices collapsed (0 if i was matched with itself 
 *  (i == maxidx), and 1 if i was matched with another vertex).
 */
static inline vtx_t __match(
    vtx_t const i, 
    vtx_t const maxidx, 
    tid_t const myid, 
    vtx_t * const * const match, 
    graph_t const * const graph)
{
  vtx_t l;
  tid_t o;

  vtx_t const mynvtxs = graph->mynvtxs[myid];

  if (maxidx == i) {
    match[myid][i] = i;
    return 0;
  } else if (maxidx < mynvtxs) {
    match[myid][i] = maxidx;
    match[myid][maxidx] = i;
    return 1;
  } else {
    o = gvtx_to_tid(maxidx,graph->dist);
    l = gvtx_to_lvtx(maxidx,graph->dist); 
    match[myid][i] = maxidx;
    match[o][l] = lvtx_to_gvtx(i,myid,graph->dist);
    return 1;
  }
}




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


/**
 * @brief Cleanup broken clusters by removing tail vertices on cluster cycles.
 *
 * @param mynvtxs The number of vertices the calling thread owns.
 * @param myid The id of the calling thread.
 * @param match The global match vector.
 * @param cmap The global coarse vertex map.
 * @param fmap The global fine vertex map.
 * @param graph The graph.
 *
 * @return The number of coarse vertices this thread owns.
 */
static vtx_t __cluster_cleanup(
    vtx_t const mynvtxs,
    tid_t const myid,
    vtx_t ** const match, 
    vtx_t ** const cmap, 
    vtx_t ** const fmap, 
    graph_t const * const graph)
{
  vtx_t i, g, maxidx, idx, l, j, gvtx, ncycle, minvtx, mycnvtxs, cnvtxs;
  tid_t o, maxtid;

  vtx_t * vtxs = vtx_init_alloc(NULL_VTX,CYCLE_MASK_SIZE);
  vtx_t * ptxs = vtx_alloc(MAX_CYCLE_SIZE);

  ncycle = 0;
  mycnvtxs = 0;
  for (i=0;i<mynvtxs;++i) {
    DL_ASSERT_EQUALS(ncycle,(vtx_t)0,PF_VTX_T);
    maxidx = match[myid][i];
    if (maxidx == NULL_VTX) { /* if we didn't get matched */
      mycnvtxs = __match_self(i,mycnvtxs,myid,match[myid],cmap[myid],
          fmap[myid],graph->dist);
    } else if (cmap[myid][i] == NULL_VTX) {
      /* check if the vertex is part of a cycle --
       * if it is, and I own the minimum lvtx (ties go to lower tid's) in that 
       * cycle, assign the cmap and the fmap
       */
      minvtx = i;
      maxtid = myid;
      gvtx = lvtx_to_gvtx(i,myid,graph->dist);
      idx = __htable_idx(gvtx,vtxs);
      vtxs[idx] = gvtx;
      ptxs[ncycle++] = idx;
      o = myid;
      do {
        DL_ASSERT(maxidx<max_gvtx(graph),"Invalid maxidx of "PF_VTX_T"\n",
            maxidx);

        if (maxidx < graph->mynvtxs[o]) {
          g = lvtx_to_gvtx(maxidx,o,graph->dist);
        } else {
          o = gvtx_to_tid(maxidx,graph->dist); 
          g = maxidx;
          maxidx = gvtx_to_lvtx(maxidx,graph->dist);
        }

        DL_ASSERT_EQUALS(gvtx_to_lvtx(g,graph->dist),maxidx,PF_VTX_T);
        DL_ASSERT_EQUALS(gvtx_to_tid(g,graph->dist),o,PF_VTX_T);
        DL_ASSERT_EQUALS(lvtx_to_gvtx(maxidx,o,graph->dist),g,PF_VTX_T);

        idx = __htable_idx(g,vtxs);
        if (vtxs[idx] == NULL_VTX) {
          /* no cycle yet -- will stay in do-while */
          if (maxidx < minvtx || (maxidx == minvtx && o > maxtid)) {
            /* set as owner of the cycle */
            minvtx = maxidx;
            maxtid = o;
          }
          /* record vertex */
          vtxs[idx] = g;
          ptxs[ncycle++] = idx;
          /* jump to next */
          maxidx = match[o][maxidx];
          DL_ASSERT(maxidx < graph->mynvtxs[o] ||
              maxidx > graph->dist.mask,"Invalid match of "PF_VTX_T"/"PF_VTX_T
              " for thread with only "PF_VTX_T" vertices\n",maxidx,
              graph->dist.mask,graph->mynvtxs[o]);
        } else {
          /* found a cycle -- will exit do-while */
          if (g == gvtx) {
            /* found the cycle gvtx is a part of */
            if (maxtid == myid) {
              fmap[myid][mycnvtxs] = i;
              cnvtxs = lvtx_to_gvtx(mycnvtxs,myid,graph->dist);
              for (j=0;j<ncycle;++j) {
                o = gvtx_to_tid(vtxs[ptxs[j]],graph->dist);
                l = gvtx_to_lvtx(vtxs[ptxs[j]],graph->dist);
                DL_ASSERT_EQUALS(cmap[o][l],NULL_VTX,PF_VTX_T);
                cmap[o][l] = cnvtxs;
              }
              ++mycnvtxs;
            }
            break;
          } else {
            /* gvtx is not part of a cycle */
            mycnvtxs = __match_self(i,mycnvtxs,myid,match[myid],cmap[myid],
                fmap[myid],graph->dist);
            break;
          }
        }
        if (ncycle >= MAX_CYCLE_SIZE) {
          /* maximum cycle size reached */
          mycnvtxs = __match_self(i,mycnvtxs,myid,match[myid],cmap[myid],
              fmap[myid],graph->dist);
          break;
        }
      } while (1);
      /* now I need to cleanup my garbage */
      if (ncycle < MAX_CYCLE_SIZE) {
        while (ncycle > 0) {
          vtxs[ptxs[--ncycle]] = NULL_VTX;
        }
      } else {
        vtx_set(vtxs,NULL_VTX,CYCLE_MASK_SIZE);
        ncycle = 0;
      }
    }
  }

  dl_free(vtxs);
  dl_free(ptxs);

  return mycnvtxs;
}


/**
 * @brief Clean up broken matches.
 *
 * @param mynvtxs The number of vertices the calling thread owns.
 * @param myid The id of the calling thread.
 * @param match The global match vector.
 * @param cmap The global coarse mertex map.
 * @param fmap The global fine vertex map.
 * @param graph The graph.
 *
 * @return  The number of coarse vertices the calling thread owns.
 */
static vtx_t __match_cleanup(
    vtx_t const mynvtxs, 
    tid_t const myid,
    vtx_t * const * const match, 
    vtx_t * const * const cmap, 
    vtx_t * const * const fmap, 
    graph_t const * const graph)
{
  vtx_t i, k, v, maxidx, mycnvtxs;
  tid_t o;

  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;

  mycnvtxs = 0;
  for (i=0;i<mynvtxs;++i) {
    maxidx = match[myid][i];
    if (maxidx == NULL_VTX) {
      match[myid][i] = i;
      fmap[myid][mycnvtxs] = i;
      cmap[myid][i] = lvtx_to_gvtx(mycnvtxs++,myid,graph->dist);
    } else if (maxidx == i) {
      fmap[myid][mycnvtxs] = i;
      cmap[myid][i] = lvtx_to_gvtx(mycnvtxs++,myid,graph->dist);
    } else if (maxidx < mynvtxs) { 
      /* a vertex I own */
      if (match[myid][maxidx] == i) { /* a good matching */
        if (i > maxidx) {
          fmap[myid][mycnvtxs] = i;
          cmap[myid][i] = cmap[myid][maxidx] = 
              lvtx_to_gvtx(mycnvtxs++,myid,graph->dist);
        }
      } else { /* a broken matching */
        match[myid][i] = i;
        fmap[myid][mycnvtxs] = i;
        cmap[myid][i] = lvtx_to_gvtx(mycnvtxs++,myid,graph->dist);
      }
    } else {
      o = gvtx_to_tid(maxidx,graph->dist);
      k = gvtx_to_lvtx(maxidx,graph->dist);
      v = lvtx_to_gvtx(i,myid,graph->dist);
      if (match[o][k] == v) { /* a good matching */
        if (is_my_cvtx(v,maxidx,xadj[myid][i+1]-xadj[myid][i],
              xadj[o][k+1] - xadj[o][k])) {
          fmap[myid][mycnvtxs] = i;
          cmap[myid][i] = cmap[o][k] =
              lvtx_to_gvtx(mycnvtxs++,myid,graph->dist);
        }
      } else { /* a broken matching */
        match[myid][i] = i;
        fmap[myid][mycnvtxs] = i;
        cmap[myid][i] = lvtx_to_gvtx(mycnvtxs++,myid,graph->dist);
      }
    }
  }
  return mycnvtxs;
}




/******************************************************************************
* PRIVATE SERIAL AGGREGATION FUNCTIONS ****************************************
******************************************************************************/


/**
 * @brief Update the aggregation statistics struct associate with the passed in 
 * objective for the current aggregation step.
 *
 * @param objective The object to update the stats of.
 * @param graph The graph.
 * @param agg The aggregation structure.
 */
static void __calc_aggstats(
    objective_t * const objective, 
    graph_t const * const graph, 
    aggregate_t const * const agg) 
{
  vtx_t i,k,v,j,mynvtxs,tcnvtxs;
  tid_t myid, o;
  aggregation_stat_t aggstat;

  tid_t const nthreads = graph->npar;

  /* aggregation stat tracking */
  init_aggregation_stat(&aggstat);
  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    for (i=0;i<mynvtxs;++i) {
      if (is_supernode(i,myid,graph,objective->snratio)) {
        ++(aggstat.nsupernodes);
      }
    }
  }

  j = 0;
  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    for (i=0;i<mynvtxs;++i) {
      if (agg->match[myid][i] == i) {
        ++(aggstat.nlonenodes);
      } else {
        v = 0;
        o = myid;
        k = i;
        do {
          ++v;
          k = agg->match[o][k];
          if (k >= graph->mynvtxs[o]) {
            o = gvtx_to_tid(k,graph->dist);
            k = gvtx_to_lvtx(k,graph->dist);
          }          
        } while (o != myid || k != i);
        dl_storemax(aggstat.maxnodesize,v);
        j += v;
      }
    }
  }
  tcnvtxs = vtx_sum(agg->nvtxs,nthreads);
  aggstat.avgnodesize = j / (real_t)(tcnvtxs-aggstat.nlonenodes);
  aggstat.coarsenrate = graph->nvtxs / (real_t)tcnvtxs;
  aggregation_stats_push(aggstat,objective->aggstats);
}


/**
 * @brief Create an aggregation using random matching.
 *
 * @param objective The objective containing aggregation parameters
 * @param mgraph The mgraph to aggregate 
 *
 * @return The new aggregation
 */
static aggregate_t * __aggregate_RM(
    objective_t * const objective, 
    mgraph_t const * const mgraph)
{
  vtx_t i,v,k,l,maxidx,collapsed,mynvtxs,last_unmatched;
  adj_t j;
  tid_t myid, o;

  adj_t const * myxadj;
  vtx_t const * myadjncy;

  /* make accessing my graph easy */
  graph_t const * const graph = mgraph->graph;
  tid_t const nthreads = graph->npar;
  vtx_t const nvtxs = graph->nvtxs; 
  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy;

  /* setup my aggregation information */
  vtx_t ** const match = r_vtx_dalloc(graph->mynvtxs,sizeof(vtx_t),nthreads);

  /* permutation array */
  vtx_t * perm = vtx_alloc(nvtxs);

  dprintf("aggregate_RM() called with a graph containing "PF_VTX_T
      " vertices\n",mgraph->graph->nvtxs);

  for (myid=0;myid<nthreads;++myid) {
    vtx_set(match[myid],NULL_VTX,graph->mynvtxs[myid]);
  }
 
  /* generate matches */
  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];

    myxadj = xadj[myid];
    myadjncy = adjncy[myid];

    vtx_incset(perm,0,1,mynvtxs);
    vtx_shuffle_r(perm,mynvtxs,&objective->seed);

    last_unmatched = collapsed = 0; 
    for (v=0;v<mynvtxs;++v) {
      if (__agg_limit_reached(mynvtxs,collapsed,objective->max_agg_rate)) {
        break;
      }
      i = perm[v];
      if (match[myid][i] == NULL_VTX) {
        maxidx = i;
        if (myxadj[i] < myxadj[i+1]) {
          for (j=myxadj[i];j<myxadj[i+1];++j) {
            k = myadjncy[j];
            if (k < mynvtxs) {
              if (match[myid][k] == NULL_VTX) {
                maxidx = k;
                break;
              }
            } else {
              o = gvtx_to_tid(k,graph->dist);
              l = gvtx_to_lvtx(k,graph->dist);
              if (match[o][l] == NULL_VTX) {
                maxidx = k;
                break;
              }
            }
          }
        } else if (objective->parttype == NERSTRAND_PARTITION_KWAY) { 
          for (;last_unmatched < mynvtxs;++last_unmatched) {
            if (match[myid][perm[last_unmatched]] == NULL_VTX) {
              maxidx = perm[last_unmatched];
              break;
            }
          }
        }
        collapsed += __match(i,maxidx,myid,match,graph);
      }
    }
  }

  dl_free(perm);

  vtx_t * const cnvtxs = vtx_alloc(nthreads);
  vtx_t ** const cmap = r_vtx_dalloc(graph->mynvtxs,sizeof(vtx_t),nthreads);
  vtx_t ** const fmap = r_vtx_dalloc(graph->mynvtxs,sizeof(vtx_t),nthreads);

  /* generate cmap and fmap */
  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    cnvtxs[myid] = __match_cleanup(mynvtxs,myid,match,cmap,fmap,graph);
  }

  aggregate_t * const agg = setup_aggregate(graph->mynvtxs,cnvtxs,match,cmap,
      fmap,nthreads);

  /* aggregation stat tracking */
  if (objective->aggstats) {
    __calc_aggstats(objective,graph,agg);
  }

  return agg;
}


/**
 * @brief Create an aggregation using random clustering.
 *
 * @param objective The objective containing aggregation parameters
 * @param mgraph The mgraph to aggregate 
 *
 * @return The new aggregation
 */
static aggregate_t * __aggregate_RC(
    objective_t * const objective, 
    mgraph_t const * const mgraph)
{
  dprintf("aggregate_RC() called with a graph containing "PF_VTX_T
      " vertices\n",mgraph->graph->nvtxs);

  vtx_t i,v,k,l,maxidx,collapsed,mynvtxs,mycnvtxs,cg;
  adj_t j;
  tid_t myid,o;
  adj_t const * myxadj;
  vtx_t const * myadjncy;

  /* make accessing my graph easy */
  graph_t const * const graph = mgraph->graph;
  tid_t const nthreads = graph->npar;
  vtx_t const nvtxs = graph->nvtxs; 
  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy;
  vtx_t * const cnvtxs = vtx_calloc(nthreads);

  /* setup my aggregation information */
  vtx_t ** const match = r_vtx_dalloc(graph->mynvtxs,sizeof(vtx_t),nthreads);
  vtx_t ** const cmap = r_vtx_dalloc(graph->mynvtxs,sizeof(vtx_t),nthreads);
  vtx_t ** const fmap = r_vtx_dalloc(graph->mynvtxs,sizeof(vtx_t),nthreads);
  
  /* permutation array */
  vtx_t * perm = vtx_alloc(nvtxs);

  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    myxadj = xadj[myid];
    myadjncy = adjncy[myid];

    vtx_incset(perm,0,1,mynvtxs);
    vtx_shuffle_r(perm,mynvtxs,&objective->seed);
    vtx_set(match[myid],NULL_VTX,mynvtxs);

    collapsed = 0;

    mycnvtxs = 0;
    /* match everyone possible on first pass */
    for (v=0;v<mynvtxs;++v) {
      if (__agg_limit_reached(mynvtxs,collapsed,objective->max_agg_rate)) {
        break;
      }
      i = perm[v];
      if (match[myid][i] == NULL_VTX) {
        maxidx = i;
        /* try and find an umatched vertex first */
        if (myxadj[i] < myxadj[i+1]) {
          for (j=myxadj[i];j<myxadj[i+1];++j) {
            k = myadjncy[j];
            if (k < mynvtxs) {
              /* a vertex I own */
              if (match[myid][k] == NULL_VTX) {
                maxidx = k;
                break;
              }
            } else {
              /* vertex someone else owns */
              o = gvtx_to_tid(k,graph->dist);
              l = gvtx_to_lvtx(k,graph->dist);
              if (match[o][l] == NULL_VTX) {
                maxidx = k;
                break;
              }
            }
          }
          /* if I didn't find a vertex an unmatched vertex, just pick my first 
           * neighbor */
          if (j == myxadj[i+1] && j > myxadj[i]) {
            maxidx = myadjncy[myxadj[i]];
          }
        } else if (objective->parttype == NERSTRAND_PARTITION_KWAY &&
            v < mynvtxs-1) {
          maxidx = perm[v+1];
        }
        /* match maker make me a match */
        if (maxidx == i) {
          /* match with self */
          fmap[myid][mycnvtxs] = i;
          cmap[myid][i] = lvtx_to_gvtx(mycnvtxs,myid,graph->dist);
          DL_ASSERT(cmap[myid][i]>=graph->dist.offset,"Generated global "
              "coarse vertex number is smaller than graph offset (gvtx = "
              PF_VTX_T", offset = "PF_VTX_T"\n",v,graph->dist.offset);
          ++mycnvtxs;
        } else {
          if (maxidx < mynvtxs) {
            o = myid;
            l = maxidx;
          } else {
            o = gvtx_to_tid(maxidx,graph->dist);
            l = gvtx_to_lvtx(maxidx,graph->dist);
          }
          if (match[o][l] == NULL_VTX) {
            /* match with unmatched foreign vertex */
            fmap[myid][mycnvtxs] = i;
            cg = lvtx_to_gvtx(mycnvtxs++,myid,graph->dist);
            DL_ASSERT(cg>=graph->dist.offset,"Generated global coarse vertex "
                "number is smaller than graph offset (gvtx = "PF_VTX_T
                ", offset = "PF_VTX_T"\n",cg,graph->dist.offset);
            cmap[myid][i] = cmap[o][l] = cg;
          } else {
            /* match with already matched foreign vertex */
            cg = cmap[myid][i] = cmap[o][l];
          }
        }
        /* assign match vector */
        collapsed += __cluster(i,maxidx,myid,match,graph);
      }
    }
    for (;v<mynvtxs;++v) {
      i = perm[v];
      if (match[myid][i] == NULL_VTX) {
        fmap[myid][mycnvtxs] = i;
        cmap[myid][i] = lvtx_to_gvtx(mycnvtxs,myid,graph->dist);
        match[myid][i] = i;
        DL_ASSERT(cmap[myid][i]>=graph->dist.offset,"Generated global "
            "coarse vertex number is smaller than graph offset (gvtx = "
            PF_VTX_T", offset = "PF_VTX_T"\n",v,graph->dist.offset);
        ++mycnvtxs;
      }
    }
    cnvtxs[myid] = mycnvtxs;
  }

  dl_free(perm);

  aggregate_t * const agg = setup_aggregate(graph->mynvtxs,cnvtxs,match,cmap,
      fmap,nthreads);

  /* aggregation stat tracking */
  if (objective->aggstats) {
    __calc_aggstats(objective,graph,agg);
  }

  return agg;
}


/**
 * @brief Create an aggregation using agglomerative modularity based matching.
 *
 * @param objective The objective containing aggregation parameters
 * @param mgraph The mgraph to aggregate 
 *
 * @return The new aggregation
 */
static aggregate_t * __aggregate_AGM(
    objective_t * const objective, 
    mgraph_t const * const mgraph)
{
  vtx_t v,i,k,l,maxidx,mynvtxs,collapsed,last_unmatched;
  adj_t j;
  wgt_t xdeg,ydeg,cwgt,xdegm;
  tid_t myid, o;

  adj_t const * myxadj;
  vtx_t const * myadjncy;

  /* set up graph stuff */
  graph_t const * graph = mgraph->graph;
  vtx_t const nvtxs = graph->nvtxs;
  tid_t const nthreads = graph->npar;
  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy;
  wgt_t const * const * const adjwgt = (wgt_t const * const *)graph->adjwgt;
  wgt_t const * const * const iadjwgt = (wgt_t const * const *)graph->iadjwgt;
  wgt_t const * const * const eadjwgt = (wgt_t const * const *)graph->eadjwgt;

  /* setup my aggregation information */
  vtx_t ** const match = r_vtx_dalloc(graph->mynvtxs,sizeof(vtx_t),nthreads);

  wgt_t const mscale = objective->degwgt/graph->gadjwgt;

  real_t cmod, maxmod;

  vtx_t * perm = vtx_alloc(nvtxs);

  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    myxadj = xadj[myid];
    myadjncy = adjncy[myid];

    vtx_incset(perm,0,1,mynvtxs);
    vtx_shuffle_r(perm,mynvtxs,&objective->seed);

    vtx_set(match[myid],NULL_VTX,mynvtxs);

    last_unmatched = collapsed = 0;
    for (v=0;v<mynvtxs;++v) {
      if (__agg_limit_reached(mynvtxs,collapsed,objective->max_agg_rate)) {
        break;
      }
      i = perm[v];
      if (match[myid][i] == NULL_VTX) {
        if (iadjwgt) {
          xdeg = eadjwgt[myid][i] + iadjwgt[myid][i];
        } else {
          xdeg = eadjwgt[myid][i];
        }
        xdegm = xdeg*mscale;
        if (objective->parttype == NERSTRAND_PARTITION_ANYWAY) {
          maxmod =  0;
        } else {
          maxmod = -graph->gadjwgt;
        }
        maxidx = i;
        if (myxadj[i] < myxadj[i+1]) {
          for (j=myxadj[i];j<myxadj[i+1];++j) {
            k = myadjncy[j];
            if (adjwgt) {
              cwgt = adjwgt[myid][j];
            } else {
              cwgt = 1.0;
            }
            if (k < mynvtxs) {
              o = myid;
              l = k;
            } else {
              o = gvtx_to_tid(k,graph->dist);
              l = gvtx_to_lvtx(k,graph->dist);
            }
            if (iadjwgt) {
              ydeg = iadjwgt[o][l] + eadjwgt[o][l];
            } else {
              ydeg = eadjwgt[o][l];
            }
            cmod = cwgt - (xdegm*ydeg);
            if (match[o][l] == NULL_VTX && maxmod < cmod) {
              maxidx = k;
              maxmod = cmod;
            }
          }
        } else if (objective->parttype == NERSTRAND_PARTITION_KWAY) { 
          for (;last_unmatched < mynvtxs;++last_unmatched) {
            if (match[myid][perm[last_unmatched]] == NULL_VTX) {
              maxidx = perm[last_unmatched];
              break;
            }
          }
        }
        collapsed += __match(i,maxidx,myid,match,graph);
      }
    }
  }

  dl_free(perm);

  vtx_t * const cnvtxs = vtx_alloc(nthreads);
  vtx_t ** const cmap = r_vtx_dalloc(graph->mynvtxs,sizeof(vtx_t),nthreads);
  vtx_t ** const fmap = r_vtx_dalloc(graph->mynvtxs,sizeof(vtx_t),nthreads);

  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    cnvtxs[myid] = __match_cleanup(mynvtxs,myid,match,cmap,fmap,graph);
  }

  aggregate_t * const agg = setup_aggregate(graph->mynvtxs,cnvtxs,match,cmap,
      fmap,nthreads);

  /* aggregation stat tracking */
  if (objective->aggstats) {
    __calc_aggstats(objective,graph,agg);
  }

  return agg;
}


/**
 * @brief Create an aggregation using agglomerative modularity based matching.
 *
 * @param objective The objective containing aggregation parameters
 * @param mgraph The mgraph to aggregate 
 *
 * @return The new aggregation
 */
static aggregate_t * __aggregate_AGH(
    objective_t * const objective, 
    mgraph_t const * const mgraph)
{
  vtx_t v,i,k,kk,l,ll,n,maxidx,mynvtxs,collapsed,mysize;
  adj_t j,jj;
  wgt_t xdeg,ydeg,cwgt,xdegm;
  tid_t myid,o,oo;

  adj_t const * myxadj;
  vtx_t const * myadjncy;

  /* set up graph stuff */
  graph_t const * graph = mgraph->graph;
  vtx_t const nvtxs = graph->nvtxs;
  tid_t const nthreads = graph->npar;
  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy;
  wgt_t const * const * const adjwgt = (wgt_t const * const *)graph->adjwgt;
  wgt_t const * const * const iadjwgt = (wgt_t const * const *)graph->iadjwgt;
  wgt_t const * const * const eadjwgt = (wgt_t const * const *)graph->eadjwgt;

  /* setup my aggregation information */
  vtx_t ** const match = r_vtx_dalloc(graph->mynvtxs,sizeof(vtx_t),nthreads);

  wgt_t const mscale = objective->degwgt/graph->gadjwgt;

  real_t cmod, maxmod;

  vtx_t * perm = vtx_alloc(nvtxs);

  for (myid=0;myid<nthreads;++myid) {
    vtx_set(match[myid],NULL_VTX,graph->mynvtxs[myid]);
  }

  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    myxadj = xadj[myid];
    myadjncy = adjncy[myid];

    vtx_incset(perm,0,1,mynvtxs);
    vtx_shuffle_r(perm,mynvtxs,&objective->seed);

    collapsed = 0;
    for (v=0;v<mynvtxs;++v) {
      if (__agg_limit_reached(mynvtxs,collapsed,objective->max_agg_rate)) {
        break;
      }
      i = perm[v];
      if (match[myid][i] == NULL_VTX) {
        if (iadjwgt) {
          xdeg = eadjwgt[myid][i] + iadjwgt[myid][i];
        } else {
          xdeg = eadjwgt[myid][i];
        }
        xdegm = xdeg*mscale;
        if (objective->parttype == NERSTRAND_PARTITION_ANYWAY) {
          maxmod =  0;
        } else {
          maxmod = -graph->gadjwgt;
        }
        maxidx = i;
        mysize = myxadj[i+1] - myxadj[i];
        for (j=myxadj[i];j<myxadj[i+1];++j) {
          k = myadjncy[j];
          if (adjwgt) {
            cwgt = adjwgt[myid][j];
          } else {
            cwgt = 1.0;
          }
          if (k < mynvtxs) {
            o = myid;
            l = k;
          } else {
            o = gvtx_to_tid(k,graph->dist);
            l = gvtx_to_lvtx(k,graph->dist);
          }
          if (iadjwgt) {
            ydeg = iadjwgt[o][l] + eadjwgt[o][l];
          } else {
            ydeg = eadjwgt[o][l];
          }
          cmod = cwgt - (xdegm*ydeg);
          if (match[o][l] == NULL_VTX && maxmod < cmod) {
            maxidx = k;
            maxmod = cmod;
          }
        }
        /* two hop */
        if (maxidx == i && mysize <= MAXDEG2HOP) {
          for (j=myxadj[i];j<myxadj[i+1];++j) {
            k = myadjncy[j];
            if (k < mynvtxs) {
              o = myid;
              l = k;
            } else {
              o = gvtx_to_tid(k,graph->dist);
              l = gvtx_to_lvtx(k,graph->dist);
            }
            n = 0;
            jj = xadj[o][l] + (i % (xadj[o][l+1]-xadj[o][l]));
            for (n=0;n<MAXSEARCH2HOP&&n<xadj[o][l+1]-xadj[o][l];++n,++jj) {
              if (jj == xadj[o][l+1]) {
                jj = xadj[o][l];
              }
              kk = adjncy[o][jj];
              if (kk < graph->mynvtxs[o]) {
                oo = o;
                ll = kk;
                kk = lvtx_to_gvtx(ll,oo,graph->dist);
              } else {
                oo = gvtx_to_tid(kk,graph->dist);
                ll = gvtx_to_lvtx(kk,graph->dist);
              }
              if (match[oo][ll] == NULL_VTX && !(ll == i && oo == myid) &&
                  xadj[oo][ll+1] - xadj[oo][ll] <= MAXDEG2HOP) {
                if (oo == myid) {
                  maxidx = gvtx_to_lvtx(kk,graph->dist);
                } else {
                  maxidx = kk;
                }
                goto ENDTWOHOP;
              }
            }
          }
          ENDTWOHOP:;
        }

        collapsed += __match(i,maxidx,myid,match,graph);
      }
    }
  }

  dl_free(perm);

  vtx_t * const cnvtxs = vtx_alloc(nthreads);
  vtx_t ** const cmap = r_vtx_dalloc(graph->mynvtxs,sizeof(vtx_t),nthreads);
  vtx_t ** const fmap = r_vtx_dalloc(graph->mynvtxs,sizeof(vtx_t),nthreads);

  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    cnvtxs[myid] = __match_cleanup(mynvtxs,myid,match,cmap,fmap,graph);
  }

  aggregate_t * const agg = setup_aggregate(graph->mynvtxs,cnvtxs,match,cmap,
      fmap,nthreads);

  /* aggregation stat tracking */
  if (objective->aggstats) {
    __calc_aggstats(objective,graph,agg);
  }

  return agg;
}


/**
 * @brief Create an aggregation using agglomerative modularity based
 * clustering.
 *
 * @param objective The objective containing aggregation parameters
 * @param mgraph The mgraph to aggregate 
 *
 * @return The new aggregation
 */
static aggregate_t * __aggregate_AGC(
    objective_t * const objective, 
    mgraph_t const * const mgraph)
{
  vtx_t v,i,k,l,cl,cg,maxidx,collapsed,mynvtxs,mycnvtxs;
  tid_t myid,o,co;
  adj_t j;
  mod_t cmod, maxmod;
  wgt_t xdeg, ydeg, cwgt, xdegm, pen;

  adj_t const * myxadj;
  vtx_t const * myadjncy;

  /* set up graph stuff */
  graph_t const * graph = mgraph->graph;
  vtx_t const nvtxs = graph->nvtxs;
  tid_t const nthreads = graph->npar;
  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy;
  wgt_t const * const * const adjwgt = (wgt_t const * const *)graph->adjwgt;
  wgt_t const * const * const iadjwgt = (wgt_t const * const *)graph->iadjwgt;
  wgt_t const * const * const eadjwgt = (wgt_t const * const *)graph->eadjwgt;

  wgt_t const mscale = objective->degwgt/graph->gadjwgt;

  aggregate_t * const agg = setup_aggregate(graph->mynvtxs,NULL,NULL,NULL,NULL,
      nthreads);

  agg->free_match = 0;
  agg->free_cmap = 0;
  agg->free_fmap = 0;
  agg->free_nvtxs = 0;

  /* setup my aggregation information */
  vtx_t ** const match = agg->match;
  vtx_t ** const cmap = agg->cmap;
  vtx_t ** const fmap = agg->fmap;

  vtx_t * perm = vtx_alloc(nvtxs);

  /* store total degree of coarse vertices */
  wgt_t ** cdeg = r_wgt_dalloc(graph->mynvtxs,sizeof(vtx_t),nthreads);
  vtx_t ** csize = r_vtx_dalloc(graph->mynvtxs,sizeof(vtx_t),nthreads);
  vw_ht_t * conn;

  /* do the real work */
  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    vtx_set(match[myid],NULL_VTX,mynvtxs);
  }

  for (myid=0;myid<nthreads;++myid) {
    mycnvtxs = 0;
    conn = vw_ht_create(graph->maxdeg[myid],graph->maxdeg[myid]/4);
    collapsed = 0;
    mynvtxs = graph->mynvtxs[myid];
    myxadj = xadj[myid];
    myadjncy = adjncy[myid];

    /* counting sort stuff */
    vtx_incset(perm,0,1,mynvtxs);
    vtx_shuffle_r(perm,mynvtxs,&objective->seed);

    for (v=0;v<mynvtxs;++v) {
      if (__agg_limit_reached(mynvtxs,collapsed,objective->max_agg_rate)) {
        break;
      }
      i = perm[v];
      if (match[myid][i] == NULL_VTX) {
        if (iadjwgt) {
          xdeg = iadjwgt[myid][i] + eadjwgt[myid][i]; 
        } else {
          xdeg = eadjwgt[myid][i];
        }
        xdegm = xdeg*mscale;
        if (objective->parttype == NERSTRAND_PARTITION_ANYWAY) {
          maxmod =  0;
        } else {
          maxmod = -graph->gadjwgt;
        }
        maxidx = i;
        if (myxadj[i] == myxadj[i+1]) {
          if (objective->parttype == NERSTRAND_PARTITION_KWAY && \
              v < mynvtxs-1 && \
              objective->max_agg_rate*graph->nvtxs > 2*objective->cnvtxs) {
            maxidx = perm[v+1];
          } 
        } else if (myxadj[i] == myxadj[i+1]+1) {
          /* if we only have one option, always match it. Brandes et al. '08
           * "On Modularity", showed that 1-degree vertices should never be in
           * their own cluster. */
          maxidx = myadjncy[myxadj[i]];
        } else if (myxadj[i] < myxadj[i+1]) {
          /* hopefully for many graphs this will fit in cache */
          for (j=myxadj[i];j<myxadj[i+1];++j) {
            k = myadjncy[j];
            if (k < mynvtxs) {
              o = myid;
              l = k;
            } else {
              o = gvtx_to_tid(k,graph->dist);
              l = gvtx_to_lvtx(k,graph->dist);
            }
            if (match[o][l] != NULL_VTX) {
              if (adjwgt) {
                vw_ht_add(cmap[o][l],adjwgt[myid][j],conn);
              } else {
                vw_ht_add(cmap[o][l],1.0,conn);
              }
            }
          }
          for (j=myxadj[i];j<myxadj[i+1];++j) {
            k = myadjncy[j];
            if (k < mynvtxs) {
              o = myid;
              l = k;
            } else {
              o = gvtx_to_tid(k,graph->dist);
              l = gvtx_to_lvtx(k,graph->dist);
            }
            if (match[o][l] == NULL_VTX) {
              if (iadjwgt) {
                ydeg = iadjwgt[o][l] + eadjwgt[o][l];
              } else {
                ydeg = eadjwgt[o][l];
              }
              if (adjwgt) {
                cwgt = adjwgt[myid][j];
              } else {
                cwgt = 1.0;
              }
              pen = 1;
            } else {
              cl = gvtx_to_lvtx(cmap[o][l],graph->dist);
              co = gvtx_to_tid(cmap[o][l],graph->dist);
              ydeg = cdeg[co][cl];
              if (csize[co][cl] >= MAX_CYCLE_SIZE) {
                continue;
              } else if (csize[co][cl] > MAX_VERTEX_SIZE) {
                pen = DESIRED_VERTEX_SIZE / (wgt_t)csize[co][cl];
              } else if (csize[co][cl] > DESIRED_VERTEX_SIZE) {
                pen = ((0.5*DESIRED_VERTEX_SIZE) / csize[co][cl]) + 0.5;
              } else { 
                pen = 1;
              }
              cwgt = vw_ht_get(cmap[o][l],conn);
            }
            cmod = (cwgt - (ydeg*xdegm))*pen; 
            if (maxidx > mynvtxs) {
              cmod *= REMOTE_PENALTY;
            }
            if (maxmod <= cmod) {
              maxidx = k;
              maxmod = cmod;
            }
          }
          /* clear the conn map */
          vw_ht_clear_chains(conn);
          for (j=myxadj[i];j<myxadj[i+1];++j) {
            k = myadjncy[j];
            if (k < mynvtxs) {
              o = myid;
              l = k;
            } else {
              o = gvtx_to_tid(k,graph->dist);
              l = gvtx_to_lvtx(k,graph->dist);
            }
            if (match[o][l] != NULL_VTX) {
              vw_ht_clear_slot(cmap[o][l],conn);
            }
          }
        }
        /* match maker make me a match */
        if (maxidx == i) {
          /* match with self */
          cdeg[myid][mycnvtxs] = xdeg;
          csize[myid][mycnvtxs] = 1;
          //fmap[myid][mycnvtxs] = i;
          cmap[myid][i] = lvtx_to_gvtx(mycnvtxs,myid,graph->dist);
          DL_ASSERT(cmap[myid][i]>=graph->dist.offset,"Generated global "
              "coarse vertex number is smaller than graph offset (gvtx = "
              PF_VTX_T", offset = "PF_VTX_T"\n",v,graph->dist.offset);
          ++mycnvtxs;
        } else {
          if (maxidx < mynvtxs) {
            o = myid;
            l = maxidx;
          } else {
            o = gvtx_to_tid(maxidx,graph->dist);
            l = gvtx_to_lvtx(maxidx,graph->dist);
          }
          if (match[o][l] == NULL_VTX) {
            /* match with unmatched foreign vertex */
            if (iadjwgt) {
              cdeg[myid][mycnvtxs] = xdeg + iadjwgt[o][l] + eadjwgt[o][l];
            } else {
              cdeg[myid][mycnvtxs] = xdeg + eadjwgt[o][l];
            }
            csize[myid][mycnvtxs] = 2;
            //fmap[myid][mycnvtxs] = i;
            cg = lvtx_to_gvtx(mycnvtxs++,myid,graph->dist);
            DL_ASSERT(cg>=graph->dist.offset,"Generated global coarse vertex "
                "number is smaller than graph offset (gvtx = "PF_VTX_T
                ", offset = "PF_VTX_T"\n",cg,graph->dist.offset);
            cmap[myid][i] = cmap[o][l] = cg;
          } else {
            /* match with already matched foreign vertex */
            cg = cmap[myid][i] = cmap[o][l];
            cl = gvtx_to_lvtx(cg,graph->dist);
            co = gvtx_to_tid(cg,graph->dist);
            cdeg[co][cl] += xdeg;
            ++csize[co][cl];
          }
        }
        /* assign match vector */
        collapsed += __cluster(i,maxidx,myid,match,graph);
      }
    }
    vw_ht_free(conn);
  }

  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    vtx_set(cmap[myid],NULL_VTX,mynvtxs);
    dl_free(cdeg[myid]);
    dl_free(csize[myid]);
  }

  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    agg->nvtxs[myid] = __cluster_cleanup(mynvtxs,myid,match,cmap,fmap,graph);
  }

  dl_free(csize);
  dl_free(cdeg);
  dl_free(perm);

  /* aggregation stat tracking */
  if (objective->aggstats) {
    __calc_aggstats(objective,graph,agg);
  }

  return agg;
}




/******************************************************************************
* PRIVATE PARALLEL AGGREGATION FUNCTIONS **************************************
******************************************************************************/


static void __par_calc_aggstats(
    objective_t * const objective, 
    graph_t const * const graph, 
    aggregate_t const * const agg)
{
  vtx_t i,j,k,v,tcnvtxs;
  tid_t o;
  aggregation_stat_t aggstat;

  tid_t const myid = dlthread_get_id(objective->comm);
  vtx_t const mynvtxs = graph->mynvtxs[myid];

  init_aggregation_stat(&aggstat);

  j = 0;
  for (i=0;i<mynvtxs;++i) {
    if (is_supernode(i,myid,graph,objective->snratio)) {
      ++(aggstat.nsupernodes);
    }
    if (agg->match[myid][i] == i) {
      ++(aggstat.nlonenodes);
    } else {
      v = 0;
      o = myid;
      k = i;
      do {
        ++v;
        k = agg->match[o][k];
        if (k >= graph->mynvtxs[o]) {
          o = gvtx_to_tid(k,graph->dist);
          k = gvtx_to_lvtx(k,graph->dist);
        }          
      } while (o != myid || k != i);
      dl_storemax(aggstat.maxnodesize,v);
      j += v;
    }
  }
  aggstat.nsupernodes = vtx_dlthread_sumreduce(aggstat.nsupernodes, \
      objective->comm);
  aggstat.nlonenodes = vtx_dlthread_sumreduce(aggstat.nlonenodes, \
      objective->comm);
  j = adj_dlthread_sumreduce(j,objective->comm);
  tcnvtxs = vtx_dlthread_sumreduce(agg->nvtxs[myid],objective->comm);
  aggstat.maxnodesize = vtx_dlthread_maxreduce_value(aggstat.maxnodesize, \
      objective->comm);
  aggstat.avgnodesize = j / (real_t)(tcnvtxs-aggstat.nlonenodes);
  aggstat.coarsenrate = graph->nvtxs / (real_t)tcnvtxs;

  if (myid == 0) {
    aggregation_stats_push(aggstat,objective->aggstats);
  }
}


/**
 * @brief Create an aggregation using random matching.
 *
 * @param objective The objective containing aggregation parameters
 * @param mgraph The mgraph to aggregate 
 *
 * @return The new aggregation
 */
static aggregate_t * __par_aggregate_RM(
    objective_t * const objective, 
    mgraph_t const * const mgraph)
{
  vtx_t i,v,k,l,maxidx,collapsed;
  adj_t j;
  tid_t o;

  tid_t const myid = dlthread_get_id(objective->comm);

  /* make accessing my graph easy */
  graph_t const * const graph = mgraph->graph;
  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy;

  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const * const myxadj = xadj[myid];
  vtx_t const * const myadjncy = adjncy[myid];

  /* setup my aggregation information */
  aggregate_t * const agg = par_setup_aggregate(mynvtxs,NULL_VTX,NULL,NULL,
      NULL,objective->comm);

  vtx_t ** const match = agg->match;
  vtx_t ** const cmap = agg->cmap;
  vtx_t ** const fmap = agg->fmap;
  
  /* permutation array */
  vtx_t * perm = vtx_alloc(mynvtxs);

  if (myid == 0) {
    dprintf("aggregate_RM() called with a graph containing "PF_VTX_T
        " vertices\n",mgraph->graph->nvtxs);
  }

  vtx_incset(perm,0,1,mynvtxs);
  vtx_shuffle_r(perm,mynvtxs,&objective->seed);
  vtx_set(match[myid],NULL_VTX,mynvtxs);

  collapsed = 0; 

  dlthread_barrier(objective->comm);

  /* generate matches */
  for (v=0;v<mynvtxs;++v) {
    if (__agg_limit_reached(mynvtxs,collapsed,objective->max_agg_rate)) {
      break;
    }
    i = perm[v];
    if (match[myid][i] == NULL_VTX) {
      maxidx = i;
      for (j=myxadj[i];j<myxadj[i+1];++j) {
        k = myadjncy[j];
        if (k < mynvtxs) {
          if (match[myid][k] == NULL_VTX) {
            maxidx = k;
            break;
          }
        } else {
          o = gvtx_to_tid(k,graph->dist);
          l = gvtx_to_lvtx(k,graph->dist);
          if (match[o][l] == NULL_VTX) {
            maxidx = k;
            break;
          }
        }
      }
      collapsed += __match(i,maxidx,myid,match,graph);
    }
  }

  dl_free(perm);

  dlthread_barrier(objective->comm);

  agg->nvtxs[myid] = __match_cleanup(mynvtxs,myid,match,cmap,fmap,graph);

  /* aggregation stat tracking */
  if (objective->aggstats) {
    __par_calc_aggstats(objective,graph,agg);
  }

  return agg;
}


/**
 * @brief Create an aggregation using random clustering.
 *
 * @param objective The objective containing aggregation parameters
 * @param mgraph The mgraph to aggregate 
 *
 * @return The new aggregation
 */
static aggregate_t * __par_aggregate_RC(
    objective_t * const objective, 
    mgraph_t const * const mgraph)
{
  vtx_t i,v,k,l,maxidx,collapsed,mycnvtxs,cg;
  adj_t j;
  tid_t o;

  tid_t const myid = dlthread_get_id(objective->comm);

  /* make accessing my graph easy */
  graph_t const * const graph = mgraph->graph;
  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy;

  adj_t const * const myxadj = xadj[myid];
  vtx_t const * const myadjncy = adjncy[myid];

  /* setup my aggregation information */
  aggregate_t * const agg = par_setup_aggregate(mynvtxs,NULL_VTX,NULL,NULL, \
      NULL,objective->comm);

  vtx_t ** const match = agg->match;
  vtx_t ** const cmap = agg->cmap;
  vtx_t ** const fmap = agg->fmap;
  
  /* permutation array */
  vtx_t * perm = vtx_alloc(mynvtxs);

  if (myid == 0) {
    dprintf("aggregate_RC() called with a graph containing "PF_VTX_T \
        " vertices\n",mgraph->graph->nvtxs);
  }

  vtx_incset(perm,0,1,mynvtxs);
  vtx_shuffle_r(perm,mynvtxs,&objective->seed);
  vtx_set(match[myid],NULL_VTX,mynvtxs);

  mycnvtxs = 0;
  collapsed = 0;

  dlthread_barrier(objective->comm);

  /* match everyone possible on first pass */
  for (v=0;v<mynvtxs;++v) {
    if (__agg_limit_reached(mynvtxs,collapsed,objective->max_agg_rate)) {
      break;
    }
    i = perm[v];
    if (match[myid][i] == NULL_VTX) {
      maxidx = i;
      /* try and find an umatched vertex first */
      for (j=myxadj[i];j<myxadj[i+1];++j) {
        k = myadjncy[j];
        if (k < mynvtxs) {
          /* a vertex I own */
          if (match[myid][k] == NULL_VTX) {
            maxidx = k;
            break;
          }
        } else {
          /* vertex someone else owns */
          o = gvtx_to_tid(k,graph->dist);
          l = gvtx_to_lvtx(k,graph->dist);
          if (match[o][l] == NULL_VTX) {
            maxidx = k;
            break;
          }
        }
      }
      /* if I didn't find a vertex an unmatched vertex, just pick my first 
       * neighbor */
      if (j == myxadj[i+1] && j > myxadj[i]) {
        maxidx = myadjncy[myxadj[i]];
      }
      /* set match */
      if (maxidx == i) {
        cmap[myid][i] = lvtx_to_gvtx(mycnvtxs++,myid,graph->dist);
        DL_ASSERT(cmap[myid][i]>=graph->dist.offset,"Generated global "
            "coarse vertex number is smaller than graph offset (gvtx = "
            PF_VTX_T", offset = "PF_VTX_T"\n",v,graph->dist.offset);
      } else {
        if (maxidx < mynvtxs) {
          o = myid;
          l = maxidx;
        } else {
          o = gvtx_to_tid(maxidx,graph->dist);
          l = gvtx_to_lvtx(maxidx,graph->dist);
        }
        if (match[o][l] == NULL_VTX) { /* duh */
          cg = lvtx_to_gvtx(mycnvtxs++,myid,graph->dist);
          DL_ASSERT(cg>=graph->dist.offset,"Generated global coarse vertex "
              "number is smaller than graph offset (gvtx = "PF_VTX_T
              ", offset = "PF_VTX_T"\n",cg,graph->dist.offset);
          cmap[myid][i] = cmap[o][l] = cg;
        } else {
          cg = cmap[myid][i] = cmap[o][l];
        }
      }
      collapsed += __cluster(i,maxidx,myid,match,graph);
    }
  }
  dl_free(perm);

  dlthread_barrier(objective->comm);

  vtx_set(cmap[myid],NULL_VTX,mynvtxs);

  dlthread_barrier(objective->comm);

  agg->nvtxs[myid] = __cluster_cleanup(mynvtxs,myid,match,cmap,fmap,graph);

  dlthread_barrier(objective->comm);

  /* aggregation stat tracking */
  if (objective->aggstats) {
    __par_calc_aggstats(objective,graph,agg);
  }

  return agg;
}


/**
 * @brief Create an aggregation using agglomerative modularity based matching.
 *
 * @param objective The objective containing aggregation parameters
 * @param mgraph The mgraph to aggregate 
 *
 * @return The new aggregation
 */
static aggregate_t * __par_aggregate_AGM(
    objective_t * const objective, 
    mgraph_t const * const mgraph)
{
  vtx_t v,i,k,l,maxidx,collapsed;
  adj_t j;
  wgt_t xdeg,ydeg,cwgt,xdegm;
  tid_t o;
  real_t cmod, maxmod;

  tid_t const myid = dlthread_get_id(objective->comm);

  /* set up graph stuff */
  graph_t const * graph = mgraph->graph;
  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy;
  wgt_t const * const * const adjwgt = (wgt_t const * const *)graph->adjwgt;
  wgt_t const * const * const iadjwgt = (wgt_t const * const *)graph->iadjwgt;
  wgt_t const * const * const eadjwgt = (wgt_t const * const *)graph->eadjwgt;

  adj_t const * myxadj = xadj[myid];
  vtx_t const * myadjncy = adjncy[myid];

  wgt_t const mscale = objective->degwgt/graph->gadjwgt;

  DL_ASSERT_EQUALS((tid_t)dlthread_get_nthreads(objective->comm),graph->npar, \
      PF_TID_T);

  /* setup my aggregation information */
  aggregate_t * const agg = par_setup_aggregate(mynvtxs,NULL_VTX,NULL,NULL, \
      NULL,objective->comm);

  vtx_t ** const match = agg->match;
  vtx_t ** const cmap = agg->cmap;
  vtx_t ** const fmap = agg->fmap;

  vtx_t * perm = vtx_alloc(mynvtxs);

  /* do the real work */
  vtx_set(match[myid],NULL_VTX,mynvtxs);

  collapsed = 0;

  vtx_incset(perm,0,1,mynvtxs);
  vtx_shuffle_r(perm,mynvtxs,&objective->seed);

  dlthread_barrier(objective->comm);

  for (v=0;v<mynvtxs;++v) {
    if (__agg_limit_reached(mynvtxs,collapsed,objective->max_agg_rate)) {
      break;
    }
    i = perm[v];
    if (match[myid][i] == NULL_VTX) {
      if (iadjwgt) {
        xdeg = eadjwgt[myid][i] + iadjwgt[myid][i];
      } else {
        xdeg = eadjwgt[myid][i];
      }
      xdegm = xdeg*mscale;
      if (objective->parttype == NERSTRAND_PARTITION_ANYWAY) {
        maxmod =  0;
      } else {
        maxmod = -graph->gadjwgt;
      }
      maxidx = i;
      for (j=myxadj[i];j<myxadj[i+1];++j) {
        k = myadjncy[j];
        if (adjwgt) {
          cwgt = adjwgt[myid][j];
        } else {
          cwgt = 1.0;
        }
        if (k < mynvtxs) {
          o = myid;
          l = k;
        } else {
          o = gvtx_to_tid(k,graph->dist);
          l = gvtx_to_lvtx(k,graph->dist);
        }
        if (iadjwgt) {
          ydeg = iadjwgt[o][l] + eadjwgt[o][l];
        } else {
          ydeg = eadjwgt[o][l];
        }
        cmod = cwgt - (xdegm*ydeg);
        if (match[o][l] == NULL_VTX && maxmod < cmod) {
          maxidx = k;
          maxmod = cmod;
        }
      }
      collapsed += __match(i,maxidx,myid,match,graph);
    }
  }

  dl_free(perm);

  dlthread_barrier(objective->comm);

  agg->nvtxs[myid] = __match_cleanup(mynvtxs,myid,match,cmap,fmap,graph);

  dlthread_barrier(objective->comm);

  /* aggregation stat tracking */
  if (objective->aggstats) {
    __par_calc_aggstats(objective,graph,agg);
  }

  return agg;
}


/**
 * @brief Create an aggregation using agglomerative modularity based matching,
 * with 2-hop matching.
 *
 * @param objective The objective containing aggregation parameters
 * @param mgraph The mgraph to aggregate 
 *
 * @return The new aggregation
 */
static aggregate_t * __par_aggregate_AGH(
    objective_t * const objective, 
    mgraph_t const * const mgraph)
{
  vtx_t v,i,k,kk,l,ll,maxidx,collapsed,mysize,n;
  adj_t j,jj;
  wgt_t xdeg,ydeg,cwgt,xdegm;
  tid_t o,oo;
  real_t cmod, maxmod;

  tid_t const myid = dlthread_get_id(objective->comm);

  /* set up graph stuff */
  graph_t const * graph = mgraph->graph;
  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy;
  wgt_t const * const * const adjwgt = (wgt_t const * const *)graph->adjwgt;
  wgt_t const * const * const iadjwgt = (wgt_t const * const *)graph->iadjwgt;
  wgt_t const * const * const eadjwgt = (wgt_t const * const *)graph->eadjwgt;

  adj_t const * myxadj = xadj[myid];
  vtx_t const * myadjncy = adjncy[myid];

  wgt_t const mscale = objective->degwgt/graph->gadjwgt;

  DL_ASSERT_EQUALS((tid_t)dlthread_get_nthreads(objective->comm),graph->npar, \
      PF_TID_T);

  /* setup my aggregation information */
  aggregate_t * const agg = par_setup_aggregate(mynvtxs,NULL_VTX,NULL,NULL, \
      NULL,objective->comm);

  vtx_t ** const match = agg->match;
  vtx_t ** const cmap = agg->cmap;
  vtx_t ** const fmap = agg->fmap;

  vtx_t * perm = vtx_alloc(mynvtxs);

  /* do the real work */
  vtx_set(match[myid],NULL_VTX,mynvtxs);

  collapsed = 0;

  /* counting sort stuff */
  vtx_incset(perm,0,1,mynvtxs);
  vtx_shuffle_r(perm,mynvtxs,&objective->seed);

  dlthread_barrier(objective->comm);

  for (v=0;v<mynvtxs;++v) {
    if (__agg_limit_reached(mynvtxs,collapsed,objective->max_agg_rate)) {
      break;
    }
    i = perm[v];
    if (match[myid][i] == NULL_VTX) {
      if (iadjwgt) {
        xdeg = eadjwgt[myid][i] + iadjwgt[myid][i];
      } else {
        xdeg = eadjwgt[myid][i];
      }
      xdegm = xdeg*mscale;
      if (objective->parttype == NERSTRAND_PARTITION_ANYWAY) {
        maxmod =  0;
      } else {
        maxmod = -graph->gadjwgt;
      }
      maxidx = i;
      mysize = myxadj[i+1] - myxadj[i];
      for (j=myxadj[i];j<myxadj[i+1];++j) {
        k = myadjncy[j];
        if (adjwgt) {
          cwgt = adjwgt[myid][j];
        } else {
          cwgt = 1.0;
        }
        if (k < mynvtxs) {
          o = myid;
          l = k;
        } else {
          o = gvtx_to_tid(k,graph->dist);
          l = gvtx_to_lvtx(k,graph->dist);
        }
        if (iadjwgt) {
          ydeg = iadjwgt[o][l] + eadjwgt[o][l];
        } else {
          ydeg = eadjwgt[o][l];
        }
        cmod = cwgt - (xdegm*ydeg);
        if (match[o][l] == NULL_VTX && maxmod < cmod) {
          maxidx = k;
          maxmod = cmod;
        }
      }
      /* two hop */
      if (maxidx == i && mysize <= MAXDEG2HOP) {
        for (j=myxadj[i];j<myxadj[i+1];++j) {
          k = myadjncy[j];
          if (k < mynvtxs) {
            o = myid;
            l = k;
          } else {
            o = gvtx_to_tid(k,graph->dist);
            l = gvtx_to_lvtx(k,graph->dist);
          }
          n = 0;
          jj = xadj[o][l] + (i % (xadj[o][l+1]-xadj[o][l]));
          for (n=0;n<MAXSEARCH2HOP&&n<xadj[o][l+1]-xadj[o][l];++n,++jj) {
            if (jj == xadj[o][l+1]) {
              jj = xadj[o][l];
            }
            kk = adjncy[o][jj];
            if (kk < graph->mynvtxs[o]) {
              ll = kk;
              kk = lvtx_to_gvtx(ll,o,graph->dist);
              oo = o;
            } else {
              oo = gvtx_to_tid(kk,graph->dist);
              ll = gvtx_to_lvtx(kk,graph->dist);
            }
            if (match[oo][ll] == NULL_VTX && !(ll == i && oo == myid) && 
                xadj[oo][ll+1] - xadj[oo][ll] <= MAXDEG2HOP) {
              if (oo == myid) {
                maxidx = gvtx_to_lvtx(kk,graph->dist);
              } else {
                maxidx = kk;
              }
              goto ENDTWOHOP;
            }
          }
        }
        ENDTWOHOP:;
      }
      collapsed += __match(i,maxidx,myid,match,graph);
    }
  }

  dl_free(perm);

  dlthread_barrier(objective->comm);

  agg->nvtxs[myid] = __match_cleanup(mynvtxs,myid,match,cmap,fmap,graph);

  dlthread_barrier(objective->comm);

  /* aggregation stat tracking */
  if (objective->aggstats) {
    __par_calc_aggstats(objective,graph,agg);
  }

  return agg;
}


wgt_t ** __pac_cdeg = NULL;
vtx_t ** __pac_csize = NULL;
/**
 * @brief Create an aggregation using agglomerative modularity based
 * clustering.
 *
 * @param objective The objective containing aggregation parameters
 * @param mgraph The mgraph to aggregate 
 *
 * @return The new aggregation
 */
static aggregate_t * __par_aggregate_AGC(
    objective_t * const objective, 
    mgraph_t const * const mgraph)
{
  vtx_t v,i,k,l,cl,cg,maxidx,collapsed,mycnvtxs;
  tid_t o,co;
  adj_t j, astart,aend;
  mod_t cmod, maxmod;
  wgt_t xdeg, xdegm, ydeg, cwgt, pen;

  vw_ht_t * conn;

  tid_t const myid = dlthread_get_id(objective->comm);
  tid_t const nthreads = dlthread_get_nthreads(objective->comm);

  /* set up graph stuff */
  graph_t const * graph = mgraph->graph;
  vtx_t const nvtxs = graph->nvtxs;
  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy;
  wgt_t const * const * const adjwgt = (wgt_t const * const *)graph->adjwgt;
  wgt_t const * const * const iadjwgt = (wgt_t const * const *)graph->iadjwgt;
  wgt_t const * const * const eadjwgt = (wgt_t const * const *)graph->eadjwgt;

  adj_t const * myxadj = xadj[myid];
  vtx_t const * myadjncy = adjncy[myid];

  wgt_t const mscale = objective->degwgt/graph->gadjwgt;

  DL_ASSERT_EQUALS(nthreads,graph->npar,PF_TID_T);

  /* setup my aggregation information */
  aggregate_t * const agg = par_setup_aggregate(mynvtxs,NULL_VTX,NULL,NULL, \
      NULL,objective->comm);

  vtx_t ** const match = agg->match;
  vtx_t ** const cmap = agg->cmap;
  vtx_t ** const fmap = agg->fmap;

  vtx_t * perm = vtx_alloc(nvtxs);

  /* do the real work */
  vtx_set(match[myid],NULL_VTX,mynvtxs);

  /* store total degree of coarse vertices */
  if (myid == 0) {
    __pac_cdeg = r_wgt_alloc(nthreads);
    __pac_csize = r_vtx_alloc(nthreads);
  }
  dlthread_barrier(objective->comm);

  __pac_cdeg[myid] = wgt_alloc(mynvtxs);
  __pac_csize[myid] = vtx_alloc(mynvtxs);

  dlthread_barrier(objective->comm);

  mycnvtxs = 0;
  conn = vw_ht_create(graph->maxdeg[myid],graph->maxdeg[myid]/4);
  collapsed = 0;

  vtx_incset(perm,0,1,mynvtxs);
  vtx_shuffle_r(perm,mynvtxs,&objective->seed);

  for (v=0;v<mynvtxs;++v) {
    if (__agg_limit_reached(mynvtxs,collapsed,objective->max_agg_rate)) {
      break;
    }
    i = perm[v];
    if (match[myid][i] == NULL_VTX) {
      if (iadjwgt) {
        xdeg = iadjwgt[myid][i] + eadjwgt[myid][i]; 
      } else {
        xdeg = eadjwgt[myid][i];
      }
      xdegm = xdeg*mscale;
      if (objective->parttype == NERSTRAND_PARTITION_ANYWAY) {
        maxmod =  0;
      } else {
        maxmod = -graph->gadjwgt;
      }
      maxidx = i;
      astart = myxadj[i];
      aend = myxadj[i+1];
      if (astart == aend) {
        if (objective->parttype == NERSTRAND_PARTITION_KWAY && \
            v < mynvtxs-1 && \
            objective->max_agg_rate*graph->nvtxs > 2*objective->cnvtxs) {
          maxidx = perm[v+1];
        } 
      } else if (astart == aend+1) {
        /* if we only have one option, always match it. Brandes et al. '08
         * "On Modularity", showed that 1-degree vertices should never be in
         * their own cluster. */
        maxidx = myadjncy[astart];
      } else {
        /* hopefully for many graphs this will fit in cache */
        for (j=astart;j<aend;++j) {
          k = myadjncy[j];
          if (k < mynvtxs) {
            o = myid;
            l = k;
          } else {
            o = gvtx_to_tid(k,graph->dist);
            l = gvtx_to_lvtx(k,graph->dist);
          }
          if (match[o][l] != NULL_VTX) {
            if (adjwgt) {
              vw_ht_add(cmap[o][l],adjwgt[myid][j],conn);
            } else {
              vw_ht_add(cmap[o][l],1.0,conn);
            }
          }
        }
        for (j=astart;j<aend;++j) {
          k = myadjncy[j];
          if (k < mynvtxs) {
            o = myid;
            l = k;
          } else {
            o = gvtx_to_tid(k,graph->dist);
            l = gvtx_to_lvtx(k,graph->dist);
          }
          if (match[o][l] == NULL_VTX) {
            if (iadjwgt) {
              ydeg = iadjwgt[o][l] + eadjwgt[o][l];
            } else {
              ydeg = eadjwgt[o][l];
            }
            if (adjwgt) {
              cwgt = adjwgt[myid][j];
            } else {
              cwgt = 1.0;
            }
            pen = 1;
          } else {
            cl = gvtx_to_lvtx(cmap[o][l],graph->dist);
            co = gvtx_to_tid(cmap[o][l],graph->dist);
            ydeg = __pac_cdeg[co][cl];
            if (co != myid && __pac_csize[co][cl] >= MAX_CYCLE_SIZE) {
              continue;
            } else if (__pac_csize[co][cl] > MAX_VERTEX_SIZE) {
              pen = DESIRED_VERTEX_SIZE / (wgt_t)__pac_csize[co][cl];
            } else if (__pac_csize[co][cl] > DESIRED_VERTEX_SIZE) {
              pen = ((0.5*DESIRED_VERTEX_SIZE) / __pac_csize[co][cl]) + 0.5;
            } else { 
              pen = 1;
            }
            cwgt = vw_ht_get(cmap[o][l],conn);
          }
          cmod = (cwgt - (ydeg*xdegm))*pen; 
          if (maxidx > mynvtxs) {
            cmod *= REMOTE_PENALTY;
          }
          if (maxmod <= cmod) {
            maxidx = k;
            maxmod = cmod;
          }
        }
        /* clear the conn map */
        vw_ht_clear_chains(conn);
        for (j=astart;j<aend;++j) {
          k = myadjncy[j];
          if (k < mynvtxs) {
            o = myid;
            l = k;
          } else {
            o = gvtx_to_tid(k,graph->dist);
            l = gvtx_to_lvtx(k,graph->dist);
          }
          if (match[o][l] != NULL_VTX) {
            vw_ht_clear_slot(cmap[o][l],conn);
          }
        }
      }
      /* match maker make me a match */
      if (maxidx == i) {
        __pac_cdeg[myid][mycnvtxs] = xdeg;
        __pac_csize[myid][mycnvtxs] = 1;
        cmap[myid][i] = lvtx_to_gvtx(mycnvtxs++,myid,graph->dist);
        DL_ASSERT(cmap[myid][i]>=graph->dist.offset,"Generated global "
            "coarse vertex number is smaller than graph offset (gvtx = "
            PF_VTX_T", offset = "PF_VTX_T"\n",v,graph->dist.offset);
      } else {
        if (maxidx < mynvtxs) {
          o = myid;
          l = maxidx;
        } else {
          o = gvtx_to_tid(maxidx,graph->dist);
          l = gvtx_to_lvtx(maxidx,graph->dist);
        }
        if (match[o][l] == NULL_VTX) { /* duh */
          if (iadjwgt) {
            __pac_cdeg[myid][mycnvtxs] = xdeg + iadjwgt[o][l] + eadjwgt[o][l];
          } else {
            __pac_cdeg[myid][mycnvtxs] = xdeg + eadjwgt[o][l];
          }
          __pac_csize[myid][mycnvtxs] = 2;
          cg = lvtx_to_gvtx(mycnvtxs++,myid,graph->dist);
          DL_ASSERT(cg>=graph->dist.offset,"Generated global coarse vertex "
              "number is smaller than graph offset (gvtx = "PF_VTX_T
              ", offset = "PF_VTX_T"\n",cg,graph->dist.offset);
          cmap[myid][i] = cmap[o][l] = cg;
        } else {
          cg = cmap[myid][i] = cmap[o][l];
          cl = gvtx_to_lvtx(cg,graph->dist);
          co = gvtx_to_tid(cg,graph->dist);
          __pac_cdeg[co][cl] += xdeg;
          ++__pac_csize[co][cl];
        }
      }
      collapsed += __cluster(i,maxidx,myid,match,graph);
    }
  }
  vw_ht_free(conn);

  dlthread_barrier(objective->comm);

  /* reset the cmap so I perform cycle clean up */
  vtx_set(cmap[myid],NULL_VTX,mynvtxs);

  dl_free(__pac_cdeg[myid]);
  dl_free(__pac_csize[myid]);
  dl_free(perm);

  dlthread_barrier(objective->comm);

  /* need to make sure that cdeg[myid] is free'd */
  if (myid == 0) {
    dl_free(__pac_cdeg);
    dl_free(__pac_csize);
  }

  agg->nvtxs[myid] = __cluster_cleanup(mynvtxs,myid,match,cmap,fmap,graph);

  dlthread_barrier(objective->comm);

  /* aggregation stat tracking */
  if (objective->aggstats) {
    __par_calc_aggstats(objective,graph,agg);
  }


  return agg;
}




/******************************************************************************
* PUBLIC SERIAL FUNCTIONS *****************************************************
******************************************************************************/


aggregate_t * init_aggregate(
    aggregate_t * const agg) 
{
  agg->cmap = NULL;
  agg->match = NULL;
  agg->free_cmap = 0;
  agg->free_match = 0;
  agg->free_fmap = 0;
  agg->free_nvtxs = 0;
  return agg;
}


aggregate_t * setup_aggregate(
    vtx_t const * const nvtxs, 
    vtx_t * const cnvtxs, 
    vtx_t ** match, 
    vtx_t ** cmap,
    vtx_t ** fmap, 
    tid_t const nthreads)
{
  aggregate_t * const agg = aggregate_calloc(1);
  agg->npar = nthreads;
  if (match) {
    agg->match = match;
    agg->free_match = 0;
  } else {
    DL_ASSERT(nvtxs != NULL,"Passed in NULL nvtxs and match\n");
    agg->match = r_vtx_dalloc(nvtxs,sizeof(vtx_t),nthreads);
    agg->free_match = 1;
  }
  if (cmap) {
    agg->cmap = cmap;
    agg->free_cmap = 0;
  } else {
    DL_ASSERT(nvtxs != NULL,"Passed in NULL nvtxs and cmap\n");
    agg->cmap = r_vtx_dalloc(nvtxs,sizeof(vtx_t),nthreads);
    agg->free_cmap = 1;
  }
  if (fmap) {
    agg->fmap = fmap;
    agg->free_fmap = 0;
  } else {
    DL_ASSERT(nvtxs != NULL || cnvtxs != NULL,
        "Passed in NULL cnvtxs, nvtxs and fmap\n");
    if (cnvtxs) {
      agg->fmap = r_vtx_dalloc(cnvtxs,sizeof(vtx_t),nthreads);
    } else {
      agg->fmap = r_vtx_dalloc(nvtxs,sizeof(vtx_t),nthreads);
    }
    agg->free_fmap = 1;
  }
  if (cnvtxs) {
    agg->nvtxs = cnvtxs;
  } else {
    agg->nvtxs = vtx_init_alloc(NULL_VTX,nthreads);
    agg->free_nvtxs = 1;
  }
  return agg;  
}


int free_aggregate(
    aggregate_t * agg)
{
  tid_t const nthreads = agg->npar;
  if (agg->free_cmap) {
    r_vtx_free(agg->cmap,nthreads);
  }
  if (agg->free_match) {
    r_vtx_free(agg->match,nthreads);
  }
  if (agg->free_fmap) {
    r_vtx_free(agg->fmap,nthreads);
  }
  if (agg->free_nvtxs) {
    dl_free(agg->nvtxs);
  }
  dl_free(agg);

  return 1;
}


aggregate_t * aggregate_graph(
    objective_t * const objective,
    mgraph_t const * const mgraph)
{
  dl_start_timer(&(objective->timers.aggregation));

  DL_ASSERT_EQUALS(check_graph(mgraph->graph),1,"%d");

  aggregate_t * agg;  

  switch (objective->aggtype) {
    case NERSTRAND_AGGREGATE_RM:
      dprintf("Using AGGREGATE_RM to aggregate graph of "PF_VTX_T" vertices\n",
          mgraph->graph->nvtxs);
      agg = __aggregate_RM(objective,mgraph);
      break;
    case NERSTRAND_AGGREGATE_AGM: 
      dprintf("Using AGGREGATE_AGM to aggregate graph of "PF_VTX_T
          " vertices\n", mgraph->graph->nvtxs);
      agg = __aggregate_AGM(objective,mgraph);
      break;
    case NERSTRAND_AGGREGATE_AGH:
      dprintf("Using AGGREGATE_AGH to aggregate graph of "PF_VTX_T
          " vertices\n",mgraph->graph->nvtxs);
      agg = __aggregate_AGH(objective,mgraph);
      break;
    case NERSTRAND_AGGREGATE_RC:
      dprintf("Using AGGREGATE_RC to aggregate graph of "PF_VTX_T
          " vertices\n", mgraph->graph->nvtxs);
      agg = __aggregate_RC(objective,mgraph);
      break;
    case NERSTRAND_AGGREGATE_AGC:
      dprintf("Using AGGREGATE_AGC to aggregate graph of "PF_VTX_T
          " vertices\n", mgraph->graph->nvtxs);
      agg = __aggregate_AGC(objective,mgraph);
      break;
    default:
      dl_error("Unknown aggregation type %d\n",objective->aggtype);
  }

  /* make sure we know what to free */
  agg->free_fmap = 1;
  agg->free_match = 1;
  agg->free_cmap = 0;

  DL_ASSERT_EQUALS(check_aggregate(agg,mgraph->graph),1,"%d");

  dl_stop_timer(&(objective->timers.aggregation));
  
  return agg;
}


int check_aggregate(
    aggregate_t const * const aggregate, 
    graph_t const * const graph) 
{
  vtx_t i,c,n,tcnvtxs;
  tid_t myid, o;

  DL_ASSERT_EQUALS(aggregate->npar,graph->npar,PF_TID_T);

  tid_t const nthreads = aggregate->npar;
  vtx_t const * const nvtxs = graph->mynvtxs;
  vtx_t const * const cnvtxs = aggregate->nvtxs;
  vtx_t const * const * const match = (vtx_t const * const *)aggregate->match;
  vtx_t const * const * const cmap = (vtx_t const * const *)aggregate->cmap;
  vtx_t const * const * const fmap = (vtx_t const * const *)aggregate->fmap;

  vtx_t ** mk = r_vtx_alloc(nthreads); 

  n = 0;
  for (myid=0;myid<nthreads;++myid) {
    mk[myid] = vtx_init_alloc(NULL_VTX,nvtxs[myid]);
  }

  tcnvtxs = vtx_sum(cnvtxs,nthreads);

  for (myid=0;myid<nthreads;++myid) {
    for (c=0;c<cnvtxs[myid];++c) {
      /* check c's chain */
      i = fmap[myid][c];

      DL_ASSERT(i < graph->mynvtxs[myid],"First vertex in coarse vertex "
          PF_VTX_T" is "PF_VTX_T"/"PF_VTX_T"\n",c,i,graph->mynvtxs[myid]);
      DL_ASSERT(cmap[myid][i] >= graph->dist.offset,"Local vertex number "
          "stored in cmap (i = "PF_VTX_T", cmap[i] = "PF_VTX_T", offset = "
          PF_VTX_T")\n",i,cmap[myid][i],graph->dist.offset);
      DL_ASSERT_EQUALS(c,gvtx_to_lvtx(cmap[myid][i],graph->dist),PF_VTX_T);

      o = myid;
      do {
        if (gvtx_to_lvtx(cmap[o][i],graph->dist) != c) {
          eprintf("Found vertex "PF_VTX_T" belonging in to coarse vertex "
              PF_VTX_T" in the chain of coarse vertex "PF_VTX_T" (g:"PF_VTX_T
              ")\n",i,cmap[o][i],c,lvtx_to_gvtx(c,myid,graph->dist));
          return 0;
        }
        mk[o][i] = i;
        ++n;
        i = match[o][i];
        if (i >= nvtxs[o]) {
          o = gvtx_to_tid(i,graph->dist);
          i = gvtx_to_lvtx(i,graph->dist);
        }
      } while (!(i == fmap[myid][c] && o == myid));
    }
  }
  for (myid=0;myid<nthreads;++myid) {
    for (i=0;i<nvtxs[myid];++i) {
      if (mk[myid][i] == NULL_VTX) {
        eprintf("Found unvisited vertex "PF_VTX_T", with cmap = "PF_VTX_T" and "
            "match = "PF_VTX_T", owned by "PF_TID_T"\n",i,cmap[myid][i],
            match[myid][i],myid);
        eprintf("Thread "PF_TID_T" owns "PF_VTX_T"/"PF_VTX_T" coarse "
            "vertices\n",myid,cnvtxs[myid],tcnvtxs);
        return 0;
      }
    }
    dl_free(mk[myid]);
  }
  dl_free(mk);

  if (n < graph->nvtxs) {
    eprintf("Visited "PF_VTX_T"/"PF_VTX_T" vertices in fmap+match arrays\n",
        n,graph->nvtxs);
    return 0;
  }

  return 1;
}




/******************************************************************************
* PUBLIC PARALLEL FUNCTIONS ***************************************************
******************************************************************************/


aggregate_t * par_setup_aggregate(
    vtx_t const mynvtxs, 
    vtx_t const mycnvtxs, 
    vtx_t * const match, 
    vtx_t * const cmap, 
    vtx_t * const fmap,
    dlthread_comm_t const comm)
{
  static aggregate_t * agg;

  tid_t const myid = dlthread_get_id(comm);
  tid_t const nthreads = dlthread_get_nthreads(comm);

  if (myid == 0) {
    agg = aggregate_calloc(1);
    agg->npar = nthreads;
    agg->match = r_vtx_alloc(nthreads);
    agg->free_match = 1;
    agg->cmap = r_vtx_alloc(nthreads);
    agg->free_cmap = 1;
    agg->fmap = r_vtx_alloc(nthreads);
    agg->free_fmap = 1;
    agg->nvtxs = vtx_alloc(nthreads);
    agg->free_nvtxs = 1;
    agg->comm = comm;
  }

  dlthread_barrier(comm);

  agg->nvtxs[myid] = mycnvtxs; 
  if (match) {
    agg->match[myid] = match;
  } else {
    agg->match[myid] = vtx_alloc(mynvtxs);
  }
  if (cmap) {
    agg->cmap[myid] = cmap;
  } else {
    agg->cmap[myid] = vtx_alloc(mynvtxs);
  }
  if (fmap) {
    agg->fmap[myid] = fmap;
  } else {
    if (mycnvtxs != NULL_VTX) {
      agg->fmap[myid] = vtx_alloc(mycnvtxs);
    } else {
      agg->fmap[myid] = vtx_alloc(mynvtxs);
    }
  }

  dlthread_barrier(comm);

  return agg;  
}


int par_free_aggregate(
    aggregate_t * agg)
{
  dlthread_comm_t comm;

  tid_t const myid = dlthread_get_id(agg->comm);

  comm = agg->comm;

  if (agg->free_cmap) {
    dl_free(agg->cmap[myid]);
  }
  if (agg->free_match) {
    dl_free(agg->match[myid]);
  }
  if (agg->free_fmap) {
    dl_free(agg->fmap[myid]);
  }

  dlthread_barrier(comm);
  if (myid == 0) {
    if (agg->free_cmap) {
      dl_free(agg->cmap);
    }
    if (agg->free_match) {
      dl_free(agg->match);
    }
    if (agg->free_fmap) {
      dl_free(agg->fmap);
    }
    if (agg->free_nvtxs) {
      dl_free(agg->nvtxs);
    }
    dl_free(agg);
  }
  dlthread_barrier(comm);

  return 1;
}


aggregate_t * par_aggregate_graph(
    objective_t * const objective,
    mgraph_t const * const mgraph)
{
  aggregate_t * agg;  

  tid_t const myid = dlthread_get_id(objective->comm);

  if (myid == 0) {
    dl_start_timer(&(objective->timers.aggregation));
  }

  DL_ASSERT_EQUALS(check_graph(mgraph->graph),1,"%d");

  switch (objective->aggtype) {
    case NERSTRAND_AGGREGATE_RM:
      if (myid == 0) {
        dprintf("Using AGGREGATE_RM to aggregate graph of "PF_VTX_T
            " vertices\n", mgraph->graph->nvtxs);
      }
      agg = __par_aggregate_RM(objective,mgraph);
      break;
    case NERSTRAND_AGGREGATE_AGM: 
      if (myid == 0) {
        dprintf("Using AGGREGATE_AGM to aggregate graph of "PF_VTX_T
            " vertices\n", mgraph->graph->nvtxs);
      }
      agg = __par_aggregate_AGM(objective,mgraph);
      break;
    case NERSTRAND_AGGREGATE_AGH:
      if (myid == 0) {
        dprintf("Using AGGREGATE_AGH to aggregate graph of "PF_VTX_T
            " vertices\n",mgraph->graph->nvtxs);
      }
      agg = __par_aggregate_AGH(objective,mgraph);
      break;
    case NERSTRAND_AGGREGATE_RC:
      if (myid == 0) {
        dprintf("Using AGGREGATE_RC to aggregate graph of "PF_VTX_T
            " vertices\n", mgraph->graph->nvtxs);
      }
      agg = __par_aggregate_RC(objective,mgraph);
      break;
    case NERSTRAND_AGGREGATE_AGC:
      if (myid == 0) {
        dprintf("Using AGGREGATE_AGC to aggregate graph of "PF_VTX_T
            " vertices\n", mgraph->graph->nvtxs);
      }
      agg = __par_aggregate_AGC(objective,mgraph);
      break;
    default:
      dl_error("Unknown aggregation type %d\n",objective->aggtype);
  }

  /* make sure we know what to free */
  agg->free_fmap = 1;
  agg->free_match = 1;
  agg->free_cmap = 0;

  dlthread_barrier(objective->comm);

  DL_ASSERT_EQUALS(check_aggregate(agg,mgraph->graph),1,"%d");
  
  if (myid == 0) {
    dl_stop_timer(&(objective->timers.aggregation));
  }

  return agg;
}


#endif
