/**
 * @file newman.c
 * @brief Functions for performing E. J. Newman's Fast Algorithm
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-10-24
 */




#ifndef NERSTRAND_NEWMAN_C
#define NERSTRAND_NEWMAN_C




#include "base.h"
#include "newman.h"




/******************************************************************************
* PRIVATE TYPES ***************************************************************
******************************************************************************/


typedef struct edge_t {
  adj_t radj;
  vtx_t u,v;
  wgt_t wgt;
} edge_t;


typedef struct edge_kv_t {
  mod_t key;
  vtx_t vtx;
  vtx_t offset; 
} edge_kv_t;


typedef struct edge_pq_t {
  edge_kv_t * elements;
  size_t size, maxsize;
  vtx_t ** index;
} edge_pq_t;




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLMEM_PREFIX edge 
#define DLMEM_TYPE_T edge_t
#define DLMEM_STATIC
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX


#define DLMEM_PREFIX edge_kv 
#define DLMEM_TYPE_T edge_kv_t
#define DLMEM_STATIC
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX


#define DLDJSET_PREFIX vtx
#define DLDJSET_TYPE_T vtx_t
#define DLDJSET_STATIC
#include "dldjset_headers.h"
#undef DLDJSET_STATIC
#undef DLDJSET_TYPE_T
#undef DLDJSET_PREFIX



/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


#define __LEFTCHILDINDEX(i) ((2*i)+1)
#define __RIGHTCHILDINDEX(i) ((2*i)+2)
#define __PARENTINDEX(i) ((i-1)>>1)


static int __check_edges(
    const edge_t * const * const edges, 
    const vtx_t * const sadj, 
    const vtx_t nvtxs)
{
  vtx_t i,j,k;
  for (i=0;i<nvtxs;++i) {
    if (sadj[i] > 0) {
      for (j=0;j<sadj[i];++j) {
        k = edges[i][j].u;
        if (edges[i][j].v != i) {
          eprintf("Bad 'v' in edge ["PF_VTX_T":"PF_VTX_T"]\n",i,j);
          return 0;
        }
        if (edges[i][j].radj >= sadj[k]) {
          eprintf("Invalid radj for ["PF_VTX_T":"PF_VTX_T"], radj: "PF_ADJ_T
              " >= sadj["PF_VTX_T"]: "PF_VTX_T"\n",i,j,edges[i][j].radj,k,
              sadj[k]);
          return 0;
        }
        if (edges[k][edges[i][j].radj].u != i) {
          eprintf("Bad radj in edge ["PF_VTX_T":"PF_VTX_T"] to "PF_VTX_T", "
              "radj = "PF_ADJ_T"\n",i,j,k,edges[i][j].radj);
          for (j=0;j<sadj[k];++j) {
            if (edges[k][j].u == i) {
              eprintf("Should be "PF_VTX_T"\n",j);
              return 0;
            }
          }
          eprintf("Couldn't find reverse edge!\n");
          return 0;
        }
        if (edges[k][edges[i][j].radj].v != k) {
          eprintf("Bad 'something' in edge ["PF_VTX_T":"PF_VTX_T"], v = "
              PF_VTX_T" instead of "PF_VTX_T"\n",i,j,
              edges[k][edges[i][j].radj].v,k);
          return 0;
        }
      }
    }
  }
  return 1;
}


static int __check_pq(
    const edge_pq_t * const q)
{
  vtx_t i,j,k;
  edge_kv_t e;

  for (i=0;i<q->size;++i) {
    e = q->elements[i];
    j = __LEFTCHILDINDEX(i);
    k = __RIGHTCHILDINDEX(i);
    if (q->index[e.vtx] == NULL) {
      eprintf("Edge with empty vertex list is in the queue: "PF_VTX_T":"
          PF_VTX_T"\n",e.vtx,e.offset);
      return 0;
    }
    if (q->index[e.vtx][e.offset] == NULL_VTX) {
      eprintf("Null edge in the priority queue at ["PF_VTX_T":"PF_VTX_T"]\n",
          e.vtx,e.offset);
      return 0;
    }
    if (j < q->size && q->elements[j].key > e.key) {
      eprintf("Priority queue not a valid heap at ["PF_VTX_T":"PF_VTX_T"]\n",
          e.vtx,e.offset);
      return 0;
    }
    if (k < q->size && q->elements[k].key > e.key) {
      eprintf("Priority queue not a valid heap at ["PF_VTX_T":"PF_VTX_T"]\n",
          e.vtx,e.offset);
      return 0;
    }
  }
  return 1;
}


static void __heap_fix_down(
    vtx_t i, 
    const edge_kv_t val, 
    edge_pq_t * const q)
{
  vtx_t k,j;
  while (1) {
    q->elements[i] = val;
    q->index[val.vtx][val.offset] = i;
    j = __LEFTCHILDINDEX(i);
    k = __RIGHTCHILDINDEX(i);
    if (j < q->size) {
      if (k < q->size && q->elements[i].key < q->elements[k].key &&
          q->elements[j].key <  q->elements[k].key) {
        q->elements[i] = q->elements[k];
        q->index[q->elements[k].vtx][q->elements[k].offset] = i;
        i = k;
      } else if (val.key < q->elements[j].key) {
        q->elements[i] = q->elements[j];
        q->index[q->elements[j].vtx][q->elements[j].offset] = i;
        i = j;
      } else {
        break;
      }
    } else {
      break;
    }
  }
}


static void __heap_fix_up(
    vtx_t i, 
    const edge_kv_t val, 
    edge_pq_t * const q)
{
  vtx_t j = __PARENTINDEX(i);
  while (i > 0 && q->elements[j].key < val.key) {
    q->elements[i] = q->elements[j];
    q->index[q->elements[j].vtx][q->elements[j].offset] = i;
    i = j;
    j = __PARENTINDEX(i);
  }
  q->elements[i] = val;
  q->index[val.vtx][val.offset] = i;
}


static edge_pq_t * __edge_pq_create(
    const size_t size, 
    const size_t nbranches)
{
  edge_pq_t * q;

  q = (edge_pq_t*)malloc(sizeof(edge_pq_t));
  q->maxsize = nbranches;
  q->size = 0;
  q->elements = edge_kv_alloc(size);
  q->index = r_vtx_calloc(nbranches);

  return q;
}


static void __edge_pq_addvertex(
    const vtx_t v, 
    const vtx_t ne, 
    edge_pq_t * const q)
{
  q->index[v] = vtx_alloc(ne);
}


static void __edge_pq_deleteedge(
    const vtx_t v, 
    const vtx_t k, 
    edge_pq_t * const q)
{
  vtx_t i;
  edge_kv_t val;

  i = q->index[v][k];

  DL_ASSERT(i != NULL_VTX,"Attempting to delete unindexed edge\n");

  val = q->elements[--q->size];
  if (val.key > q->elements[i].key) {
    __heap_fix_up(i,val,q);
  } else {
    __heap_fix_down(i,val,q);
  }
  q->index[v][k] = NULL_VTX;

  DL_ASSERT(__check_pq(q),"Bad priority queue after edge deletion\n");
}


static void __edge_pq_deletevertex(
    const vtx_t v, 
    const vtx_t ne,
    edge_pq_t * const q)  
{
  vtx_t i, k;
  for (k=0;k<ne;++k) {
    i = q->index[v][k];
    if (i != NULL_VTX) {
      __edge_pq_deleteedge(v,k,q);
    }
  }
  dl_free(q->index[v]);
  q->index[v] = NULL;
}


static int __edge_pq_push(
    const mod_t mod, 
    const vtx_t v, 
    const vtx_t i,
    edge_pq_t * const q)
{
  vtx_t j;
  edge_kv_t kv;
  kv.key = mod;
  kv.vtx = v;
  kv.offset = i;
  j = (vtx_t)q->size++;
  DL_ASSERT(q->index[v] != NULL,"Trying to push an edge to a deleted "
      "vertex\n");
  __heap_fix_up(j,kv,q);

  DL_ASSERT(__check_pq(q),"Bad priority queue after push\n");

  return 1;
}


static edge_kv_t __edge_pq_pop(
    edge_pq_t * q)
{
  vtx_t i = 0;
  edge_kv_t kv = q->elements[i];
  edge_kv_t val = q->elements[--q->size];
  __heap_fix_down(i,val,q);
  q->index[kv.vtx][kv.offset] = NULL_VTX;

  DL_ASSERT(__check_pq(q),"Bad priority queue after pop\n");

  return kv;
}


static void __edge_pq_move(
    const vtx_t v, 
    const vtx_t i, 
    const vtx_t j, 
    edge_pq_t * const q)
{
  vtx_t k = (vtx_t)q->index[v][i];

  DL_ASSERT(k != NULL_VTX,"Attempting to move "PF_VTX_T" to "PF_VTX_T" for "
      "vertex "PF_VTX_T"\n",i,j,v);

  q->index[v][j] = k;
  q->elements[k].offset = j;
}


static void __edge_pq_free(
    edge_pq_t * q) 
{
  size_t i;
  dl_free(q->elements);
  for (i=0;i<q->maxsize;++i) {
    if (q->index[i]) {
      dl_free(q->index[i]);
    }
  }
  dl_free(q->index);
  dl_free(q);
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


clustering_t * generate_newman_clustering(
    const graph_t * const graph)
{
  vtx_t i, v, u, k, l, mynvtxs, nv, nne, nm, cm;
  mod_t mod, curmod, maxmod;
  adj_t j;
  tid_t t, o;
  cid_t nclusters, c;
  cid_t * ctable;
  vtx_t ** lbl, * sadj, * htable, * mb, *ma;
  edge_t ** edges, * edge, * newedges;
  clustering_t * clustering;
  edge_pq_t * q;
  cid_t ** where, * clbl;
  wgt_t * deg, * ideg;
  vtx_djset_t * set, * maxset;
  edge_kv_t kv;

  const tid_t nthreads = graph->npar;
  const vtx_t nvtxs = graph->nvtxs;
  const adj_t nedges = graph->nedges;
  const twgt_t gadjwgt = graph->gadjwgt;
  const adj_t * const * const xadj = (const adj_t * const *)graph->xadj;
  const vtx_t * const * const adjncy = (const vtx_t * const *)graph->adjncy;
  const wgt_t * const * const adjwgt = (const wgt_t * const *)graph->adjwgt;
  const wgt_t * const * const iadjwgt = (const wgt_t * const *)graph->iadjwgt;
  const wgt_t * const * const eadjwgt = (const wgt_t * const *)graph->eadjwgt;

  lbl = r_vtx_dalloc(graph->mynvtxs,sizeof(vtx_t),nthreads);
  edges = r_edge_alloc(nvtxs);
  sadj = vtx_calloc(nvtxs);
  deg = wgt_alloc(nvtxs);
  ideg = wgt_alloc(nvtxs);
  htable = vtx_init_alloc(NULL_VTX,nvtxs);

  q = __edge_pq_create(nedges,nvtxs);

  maxmod = 0;
  nv = 0;
  /* build up the vertices */
  for (t=0;t<nthreads;++t) {
    mynvtxs = graph->mynvtxs[t];
    for (i=0;i<mynvtxs;++i) {
      lbl[t][i] = nv;
      edges[nv] = edge_alloc(xadj[t][i+1]-xadj[t][i]);
      __edge_pq_addvertex(nv,(vtx_t)(xadj[t][i+1]-xadj[t][i]),q);
      ideg[nv] = iadjwgt[t][i];
      deg[nv] = eadjwgt[t][i] + ideg[nv];
      ++nv;
    } 
  }

  /* build the edges and insert them in the priority queue */
  for (t=0;t<nthreads;++t) {
    mynvtxs = graph->mynvtxs[t];
    for (i=0;i<mynvtxs;++i) {
      v = lbl[t][i];
      for (j=xadj[t][i];j<xadj[t][i+1];++j) {
        k = adjncy[t][j];
        if (k < mynvtxs) {
          o = t;
          l = k;
        } else {
          o = gvtx_to_tid(k,graph->dist);
          l = gvtx_to_lvtx(k,graph->dist);
        }
        u = lbl[o][l];
        if (v > u) {
          /* insert edges into the priority queue */
          mod = (mod_t)(adjwgt[t][j] - ((deg[v]*deg[u])/gadjwgt));
          __edge_pq_push(mod,v,sadj[v],q);
          __edge_pq_push(mod,u,sadj[u],q);
          /* handle v's edge */
          edges[v][sadj[v]].v = v;
          edges[v][sadj[v]].u = u;
          edges[v][sadj[v]].radj = sadj[u];
          edges[v][sadj[v]].wgt = adjwgt[t][j];
          /* handle u's edge */
          edges[u][sadj[u]].v = u;
          edges[u][sadj[u]].u = v;
          edges[u][sadj[u]].radj = sadj[v];
          edges[u][sadj[u]].wgt = adjwgt[t][j];
          ++sadj[v];
          ++sadj[u];
        }
      }
    }
  }

  DL_ASSERT(__check_edges((const edge_t * const *)edges,sadj,nvtxs),"Bad edge "
      "array after construction\n");

  set = vtx_djset_create(0,nvtxs);
  maxset = vtx_djset_create(0,nvtxs);
  ma = vtx_alloc(nvtxs);
  mb = vtx_alloc(nvtxs);

  curmod = 0.0;
  nm = cm = 0;
  nclusters = graph->nvtxs;
  /* collapse every edge keeping track of the maximum modularity */
  while (q->size > 0) {
    /* retrieve the collapsing edge */
    kv = __edge_pq_pop(q);
    mod = kv.key;
    edge = edges[kv.vtx]+kv.offset;
    v = edge->v;
    u = edge->u;
    DL_ASSERT(u != v,"Encountered a loop edge!\n");

    /* update the modularity and stuff */
    curmod += 2*mod; 
    ideg[v] += 2*edge->wgt + ideg[u];
    deg[v] += deg[u];

    /* update the number of clusters and merge the vertices */
    ma[nm] = v;
    mb[nm] = u;
    ++nm;
    vtx_djset_join(v,u,set);
    --nclusters;

    /* delete edges from priority queue */
    __edge_pq_deletevertex(v,sadj[v],q);
    __edge_pq_deletevertex(u,sadj[u],q);

    /* update edges */
    nne = 0;
    newedges = edge_alloc(sadj[v]+sadj[u]); 
    for (i=0;i<sadj[v];++i) {
      k = edges[v][i].u;
      if (k != u) {
        /* remove remote edge ends from priority queue */
        __edge_pq_deleteedge(k,(vtx_t)edges[v][i].radj,q);
        htable[k] = nne;
        newedges[nne] = edges[v][i];
        /* update radj on remote vertices */
        edges[k][edges[v][i].radj].radj = (adj_t)nne;
        ++nne;
      }
    }
    /* replace the old edges, so updates can be performed on the new ones */
    dl_free(edges[v]);
    edges[v] = newedges;
    for (i=0;i<sadj[u];++i) {
      k = edges[u][i].u;
      if (k != v) {
        /* remove remote edge ends from priority queue */
        __edge_pq_deleteedge(k,(vtx_t)edges[u][i].radj,q);
        if (htable[k] == NULL_VTX) {
          htable[k] = nne;
          newedges[nne] = edges[u][i];
          newedges[nne].v = v;
          /* update radj and u on remote vertices */
          edges[k][newedges[nne].radj].radj = nne;
          edges[k][newedges[nne].radj].u = v;
          ++nne;
        } else {
          newedges[htable[k]].wgt += edges[u][i].wgt;
          edges[k][newedges[htable[k]].radj].wgt += edges[u][i].wgt;
          --sadj[k];
          if (edges[u][i].radj != sadj[k]) {
            if (q->index[k][sadj[k]] != NULL_VTX) {
              /* update priority queue index if it exists */
              __edge_pq_move(k,sadj[k],(vtx_t)edges[u][i].radj,q);
            }
            /* remove the edge to u and update it's radj */
            edges[k][edges[u][i].radj] = edges[k][sadj[k]];
            edge = edges[k]+edges[u][i].radj;
            edges[edge->u][edge->radj].radj = edges[u][i].radj; 
          }
        }
      }
    }
    sadj[v] = nne;
    sadj[u] = 0;

    /* update merge priorities of edges */
    __edge_pq_addvertex(v,nne,q);
    for (i=0;i<nne;++i) {
      /* only update the priorities in one direction */
      k = edges[v][i].u;
      mod = (mod_t)(edges[v][i].wgt - ((deg[v]*deg[k])/gadjwgt));
      __edge_pq_push(mod,v,i,q);
      __edge_pq_push(mod,k,(vtx_t)edges[v][i].radj,q);

      /* clear the htable */
      htable[k] = NULL_VTX;
    }

    /* check if we've reached a new maximum modularity */
    if (curmod >= maxmod) {
      maxmod = curmod;
      /* update the maxset to reflect the new maxmod */
      for (;cm<nm;++cm) {
        vtx_djset_join(ma[cm],mb[cm],maxset);
      }
    }

    DL_ASSERT(__check_edges((const edge_t * const *)edges,sadj,nvtxs),"Bad "
        "edge array after collapsing "PF_VTX_T" and "PF_VTX_T" (iteration "
        PF_VTX_T")\n",v,u,nvtxs-(vtx_t)nclusters);
  }

  /* free our intermediate data vectors */
  dl_free(ma);
  dl_free(mb);
  dl_free(htable);
  dl_free(deg);
  dl_free(ideg);
  __edge_pq_free(q);
  vtx_djset_free(set);
  dl_free(sadj);
  r_edge_free(edges,nvtxs);

  /* use maxset to create the where vector */
  ctable = cid_init_alloc(NULL_CID,nvtxs);
  clbl = cid_alloc(maxset->nsets);
  nclusters = 0;
  where = r_cid_dalloc(graph->mynvtxs,sizeof(vtx_t),nthreads);
  for (t=0;t<nthreads;++t) {
    mynvtxs = graph->mynvtxs[t];
    for (i=0;i<mynvtxs;++i) {
      v = lbl[t][i];
      l = vtx_djset_find(v,maxset);
      if ((c = ctable[l]) == NULL_CID) {
        c = ctable[l] = nclusters++; 
      }
      where[t][i] = c;
    }
  }

  /* free last few things */
  dl_free(ctable);
  dl_free(clbl);
  vtx_djset_free(maxset);
  r_vtx_free(lbl,nthreads);

  /* create the clustering */
  clustering = setup_clustering(nclusters,where,graph);

  return clustering;
}




#endif
