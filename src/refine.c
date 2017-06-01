/**
 * @file refine.c
 * @brief Functions for refinement operations
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-10-08
 */




#ifndef NERSTRAND_REFINE_C
#define NERSTRAND_REFINE_C




#include "cluster.h"
#include "refine.h"




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLMEM_PREFIX move
#define DLMEM_TYPE_T move_t
#define DLMEM_STATIC
#define DLMEM_DLTYPE DLTYPE_STRUCT
#include "dlmem_headers.h"
#undef DLMEM_DLTYPE
#undef DLMEM_STATIC
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX


#define DLMEM_PREFIX update
#define DLMEM_TYPE_T update_t
#define DLMEM_STATIC
#define DLMEM_DLTYPE DLTYPE_STRUCT
#include "dlmem_headers.h"
#undef DLMEM_DLTYPE
#undef DLMEM_STATIC
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX


#define DLBUFFER_PREFIX move
#define DLBUFFER_TYPE_T move_t
#define DLBUFFER_STATIC
#include "dlbuffer_headers.h"
#undef DLBUFFER_STATIC
#undef DLBUFFER_TYPE_T
#undef DLBUFFER_PREFIX


#define DLBUFFER_PREFIX update
#define DLBUFFER_TYPE_T update_t
#define DLBUFFER_STATIC
#include "dlbuffer_headers.h"
#undef DLBUFFER_STATIC
#undef DLBUFFER_TYPE_T
#undef DLBUFFER_PREFIX


#define DLCB_PREFIX update
#define DLCB_TYPE_T update_t
#define DLCB_BUFFER_TYPE_T update_buffer_t
#define DLCB_BUFFER_PREFIX update_buffer
#define DLCB_STATIC
#include "dlcb_headers.h"
#undef DLCB_BUFFER_TYPE_T
#undef DLCB_BUFFER_PREFIX
#undef DLCB_STATIC
#undef DLCB_TYPE_T
#undef DLCB_PREFIX




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


static size_t const CLUSTER_CHUNK = 32; 




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static inline int __valid_move(
    cid_t const from, 
    cid_t const to, 
    int p)
{
  if (p == 0) {
    return from > to;
  } else {
    return from < to;
  }
}


static inline int __update_vertex(
    vtx_t const v, 
    vtx_t const g,
    cid_t const home, 
    wgt_t const ewgt, 
    cid_t const from, 
    cid_t const to, 
    wgt_t const invgadjwgt, 
    ucinfo_t * const ucinfo, 
    cluster_t * const clusters,
    vtx_iset_t * const mk, 
    rv_pq_t * const queue, 
    adj_t const size,
    wgt_t const iwgt)
{
  cid_t ito,ifrom,c;
  real_t priority;
  int eg;

  nbrinfo_t * const myinfo = ucinfo->nbrinfo+v;
  vtx_iset_t * const bnd = ucinfo->bnd;

  /* check if we need a spot for k's nbrinfo */
  if (myinfo->nbrstart == NO_NBRS) {
    DL_ASSERT_EQUALS(myinfo->nnbrs,0,PF_CID_T);
    myinfo->nbrstart = next_nbr_adj(ucinfo,size);
  }

  /* point nbradj and nbrwgt at k's index */
  cid_t * const nbradj = ucinfo->adjncy + myinfo->nbrstart;
  wgt_t * const nbrwgt = ucinfo->adjwgt + myinfo->nbrstart;
  vtx_t * const nbrnum = ucinfo->adjnum + myinfo->nbrstart;

  /* find the index of to and from clusters */
  ito = NULL_CID;
  ifrom = NULL_CID;
  for (c=0;c<myinfo->nnbrs;++c) {
    if (nbradj[c] == to) {
      ito = c;
    } else if (nbradj[c] == from) {
      ifrom = c;
    }
  }

  /* update nbrinfo of k */
  if (home == to) {
    DL_ASSERT(ifrom != NULL_CID,"ifrom did not get set when the neighbor is "
        "in the 'to' partition\n");
    /* ito will be unset */
    myinfo->iw += ewgt;
    ++myinfo->in;
    myinfo->ew -= ewgt;
    nbrwgt[ifrom] -= ewgt;
    --nbrnum[ifrom];
    if (nbrnum[ifrom] == 0) {
      --myinfo->nnbrs;
      nbradj[ifrom] = nbradj[myinfo->nnbrs];
      nbrwgt[ifrom] = nbrwgt[myinfo->nnbrs];
      nbrnum[ifrom] = nbrnum[myinfo->nnbrs];
    }
    eg = 1;
  } else if (home == from) {
    /* ifrom will be unset */
    myinfo->ew += ewgt;
    myinfo->iw -= ewgt;
    --myinfo->in;
    if (ito == NULL_CID) { /* no prior connection to 'to' */
      ito = myinfo->nnbrs++;
      nbradj[ito] = to;
      nbrwgt[ito] = ewgt;
      nbrnum[ito] = 1;
    } else {
      nbrwgt[ito] += ewgt;
      ++nbrnum[ito];
    }
    eg = -1;
  } else {
    DL_ASSERT(ifrom != NULL_CID,"Vertex "PF_VTX_T" in "PF_CID_T" did not get "
          "ifrom set when the neighbor (ewgt = "PF_WGT_T")\nis neither in the "
          "'to' ("PF_CID_T") or 'from' ("PF_CID_T") partitions\n",v,home,ewgt,
          to,from);

    /* ifrom will be set, ito maybe set */
    nbrwgt[ifrom] -= ewgt;
    --nbrnum[ifrom];
    if (nbrnum[ifrom] == 0) {
      --myinfo->nnbrs;
      nbradj[ifrom] = nbradj[myinfo->nnbrs];
      nbrwgt[ifrom] = nbrwgt[myinfo->nnbrs];
      nbrnum[ifrom] = nbrnum[myinfo->nnbrs];
      if (ito == myinfo->nnbrs) {
        ito = ifrom;
      }
    }
    if (ito == NULL_CID) { /* no prior connection to 'to' */
      ito = myinfo->nnbrs++;
      nbradj[ito] = to;
      nbrwgt[ito] = ewgt;
      nbrnum[ito] = 1;
    } else {
      nbrwgt[ito] += ewgt;
      ++nbrnum[ito];
    }
    eg = 0;
  }

  if (myinfo->nnbrs == 0) {
    myinfo->ew = 0;
  }

  /* update boundary and queue status of k */
  priority = calc_mod_difference(iwgt,myinfo->iw,
      myinfo->ew+myinfo->iw,clusters[home].eadjwgt+clusters[home].viadjwgt+
      clusters[home].iadjwgt,invgadjwgt);
  if (vtx_iset_contains(v,bnd)) {
    if (!is_bnd(myinfo->iw,myinfo->ew,iwgt)) {
      vtx_iset_remove(v,bnd);
      if (queue && rv_pq_contains(g,queue)) {
        rv_pq_remove(g,queue);
      }
    } else if (queue && rv_pq_contains(g,queue)) {
      rv_pq_update(priority,g,queue);
    }
  } else {
    if (is_bnd(myinfo->iw,myinfo->ew,iwgt)) {
      vtx_iset_add(v,ucinfo->bnd);
      /* only add it if it hasn't been in the queue yet */
      if (queue && !vtx_iset_contains(v,mk)) {
        rv_pq_push(priority,g,queue);
      }
    }
  }

  return eg;
}


static inline int __par_update_vertex(
    vtx_t const v, 
    cid_t const home, 
    wgt_t const ewgt, 
    cid_t const from, 
    cid_t const to, 
    wgt_t const invgadjwgt, 
    ucinfo_t * const ucinfo, 
    wgt_t * const ceadjwgt, 
    wgt_t * const ciadjwgt,
    vtx_iset_t * const mk, 
    rv_pq_t * const queue, 
    adj_t const size,
    wgt_t const iwgt)
{
  cid_t ito,ifrom,c;
  real_t priority;
  int eg;

  nbrinfo_t * const myinfo = ucinfo->nbrinfo+v;
  vtx_iset_t * const bnd = ucinfo->bnd;

  /* check if we need a spot for k's nbrinfo */
  if (myinfo->nbrstart == NO_NBRS) {
    DL_ASSERT_EQUALS(myinfo->nnbrs,0,PF_CID_T);
    myinfo->nbrstart = next_nbr_adj(ucinfo,size);
  }

  /* point nbradj and nbrwgt at k's index */
  cid_t * const nbradj = ucinfo->adjncy + myinfo->nbrstart;
  wgt_t * const nbrwgt = ucinfo->adjwgt + myinfo->nbrstart;
  vtx_t * const nbrnum = ucinfo->adjnum + myinfo->nbrstart;

  /* find the index of to and from clusters */
  ito = NULL_CID;
  ifrom = NULL_CID;
  for (c=0;c<myinfo->nnbrs;++c) {
    if (nbradj[c] == to) {
      ito = c;
    } else if (nbradj[c] == from) {
      ifrom = c;
    }
  }

  /* update nbrinfo of k */
  if (home == to) {
    DL_ASSERT(ifrom != NULL_CID,"ifrom did not get set when the neighbor is "
        "in the 'to' partition\n");
    /* ito will be unset */
    myinfo->iw += ewgt;
    ++myinfo->in;
    myinfo->ew -= ewgt;
    nbrwgt[ifrom] -= ewgt;
    --nbrnum[ifrom];
    if (nbrnum[ifrom] == 0) {
      --myinfo->nnbrs;
      nbradj[ifrom] = nbradj[myinfo->nnbrs];
      nbrwgt[ifrom] = nbrwgt[myinfo->nnbrs];
      nbrnum[ifrom] = nbrnum[myinfo->nnbrs];
    }
    eg = 1;
    if (ceadjwgt != NULL) {
      ceadjwgt[to] -= ewgt / 2.0;
      ceadjwgt[from] -= ewgt / 2.0;
    }
    if (ciadjwgt != NULL) {
      ciadjwgt[to] += ewgt; /* internal edge weights are counted twice */
    }
  } else if (home == from) {
    /* ifrom will be unset */
    myinfo->ew += ewgt;
    myinfo->iw -= ewgt;
    --myinfo->in;
    if (ito == NULL_CID) { /* no prior connection to 'to' */
      ito = myinfo->nnbrs++;
      nbradj[ito] = to;
      nbrwgt[ito] = ewgt;
      nbrnum[ito] = 1;
    } else {
      nbrwgt[ito] += ewgt;
      ++nbrnum[ito];
    }
    eg = -1;
    if (ceadjwgt != NULL) {
      ceadjwgt[from] += ewgt / 2.0;
      ceadjwgt[to] += ewgt / 2.0;
    }
    if (ciadjwgt != NULL) {
      ciadjwgt[from] -= ewgt; /* internal edge weights are counted twice */
    }
  } else {
    DL_ASSERT(ifrom != NULL_CID,"Vertex "PF_VTX_T" in "PF_CID_T" did not get "
          "ifrom set when the neighbor (ewgt = "PF_WGT_T")\nis neither in the "
          "'to' ("PF_CID_T") or 'from' ("PF_CID_T") partitions\n",v,home,ewgt,
          to,from);

    /* ifrom will be set, ito maybe set */
    nbrwgt[ifrom] -= ewgt;
    --nbrnum[ifrom];
    if (nbrnum[ifrom] == 0) {
      --myinfo->nnbrs;
      nbradj[ifrom] = nbradj[myinfo->nnbrs];
      nbrwgt[ifrom] = nbrwgt[myinfo->nnbrs];
      nbrnum[ifrom] = nbrnum[myinfo->nnbrs];
      if (ito == myinfo->nnbrs) {
        ito = ifrom;
      }
    }
    if (ito == NULL_CID) { /* no prior connection to 'to' */
      ito = myinfo->nnbrs++;
      nbradj[ito] = to;
      nbrwgt[ito] = ewgt;
      nbrnum[ito] = 1;
    } else {
      nbrwgt[ito] += ewgt;
      ++nbrnum[ito];
    }
    eg = 0;
    if (ceadjwgt != NULL) {
      ceadjwgt[from] -= ewgt / 2.0;
      ceadjwgt[to] += ewgt / 2.0;
    }
  }

  if (myinfo->nnbrs == 0) {
    myinfo->ew = 0;
  }

  /* update boundary and queue status of k */
  priority = calc_mod_difference(iwgt,myinfo->iw,
      myinfo->ew+myinfo->iw,ceadjwgt[home]+ciadjwgt[home],invgadjwgt);
  if (vtx_iset_contains(v,bnd)) {
    if (!is_bnd(myinfo->iw,myinfo->ew,iwgt)) {
      vtx_iset_remove(v,bnd);
      if (queue && rv_pq_contains(v,queue)) {
        rv_pq_remove(v,queue);
      }
    } else if (queue && rv_pq_contains(v,queue)) {
      rv_pq_update(priority,v,queue);
    }
  } else {
    if (is_bnd(myinfo->iw,myinfo->ew,iwgt)) {
      vtx_iset_add(v,ucinfo->bnd);
      /* only add it if it hasn't been in the queue yet */
      if (queue && !vtx_iset_contains(v,mk)) {
        rv_pq_push(priority,v,queue);
      }
    }
  }

  return eg;
}


static inline void __move_vertex(
    vtx_t const v, 
    cid_t const from, 
    cid_t const to, 
    cid_t * const where,
    ucinfo_t * const ucinfo, 
    vtx_iset_t * const mk, 
    adj_t const size,
    wgt_t const iwgt)
{
  wgt_t oid;
  vtx_t oin;
  cid_t c;
  cid_t * nbradj;
  wgt_t * nbrwgt;
  vtx_t * nbrnum;

  /* expose uncoarsening stuff */
  nbrinfo_t * const myinfo = ucinfo->nbrinfo+v;
  vtx_iset_t * const bnd = ucinfo->bnd;

  /* if this guy doesn't have neighbors, he will */
  if (myinfo->nbrstart == NO_NBRS) {
    myinfo->nbrstart = next_nbr_adj(ucinfo,size);
  }

  /* point nbradj and nbrwgt at v's index */
  nbradj = ucinfo->adjncy + myinfo->nbrstart;
  nbrwgt = ucinfo->adjwgt + myinfo->nbrstart;
  nbrnum = ucinfo->adjnum + myinfo->nbrstart;

  /* find the cluster location */
  for (c=0;c<myinfo->nnbrs;++c) {
    if (nbradj[c] == to) {
      break;
    }
  }
  if (myinfo->in > 0) {
    /* if this vertex has internal connections */
    if (c==myinfo->nnbrs) {
      /* if v is not connected to 'to', create a space for the missing edge to
       be handled below */
      nbrwgt[myinfo->nnbrs] = 0;
      nbrnum[myinfo->nnbrs] = 0;
      ++myinfo->nnbrs;
    }
    DL_ASSERT(myinfo->nnbrs <= size,"Number of neighbors ("PF_CID_T") is "
        "greater than size ("PF_ADJ_T")\n",myinfo->nnbrs,size);

    /* update v's info */
    oid = myinfo->iw;
    oin = myinfo->in;
    myinfo->in = nbrnum[c];
    myinfo->iw = nbrwgt[c];
    myinfo->ew += oid - nbrwgt[c];
    nbradj[c] = from;
    nbrwgt[c] = oid;
    nbrnum[c] = oin;
    if (nbrnum[c] == 0) {
      --myinfo->nnbrs;
      nbradj[c] = nbradj[myinfo->nnbrs];
      nbrwgt[c] = nbrwgt[myinfo->nnbrs];
      nbrnum[c] = nbrnum[myinfo->nnbrs];
    }
  } else if (c < myinfo->nnbrs) {
    /* otherwise if this vertex is connected to 'to' */
    myinfo->iw = nbrwgt[c];
    myinfo->in = nbrnum[c];
    myinfo->ew -= nbrwgt[c];
    --myinfo->nnbrs;
    nbradj[c] = nbradj[myinfo->nnbrs];
    nbrwgt[c] = nbrwgt[myinfo->nnbrs];
    nbrnum[c] = nbrnum[myinfo->nnbrs];
  }   
  DL_ASSERT_EQUALS(where[v],from,PF_CID_T);
  where[v] = to;

  /* check if v should still be in the bnd */
  if (vtx_iset_contains(v,bnd)) {
    if (!is_bnd(myinfo->iw,myinfo->ew,iwgt)) {
      vtx_iset_remove(v,bnd);
    } 
  } else {
    if (is_bnd(myinfo->iw,myinfo->ew,iwgt)) {
      vtx_iset_add(v,bnd);
    }
  }

  /* mark the vertex as moved */
  if (mk) {
    vtx_iset_add(v,mk);
  }
}




/******************************************************************************
* PRIVATE SERIAL FUNCTIONS ****************************************************
******************************************************************************/


/**
 * @brief Remove empty clusters serially
 *
 * @param clustering The clustering to remove clusters from
 * @param mgraph The clustered mgraph
 *
 * @return The modified clustering
 */
static clustering_t * __remove_empty_clusters(
    clustering_t * const clustering,
    mgraph_t const * const mgraph)
{
  vtx_t i, mynvtxs, myid;
  cid_t me, c, nkept;
  nbrinfo_t myinfo;

  nbrinfo_t * nbrinfo;
  cid_t * adjncy;

  /* graph */
  graph_t const * const graph = mgraph->graph;
  tid_t const nthreads = mgraph->graph->npar;


  /* clustering */
  cid_t const nclusters = clustering->nclusters;
  cid_t * const * const where = (cid_t * const *)clustering->where;
  cluster_t * clusters = clustering->clusters;

  cid_t * rename = cid_alloc(nclusters);
  clustering->clusters = cluster_calloc(nclusters);

  nkept = 0;
  for (c=0;c<nclusters;++c) {
    if (clusters[c].nvtxs > 0) {
      rename[c] = nkept;
      clustering->clusters[nkept++] = clusters[c];
    }
  }

  clustering->nclusters = nkept;

  /* renumber vertex locations */
  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];

    for (i=0;i<mynvtxs;++i) {
      me = where[myid][i];
      DL_ASSERT(clusters[me].nvtxs > 0,"Found vertex "PF_VTX_T" in empty "
          "cluster "PF_CID_T"\n",i,c);
      where[myid][i] = rename[me];
    }
  }

  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    nbrinfo = mgraph->ucinfo[myid]->nbrinfo;
    adjncy = mgraph->ucinfo[myid]->adjncy;

    /* renumber ucinfo locations */
    for (i=0;i<mynvtxs;++i) {
      myinfo = nbrinfo[i];
      if (myinfo.nbrstart != NO_NBRS) {
        for (c=myinfo.nbrstart;c<myinfo.nbrstart+myinfo.nnbrs;++c) {
          me = adjncy[c];
          adjncy[c] = rename[me];
        }
      }
    }
  }

  dl_free(clusters);
  dl_free(rename);

  return clustering;
}


/**
 * @brief Refine the clustering on a graph using greedy refinement to maximize
 * teh modularity of the clustering
 *
 * @param objective The objective contianing the current random seed
 * @param clustering The clustering to refine
 * @param mgraph The clustered graph
 *
 * @return The updated clsutering
 */
static clustering_t * __refine_GREEDY(
    clustering_t * const clustering, 
    mgraph_t const * const mgraph,
    size_t const npasses,
    const int empty,
    refinement_stats_t * const refstats)
{
  size_t pass;
  vtx_t g,v,i,k,nmoved,mynvtxs,tnbnd;
  adj_t j;
  cid_t from, to,c,tmp_to;
  real_t priority;
  wgt_t maxgain, gain,vd,da,db,aconn,bconn,ewgt,oiw,iwgt;
  tid_t myid, o;
  cid_t * nbradj;
  wgt_t * nbrwgt;
  nbrinfo_t * myinfo;
  vtx_iset_t * mybnd;
  refinement_stat_t refstat;

  /* expose graph stuff */
  graph_t const * const graph = mgraph->graph;    
  tid_t const nthreads = graph->npar;
  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy;
  wgt_t const * const * const adjwgt = (wgt_t const * const *)graph->adjwgt;
  wgt_t const * const * const iadjwgt = (wgt_t const * const *)graph->iadjwgt;
  wgt_t const * const * const eadjwgt = (wgt_t const * const *)graph->eadjwgt;
  wgt_t const gadjwgt = graph->gadjwgt;
  wgt_t const invgadjwgt = 1.0/gadjwgt;

  /* expose mgraph stuff */
  ucinfo_t * const * const ucinfo = (ucinfo_t * const *)mgraph->ucinfo;

  /* expose cluster stuff */
  cid_t const nclusters = clustering->nclusters;
  cid_t * const * const where = clustering->where;
  cluster_t * clusters = clustering->clusters;

  if (refstats) {
    init_refinement_stat(&refstat);
  }

  /* bail out if there's nothing to refine */
  if (nclusters <= 1) {
    return clustering;
  }

  rv_pq_t * queue = rv_pq_create(0,max_gvtx(graph));

  vtx_iset_t ** mk = (vtx_iset_t**)malloc(nthreads*sizeof(vtx_iset_t*));
  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    mk[myid] = vtx_iset_create(0,mynvtxs);
  }

  /* out pass loop */
  for (pass=0;pass<npasses;++pass) {
    nmoved = 0;
    /* clear mark */
    for (myid=0;myid<nthreads;++myid) {
      vtx_iset_clear(mk[myid]);
      mybnd = ucinfo[myid]->bnd;
      /* Add boundary vertices to the queue */
      for (i=0;i<(vtx_t)mybnd->size;++i) {
        v = vtx_iset_get(i,mybnd);
        myinfo = ucinfo[myid]->nbrinfo+v;
        DL_ASSERT_EQUALS(vtx_iset_indexof(v,mybnd),i,PF_VTX_T);
        from = where[myid][v];
        priority = calc_mod_difference(iadjwgt ? iadjwgt[myid][v] : 0,
            myinfo->iw, eadjwgt[myid][v], clusters[from].viadjwgt +
            clusters[from].eadjwgt+clusters[from].iadjwgt,invgadjwgt);
        rv_pq_push(priority,lvtx_to_gvtx(v,myid,graph->dist),queue);
      }
    }
    while (queue->size > 0) {
      g = rv_pq_pop(queue);
      myid = gvtx_to_tid(g,graph->dist);
      v = gvtx_to_lvtx(g,graph->dist);

      mynvtxs = graph->mynvtxs[myid];
      mybnd = ucinfo[myid]->bnd;
      myinfo = ucinfo[myid]->nbrinfo+v;

      DL_ASSERT(vtx_iset_contains(v,mybnd),"Pulled a vertex from the "
          "priority queue that is not in the boundary (i"PF_WGT_T",e"
          PF_WGT_T")\n",myinfo->iw,myinfo->ew);

      /* set all variables relevant to v */
      from = where[myid][v];

      /* don't empty a cluster */
      if (!empty && clusters[from].nvtxs == 1) {
        continue;
      }

      nbradj = ucinfo[myid]->adjncy + myinfo->nbrstart;
      nbrwgt = ucinfo[myid]->adjwgt + myinfo->nbrstart;

      to = from;
      maxgain = 0;

      aconn = myinfo->iw;

      if (iadjwgt) {
        vd = (eadjwgt[myid][v]+iadjwgt[myid][v]);
      } else {
        vd = eadjwgt[myid][v];
      }
      da = (clusters[from].viadjwgt + clusters[from].eadjwgt +
          clusters[from].iadjwgt) - vd;
      vd *= invgadjwgt; /* scale m out ahead of time */
      /* decide who v is most connected to */
      for (c=0;c<myinfo->nnbrs;++c) {
        tmp_to = nbradj[c];
        bconn = nbrwgt[c];
        db = clusters[tmp_to].viadjwgt + clusters[tmp_to].eadjwgt +
            clusters[tmp_to].iadjwgt;
        gain = (bconn - aconn) + (vd*(da-db));
        if (maxgain < gain) { 
          maxgain = gain;
          to = nbradj[c];
        }
      }

      /* if we v is in a good spot, leave it there */
      if (to == from) {
        continue;
      }

      /* record the move */
      ++nmoved;

      /* save info */
      oiw = myinfo->iw;

      if (iadjwgt) {
        iwgt = iadjwgt[myid][v];
      } else {
        iwgt = 0;
      }

      /* actually add the vertex */
      __move_vertex(v,from,to,where[myid],ucinfo[myid],mk[myid],
          dl_min(nclusters,xadj[myid][v+1]-xadj[myid][v]),iwgt);

      /* vertex-internal edges only get counted by the owning thread */
      if (iwgt) {
        clusters[from].viadjwgt -= iwgt;
        clusters[to].viadjwgt += iwgt;
      }

      /* vertex and edges */
      --clusters[from].nvtxs;
      ++clusters[to].nvtxs;

      /* internal edges get counted twice -- once by each thread */
      clusters[from].iadjwgt -= oiw*2;
      clusters[to].iadjwgt += myinfo->iw*2;

      /* external edges get counted once -- half by each thread */
      clusters[from].eadjwgt -= (eadjwgt[myid][v]-oiw)-oiw;
      clusters[to].eadjwgt += (eadjwgt[myid][v]-myinfo->iw)-myinfo->iw;

      for (j=xadj[myid][v];j<xadj[myid][v+1];++j) {
        k = adjncy[myid][j];
        if (k < mynvtxs) {
          g = lvtx_to_gvtx(k,myid,graph->dist);
          o = myid;
        } else {
          g = k;
          k = gvtx_to_lvtx(g,graph->dist);
          o = gvtx_to_tid(g,graph->dist);
        }
        if (adjwgt) {
          ewgt = adjwgt[myid][j];
        } else {
          ewgt = 1.0;
        }
        if (iadjwgt) {
          iwgt = iadjwgt[o][k];
        } else {
          iwgt = 0;
        }

        __update_vertex(k,g,where[o][k],ewgt,from,to,invgadjwgt,
            ucinfo[o],clusters,mk[o],queue,
            dl_min(nclusters,xadj[o][k+1]-xadj[o][k]),iwgt);
      }
    }
    
    DL_ASSERT(check_ucinfo(mgraph->graph,
          (ucinfo_t const * const *)mgraph->ucinfo,
          (cid_t const * const *)where,nclusters) == 1, 
        "Bad ucinfo after pass "PF_SIZE_T"\n",pass);

    if (refstats) {
      refstat.nmoves += nmoved;
    }

    if (nmoved == 0) {
      dprintf("No moves made in pass "PF_SIZE_T"\n",pass);
      ++pass;
      break;
    } 
    dprintf("Made "PF_VTX_T" moves in pass "PF_SIZE_T"\n",nmoved,pass);
  }

  if (refstats) {
    tnbnd = 0;
    for (myid=0;myid<nthreads;++myid) {
      tnbnd += ucinfo[myid]->bnd->size; 
    }
    refstat.npasses = pass;
    refstat.nbnd = tnbnd;
    refinement_stats_push(refstat,refstats);
  }

  for (myid=0;myid<nthreads;++myid) {
    vtx_iset_free(mk[myid]);
  }
  dl_free(mk);
  rv_pq_free(queue);

  if (empty) {
    __remove_empty_clusters(clustering,mgraph);
  }

  dprintf("Finished refinement of graph at level "PF_SIZE_T"\n",mgraph->level);

  return clustering;
}


/**
 * @brief Refine the clustering on a graph using iterative refinement to 
 * maximize the modularity of the clustering
 *
 * @param objective The objective contianing the current random seed
 * @param clustering The clustering to refine
 * @param mgraph The clustered graph
 *
 * @return The updated clsutering
 */
static clustering_t * __refine_RANDOM(
    clustering_t * const clustering, 
    mgraph_t const * const mgraph,
    size_t const npasses,
    const int empty,
    unsigned int seed,
    refinement_stats_t * const refstats)
{
  size_t pass;
  vtx_t g,v,i,k,nmoved,mynvtxs,tnbnd,sq;
  adj_t j;
  cid_t from, to,c,tmp_to;
  wgt_t maxgain, gain,vd,da,db,aconn,bconn,ewgt,oiw,iwgt;
  tid_t myid, o;
  cid_t * nbradj;
  wgt_t * nbrwgt;
  nbrinfo_t * myinfo;
  vtx_iset_t * mybnd;
  refinement_stat_t refstat;
  vtx_t * q;

  /* expose graph stuff */
  graph_t const * const graph = mgraph->graph;    
  tid_t const nthreads = graph->npar;
  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy;
  wgt_t const * const * const adjwgt = (wgt_t const * const *)graph->adjwgt;
  wgt_t const * const * const iadjwgt = (wgt_t const * const *)graph->iadjwgt;
  wgt_t const * const * const eadjwgt = (wgt_t const * const *)graph->eadjwgt;
  wgt_t const gadjwgt = graph->gadjwgt;

  /* expose mgraph stuff */
  ucinfo_t * const * const ucinfo = (ucinfo_t * const *)mgraph->ucinfo;

  /* expose cluster stuff */
  cid_t const nclusters = clustering->nclusters;
  cid_t * const * const where = clustering->where;
  cluster_t * clusters = clustering->clusters;

  if (refstats) {
    init_refinement_stat(&refstat);
  }

  /* bail out if there's nothing to refine */
  if (nclusters <= 1) {
    return clustering;
  }

  q = vtx_alloc(max_gvtx(graph)+1);

  /* out pass loop */
  for (pass=0;pass<npasses;++pass) {
    nmoved = 0;
    sq = 0;
    /* clear mark */
    for (myid=0;myid<nthreads;++myid) {
      mybnd = ucinfo[myid]->bnd;
      /* Add boundary vertices to the queue */
      for (i=0;i<(vtx_t)mybnd->size;++i) {
        v = vtx_iset_get(i,mybnd);
        DL_ASSERT_EQUALS(vtx_iset_indexof(v,mybnd),i,PF_VTX_T);
        q[sq++] = lvtx_to_gvtx(v,myid,graph->dist);
      }
    }
    vtx_pseudo_shuffle_r(q,sq/8,sq,&seed);
    while (sq > 0) {
      g = q[--sq];
      myid = gvtx_to_tid(g,graph->dist);
      v = gvtx_to_lvtx(g,graph->dist);

      mynvtxs = graph->mynvtxs[myid];
      mybnd = ucinfo[myid]->bnd;
      myinfo = ucinfo[myid]->nbrinfo+v;

      if (iadjwgt)  {
        iwgt = iadjwgt[myid][v];
      } else {
        iwgt = 0;
      }

      if (!is_bnd(myinfo->iw,myinfo->ew,iwgt)) {
        continue;
      }

      /* set all variables relevant to v */
      from = where[myid][v];

      /* don't empty a cluster in kway clustering */
      if (!empty && clusters[from].nvtxs == 1) {
        continue;
      }

      nbradj = ucinfo[myid]->adjncy + myinfo->nbrstart;
      nbrwgt = ucinfo[myid]->adjwgt + myinfo->nbrstart;

      to = from;
      maxgain = 0;

      aconn = myinfo->iw;

      if (iadjwgt) {
        vd = (eadjwgt[myid][v]+iadjwgt[myid][v]);
      } else {
        vd = eadjwgt[myid][v];
      }
      da = (clusters[from].viadjwgt + clusters[from].eadjwgt +
          clusters[from].iadjwgt) - vd;
      vd /= gadjwgt; /* scale m out ahead of time */
      /* decide who v is most connected to */
      for (c=0;c<myinfo->nnbrs;++c) {
        tmp_to = nbradj[c];
        bconn = nbrwgt[c];
        db = clusters[tmp_to].viadjwgt + clusters[tmp_to].eadjwgt +
            clusters[tmp_to].iadjwgt;
        gain = (bconn - aconn) + (vd*(da-db));
        if (maxgain < gain) { 
          maxgain = gain;
          to = nbradj[c];
        }
      }

      /* if we v is in a good spot, leave it there */
      if (to == from) {
        continue;
      }

      /* record the move */
      ++nmoved;

      /* save info */
      oiw = myinfo->iw;

      if (iadjwgt) {
        iwgt = iadjwgt[myid][v];
      } else {
        iwgt = 0;
      }

      /* actually add the vertex */
      __move_vertex(v,from,to,where[myid],ucinfo[myid],NULL,
          dl_min(nclusters,xadj[myid][v+1]-xadj[myid][v]),iwgt);

      /* vertex-internal edges only get counted by the owning thread */
      if (iwgt) {
        clusters[from].viadjwgt -= iwgt;
        clusters[to].viadjwgt += iwgt;
      }

      /* vertex and edges */
      --clusters[from].nvtxs;
      ++clusters[to].nvtxs;

      /* internal edges get counted twice -- once by each thread */
      clusters[from].iadjwgt -= oiw*2;
      clusters[to].iadjwgt += myinfo->iw*2;

      /* external edges get counted once -- half by each thread */
      clusters[from].eadjwgt -= (eadjwgt[myid][v]-oiw)-oiw;
      clusters[to].eadjwgt += (eadjwgt[myid][v]-myinfo->iw)-myinfo->iw;

      for (j=xadj[myid][v];j<xadj[myid][v+1];++j) {
        k = adjncy[myid][j];
        if (k < mynvtxs) {
          g = lvtx_to_gvtx(k,myid,graph->dist);
          o = myid;
        } else {
          g = k;
          k = gvtx_to_lvtx(g,graph->dist);
          o = gvtx_to_tid(g,graph->dist);
        }
        if (adjwgt) {
          ewgt = adjwgt[myid][j];
        } else {
          ewgt = 1.0;
        }
        if (iadjwgt) {
          iwgt = iadjwgt[o][k];
        } else {
          iwgt = 0;
        }

        __update_vertex(k,g,where[o][k],ewgt,from,to,0, /* no priority */
            ucinfo[o],clusters,NULL,NULL,
            dl_min(nclusters,xadj[o][k+1]-xadj[o][k]),iwgt);

        myinfo = ucinfo[o]->nbrinfo+k;

        /* manually do updates to the priority queue */
        c = where[o][k];
      }
    }

    dprintf("Made "PF_VTX_T" moves in pass "PF_SIZE_T"\n",nmoved,pass);
    
    DL_ASSERT(check_ucinfo(mgraph->graph,
          (ucinfo_t const * const *)mgraph->ucinfo,
          (cid_t const * const *)where,nclusters) == 1, 
        "Bad ucinfo after pass "PF_SIZE_T"\n",pass);

    if (refstats) {
      refstat.nmoves += nmoved;
    }

    if (nmoved == 0) {
      ++pass;
      break;
    } 
  }

  if (refstats) {
    tnbnd = 0;
    for (myid=0;myid<nthreads;++myid) {
      tnbnd += ucinfo[myid]->bnd->size; 
    }
    refstat.npasses = pass;
    refstat.nbnd = tnbnd;
    refinement_stats_push(refstat,refstats);
  }

  dl_free(q);

  if (empty) {
    __remove_empty_clusters(clustering,mgraph);
  }

  dprintf("Finished refinement of graph at level "PF_SIZE_T"\n",mgraph->level);

  return clustering;
}




/******************************************************************************
* PRIVATE PARALLEL FUNCTIONS **************************************************
******************************************************************************/


/**
 * @brief Renumber to clustering and delete empty clusters
 *
 * @param clustering The clustering cleanup
 * @param mgraph THe clustered graph
 *
 * @return The cleaned up clustering  or NULL on error
 */
static clustering_t * __par_remove_empty_clusters(
    clustering_t * const clustering, 
    mgraph_t const * const mgraph,
    dlthread_comm_t const comm)
{
  vtx_t i;
  cid_t me, c, nkept, mystart, myend;
  nbrinfo_t myinfo;
  cid_t * chunk, * rename;

  tid_t const myid = dlthread_get_id(comm);
  tid_t const nthreads = dlthread_get_nthreads(comm);

  /* graph */
  graph_t const * const graph = mgraph->graph;
  vtx_t const mynvtxs = graph->mynvtxs[myid];

  /* ucinfo */
  cid_t * const adjncy = mgraph->ucinfo[myid]->adjncy;
  nbrinfo_t const * const nbrinfo = mgraph->ucinfo[myid]->nbrinfo;

  /* clustering */
  cid_t const nclusters = clustering->nclusters;
  cid_t * const where = clustering->where[myid];
  cluster_t * clusters = clustering->clusters;

  rename = dlthread_get_shmem((sizeof(rename)*nclusters) + \
      (sizeof(chunk)*(nthreads+1)),comm);
  chunk = rename + nclusters;

  mystart = (nclusters / nthreads)*myid + dl_min(myid,nclusters % nthreads);
  myend = mystart + (nclusters/nthreads) + (myid<nclusters%nthreads ? 1 : 0);

  nkept = 0;
  for (c=mystart;c<myend;++c) {
    if (clustering->clusters[c].nvtxs > 0) {
      ++nkept;
    }
  }
  chunk[myid] = nkept;

  dlthread_barrier(comm);

  if (myid == 0) {
    /* won't scale well beyond a few hundred threads */
    chunk[nthreads] = 0;
    cid_prefixsum_exc(chunk,nthreads+1);
    clustering->nclusters = chunk[nthreads];
    clustering->clusters = cluster_calloc(clustering->nclusters);
  }
  dlthread_barrier(comm);

  nkept = chunk[myid];
  for (c=mystart;c<myend;++c) {
    if (clusters[c].nvtxs > 0) {
      rename[c] = nkept;
      clustering->clusters[nkept++] = clusters[c];
    }
  }
  dlthread_barrier(comm);

  /* renumber vertex locations */
  for (i=0;i<mynvtxs;++i) {
    me = where[i];
    where[i] = rename[me];
  }

  /* renumber ucinfo locations */
  for (i=0;i<mynvtxs;++i) {
    myinfo = nbrinfo[i];
    if (myinfo.nbrstart != NO_NBRS) {
      for (c=myinfo.nbrstart;c<myinfo.nbrstart+myinfo.nnbrs;++c) {
        me = adjncy[c];
        adjncy[c] = rename[me];
      }
    }
  }

  dlthread_free_shmem(rename,comm);
  if (myid == 0) {
    dl_free(clusters);
    dprintf("Deleted "PF_CID_T"/"PF_CID_T" clusters\n",
        nclusters-clustering->nclusters,nclusters);
  }

  return clustering;
}


size_t __prg_done;
vtx_t ** __prg_cnvtxs;
wgt_t ** __prg_ceadjwgt;
wgt_t ** __prg_ciadjwgt;
wgt_t ** __prg_cviadjwgt;
/**
 * @brief Refine a clustering on a graph using greedy refinement to maximize
 * the modularity of the clustering -- to be called by each thread in a
 * parallel region
 *
 * @param objective The objective containing refinement parameters
 * @param clustering The clustering to refine
 * @param mgraph The clustered graph
 *
 * @return The updated clustering 
 */
static clustering_t * __par_refine_GREEDY(
    clustering_t * const clustering, 
    mgraph_t const * const mgraph,
    size_t const npasses,
    const int empty,
    refinement_stats_t * const refstats,
    dlthread_comm_t const comm)
{
  int dir;
  size_t pass;
  vtx_t v,i,k,nmoved,nunmoved,tnbnd;
  adj_t j;
  cid_t from, to,c,tmp_to,e, ec, esize;
  real_t priority;
  wgt_t maxgain, gain,vd,da,db,aconn,bconn,ewgt,oiw,iwgt;
  tid_t o, t;
  move_t move;
  update_t * uptr;
  update_t up;
  cid_t * nbradj;
  wgt_t * nbrwgt;
  nbrinfo_t * myinfo;
  update_combuffer_t * combuffer;
  update_buffer_t * updates;
  move_buffer_t * moves;
  refinement_stat_t refstat;

  tid_t const myid = dlthread_get_id(comm);
  tid_t const nthreads = dlthread_get_nthreads(comm);

  if (myid == 0) {
    dprintf("Using Modularity Refinement\n");
  }

  /* expose graph stuff */
  graph_t const * const graph = mgraph->graph;    
  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const * const xadj = graph->xadj[myid];
  vtx_t const * const adjncy = graph->adjncy[myid];
  wgt_t const * const adjwgt = (graph->adjwgt ? graph->adjwgt[myid] : NULL);
  wgt_t const * const iadjwgt = (graph->iadjwgt ? graph->iadjwgt[myid] : NULL);
  wgt_t const * const eadjwgt = (graph->eadjwgt ? graph->eadjwgt[myid] : NULL);
  wgt_t const gadjwgt = graph->gadjwgt;
  wgt_t const invgadjwgt = 1.0/gadjwgt;

  DL_ASSERT_EQUALS(nthreads,graph->npar,PF_TID_T);

  /* expose cluster stuff */
  cid_t const nclusters = clustering->nclusters;
  cid_t * const where = clustering->where[myid];
  cluster_t * clusters = clustering->clusters;

  if (refstats) {
    init_refinement_stat(&refstat);
  }

  /* bail out if there's nothing to refine */
  if (nclusters <= 1) {
    return clustering;
  }

  rv_pq_t * queue = rv_pq_create(0,graph->nvtxs);

  /* expose mgraph stuff */
  ucinfo_t * const ucinfo = mgraph->ucinfo[myid];
  vtx_iset_t * const mybnd = ucinfo->bnd;

  if (myid == 0) {
    __prg_cnvtxs = r_vtx_alloc(nthreads);
    __prg_ciadjwgt = r_wgt_alloc(nthreads);
    __prg_ceadjwgt = r_wgt_alloc(nthreads);
    __prg_cviadjwgt = r_wgt_alloc(nthreads);
  }
  dlthread_barrier(comm);

  /* private cluster vectors */
  wgt_t * cdegs = wgt_alloc(nclusters);
  vtx_t * cnvtxs = __prg_cnvtxs[myid] = vtx_calloc(nclusters);

  /* use one array for reduction purposes later */
  wgt_t * ceadjwgt = __prg_ceadjwgt[myid] = wgt_calloc(nclusters);
  wgt_t * ciadjwgt = __prg_ciadjwgt[myid] = wgt_calloc(nclusters);
  wgt_t * cviadjwgt = __prg_cviadjwgt[myid] = wgt_calloc(nclusters);

  vtx_iset_t * mk = vtx_iset_create(0,mynvtxs);

  /* create com buffers */
  combuffer = update_combuffer_create(comm);

  /* create move buffers */
  moves = move_buffer_create(dl_max(mynvtxs/4,16));

  /* out pass loop */
  for (pass=0;pass<npasses;++pass) {
    nmoved = 0;
    nunmoved = 0;
    /* clear mark */
    vtx_iset_clear(mk);
    for (dir=0;dir<2;++dir) {
      if (myid==0) {
        __prg_done = 0;
      }
      for (c=0;c<nclusters;++c) {
        cdegs[c] = clusters[c].viadjwgt + clusters[c].eadjwgt 
            + clusters[c].iadjwgt;
      }

      /* Add boundary vertices to the queue */
      for (i=0;i<(vtx_t)mybnd->size;++i) {
        v = vtx_iset_get(i,mybnd);
        if (vtx_iset_contains(v,mk)) {
            continue;
        } else {
          myinfo = ucinfo->nbrinfo+v;
          DL_ASSERT_EQUALS(vtx_iset_indexof(v,mybnd),i,PF_VTX_T);
          priority = calc_mod_difference(iadjwgt ? iadjwgt[v] : 0, myinfo->iw,
              eadjwgt[v],cdegs[where[v]],invgadjwgt);
          rv_pq_push(priority,v,queue);
        }
      }

      /* process vertices in my queue */
      while (queue->size > 0) {
        v = rv_pq_pop(queue);
        myinfo = ucinfo->nbrinfo+v;

        DL_ASSERT(vtx_iset_contains(v,mybnd),"Pulled a vertex from the "
            "priority queue that is not in the boundary (i"PF_WGT_T",e"
            PF_WGT_T")\n",myinfo->iw,myinfo->ew);

        /* set all variables relevant to v */
        from = where[v];

        /* don't empty a cluster in kway clustering */
        if (!empty && clusters[from].nvtxs + cnvtxs[from] <= 1) {
          continue;
        }

        nbradj = ucinfo->adjncy + myinfo->nbrstart;
        nbrwgt = ucinfo->adjwgt + myinfo->nbrstart;

        to = from;
        maxgain = 0;

        aconn = myinfo->iw;

        if (iadjwgt) {
          vd = (eadjwgt[v]+iadjwgt[v]);
        } else {
          vd = eadjwgt[v];
        }
        da = (cviadjwgt[from] + ceadjwgt[from] +
            ciadjwgt[from]) + cdegs[from] - vd;
        vd *= invgadjwgt; /* scale m out ahead of time */
        /* decide who v is most connected to */
        for (c=0;c<myinfo->nnbrs;++c) {
          tmp_to = nbradj[c];
          if (__valid_move(from,tmp_to,dir)) {
            bconn = nbrwgt[c];
            db = cviadjwgt[tmp_to] + ceadjwgt[tmp_to] +
                ciadjwgt[tmp_to] + cdegs[tmp_to];
            gain = (bconn - aconn) + (vd*(da-db));
            if (maxgain < gain) { 
              maxgain = gain;
              to = nbradj[c];
            }
          }
        }

        /* if we v is in a good spot, leave it there */
        if (to == from) {
          continue;
        }

        /* record the move */
        ++nmoved;
        move.vid = v;
        move.to = to;
        move.from = from;
        move_buffer_add(move,moves);

        /* save info */
        oiw = myinfo->iw;

        if (iadjwgt) {
          iwgt = iadjwgt[v];
        } else {
          iwgt = 0;
        }

        /* actually add the vertex */
        __move_vertex(v,from,to,where,ucinfo,mk,
            dl_min(nclusters,xadj[v+1]-xadj[v]),iwgt);

        /* vertex and edges */
        --cnvtxs[from];
        ++cnvtxs[to];

        /* vertex-internal edges only get counted by the owning thread */
        if (iadjwgt) {
          cviadjwgt[from] -= iadjwgt[v];
          cviadjwgt[to] += iadjwgt[v];
        }

        /* internal edges get counted twice -- once by each thread */
        ciadjwgt[from] -= oiw;
        ciadjwgt[to] += myinfo->iw;

        /* external edges get counted once -- half by each thread */
        ceadjwgt[from] -= ((eadjwgt[v]-oiw)-oiw)/2.0;
        ceadjwgt[to] += ((eadjwgt[v]-myinfo->iw)-myinfo->iw)/2.0;

        /* update neighbors I own */
        for (j=xadj[v];j<xadj[v+1];++j) {
          k = adjncy[j];
          if (k < mynvtxs) {
            if (adjwgt) {
              ewgt = adjwgt[j];
            } else {
              ewgt = 1.0;
            }
            c = where[k];
            if (iadjwgt) {
              iwgt = iadjwgt[k];
            } else {
              iwgt = 0;
            }
            __par_update_vertex(k,where[k],ewgt,from,to,invgadjwgt,
                ucinfo,ceadjwgt,ciadjwgt,mk,queue,
                dl_min(nclusters,xadj[k+1]-xadj[k]),iwgt);
          } else {
            /* tell the other thread about it once we decide to keep it */
          }
        }
      }
      rv_pq_clear(queue); 

      /* need to remove this to be compatible with pthreads */
      dlthread_atomic_add(__prg_done,comm);

      dlthread_barrier(comm);

      for (e=myid;e<(nclusters/CLUSTER_CHUNK)+1;e+=nthreads) {
        ec = e*CLUSTER_CHUNK;
        esize = dl_min(CLUSTER_CHUNK,nclusters-ec);
        /* cnvtxs */
        for (t=(myid+1)%nthreads;t!=myid;t=(t+1)%nthreads) {
          for (c=0;c<esize;++c) {
            cnvtxs[c+ec] += __prg_cnvtxs[t][c+ec];
          }
        }
        for (c=0;c<esize;++c) {
          clusters[c+ec].nvtxs += cnvtxs[c+ec];
        }
        /* eadjwgt */
        for (t=(myid+1)%nthreads;t!=myid;t=(t+1)%nthreads) {
          for (c=0;c<esize;++c) {
            ceadjwgt[c+ec] += __prg_ceadjwgt[t][c+ec];
          }
        }
        for (c=0;c<esize;++c) {
          clusters[c+ec].eadjwgt += ceadjwgt[c+ec];
        }
        /* iadjwgt */
        for (t=(myid+1)%nthreads;t!=myid;t=(t+1)%nthreads) {
          for (c=0;c<esize;++c) {
            ciadjwgt[c+ec] += __prg_ciadjwgt[t][c+ec];
          }
        }
        for (c=0;c<esize;++c) {
          clusters[c+ec].iadjwgt += ciadjwgt[c+ec];
        }
        /* viadjwgt */
        for (t=(myid+1)%nthreads;t!=myid;t=(t+1)%nthreads) {
          for (c=0;c<esize;++c) {
            cviadjwgt[c+ec] += __prg_cviadjwgt[t][c+ec];
          }
        }
        for (c=0;c<esize;++c) {
          clusters[c+ec].viadjwgt += cviadjwgt[c+ec];
        }
      }

      dlthread_barrier(comm);

      vtx_set(cnvtxs,0,nclusters);
      wgt_set(ceadjwgt,0,nclusters);
      wgt_set(ciadjwgt,0,nclusters);
      wgt_set(cviadjwgt,0,nclusters);

      for (c=0;c<nclusters;++c) {
        cdegs[c] = clusters[c].viadjwgt + clusters[c].eadjwgt 
            + clusters[c].iadjwgt;
      }

      /* undo bad moves */
      for (i=moves->size;i>0;) {
        move = moves->elements[--i];
        v = move.vid;
        to = move.to;
        from = move.from;

        myinfo = ucinfo->nbrinfo+v;
        nbradj = ucinfo->adjncy + myinfo->nbrstart;
        nbrwgt = ucinfo->adjwgt + myinfo->nbrstart;

        /* calculate the gain of undoing the move */
        aconn = myinfo->iw;
        if (iadjwgt) {
          vd = (eadjwgt[v]+iadjwgt[v]);
        } else {
          vd = eadjwgt[v];
        }
        da = (cviadjwgt[to] + ceadjwgt[to] +
            ciadjwgt[to]) + cdegs[to] - vd;
        vd /= gadjwgt; /* scale m out ahead of time */
        for (c=0;c<myinfo->nnbrs;++c) {
          if (nbradj[c] == from) {
            break;
          }
        }
        if (c < myinfo->nnbrs) {
          bconn = nbrwgt[c];
        } else {
          bconn = 0;
        }
        db = cviadjwgt[from] + ceadjwgt[from] +
            ciadjwgt[from] + cdegs[from];
        gain = (bconn - aconn) + (vd*(da-db));

        if (gain > 0 || (!empty && clusters[from].nvtxs + cnvtxs[from] == 0)) {
          /* undo move */
          --nmoved;
          ++nunmoved;

          /* save info */
          oiw = myinfo->iw;

          if (iadjwgt) {
            iwgt = iadjwgt[v];
          } else {
            iwgt = 0;
          }

          /* actually add the vertex */
          __move_vertex(v,to,from,where,ucinfo,mk,
              dl_min(nclusters,xadj[v+1]-xadj[v]),iwgt);

          /* vertex and edges */
          --cnvtxs[to];
          ++cnvtxs[from];

          /* vertex-internal edges only get counted by the owning thread */
          if (iadjwgt) {
            cviadjwgt[to] -= iadjwgt[v];
            cviadjwgt[from] += iadjwgt[v];
          }

          /* internal edges get counted twice -- once by each thread */
          ciadjwgt[to] -= oiw;
          ciadjwgt[from] += myinfo->iw;

          /* external edges get counted once -- half by each thread */
          ceadjwgt[to] -= (((eadjwgt[v]-oiw)-oiw)/2.0);
          ceadjwgt[from] += (((eadjwgt[v]-myinfo->iw)-myinfo->iw)/2.0);

          /* update neighbors I own */
          for (j=xadj[v];j<xadj[v+1];++j) {
            k = adjncy[j];
            if (k < mynvtxs) {
              if (adjwgt) {
                ewgt = adjwgt[j];
              } else {
                ewgt = 1.0;
              }

              if (iadjwgt) {
                iwgt = iadjwgt[k];
              } else {
                iwgt = 0;
              }

              __par_update_vertex(k,where[k],ewgt,to,from,
                  invgadjwgt,ucinfo,ceadjwgt,ciadjwgt,mk,NULL,
                  dl_min(nclusters,xadj[k+1]-xadj[k]),iwgt);
            } else {
              /* The other thread never needs to know */
            }
          }
        } else {
          /* tell other threads since we're keeping the move */
          for (j=xadj[v];j<xadj[v+1];++j) {
            k = adjncy[j];
            if (k >= mynvtxs) {
              /* tell the other thread about it later */
              o = gvtx_to_tid(k,graph->dist);
              k = gvtx_to_lvtx(k,graph->dist);
              up.vid = k;
              if (adjwgt) {
                up.wgt = adjwgt[j];
              } else {
                up.wgt = 1.0f;
              }
              up.to = to;
              up.from = from;
              update_combuffer_add(o,up,combuffer);
            }
          }
        }
      }
      move_buffer_clear(moves);

      /* send updates */
      update_combuffer_send(combuffer);
        
      /* process updates */
      for (o=(myid+1)%nthreads;o!=myid;o=(o+1)%nthreads) {
        updates = update_combuffer_get(o,combuffer);
        uptr = updates->elements;
        for (i=0;i<(vtx_t)updates->size;++i) {
          k = uptr[i].vid;
          ewgt = uptr[i].wgt;
          to = uptr[i].to;
          from = uptr[i].from;

          if (iadjwgt) {
            iwgt = iadjwgt[k];
          } else {
            iwgt = 0;
          }

          __par_update_vertex(k,where[k],ewgt,from,to,invgadjwgt,ucinfo,
              ceadjwgt,ciadjwgt,mk,NULL,
              dl_min(nclusters,xadj[k+1]-xadj[k]),iwgt);
        }
      }
      update_combuffer_clear(combuffer);

      dlthread_barrier(comm);
      for (e=myid;e<(nclusters/CLUSTER_CHUNK)+1;e+=nthreads) {
        ec = e*CLUSTER_CHUNK;
        esize = dl_min(CLUSTER_CHUNK,nclusters-ec);
        /* cnvtxs */
        for (t=(myid+1)%nthreads;t!=myid;t=(t+1)%nthreads) {
          for (c=0;c<esize;++c) {
            cnvtxs[c+ec] += __prg_cnvtxs[t][c+ec];
          }
        }
        for (c=0;c<esize;++c) {
          clusters[c+ec].nvtxs += cnvtxs[c+ec];
        }
        /* eadjwgt */
        for (t=(myid+1)%nthreads;t!=myid;t=(t+1)%nthreads) {
          for (c=0;c<esize;++c) {
            ceadjwgt[c+ec] += __prg_ceadjwgt[t][c+ec];
          }
        }
        for (c=0;c<esize;++c) {
          clusters[c+ec].eadjwgt += ceadjwgt[c+ec];
        }
        /* iadjwgt */
        for (t=(myid+1)%nthreads;t!=myid;t=(t+1)%nthreads) {
          for (c=0;c<esize;++c) {
            ciadjwgt[c+ec] += __prg_ciadjwgt[t][c+ec];
          }
        }
        for (c=0;c<esize;++c) {
          clusters[c+ec].iadjwgt += ciadjwgt[c+ec];
        }
        /* viadjwgt */
        for (t=(myid+1)%nthreads;t!=myid;t=(t+1)%nthreads) {
          for (c=0;c<esize;++c) {
            cviadjwgt[c+ec] += __prg_cviadjwgt[t][c+ec];
          }
        }
        for (c=0;c<esize;++c) {
          clusters[c+ec].viadjwgt += cviadjwgt[c+ec];
        }
      }

      dlthread_barrier(comm);

      vtx_set(cnvtxs,0,nclusters);
      wgt_set(ceadjwgt,0,nclusters);
      wgt_set(ciadjwgt,0,nclusters);
      wgt_set(cviadjwgt,0,nclusters);
    }

    /* sum moves */
    nmoved = vtx_dlthread_sumreduce(nmoved,comm);
    if (refstats) {
      nunmoved = vtx_dlthread_sumreduce(nunmoved,comm);
      if (myid==0) {
        refstat.nmoves += nmoved;
        refstat.nunmoves += nunmoved;
      }
    }

    if (myid == 0) {
      dprintf("Made "PF_VTX_T" moves in pass "PF_SIZE_T"\n",nmoved,pass);
    }

    #ifdef DEBUG
    DL_ASSERT(check_ucinfo(mgraph->graph,
          (ucinfo_t const * const *)mgraph->ucinfo,
          (cid_t const * const *)clustering->where,nclusters) == 1, 
        "Bad ucinfo after pass "PF_SIZE_T"\n",pass);
    DL_ASSERT(check_clustering(clustering,mgraph->graph) == 1,
        "Bad clustering after pass "PF_SIZE_T"\n",pass);
    dlthread_barrier(comm);
    #endif

    if (nmoved == 0) {
      ++pass;
      break;
    } 
  }

  if (myid==0) {
    dl_free(__prg_cnvtxs);
    dl_free(__prg_ciadjwgt);
    dl_free(__prg_ceadjwgt);
    dl_free(__prg_cviadjwgt);
  }

  if (refstats) {
    tnbnd = vtx_dlthread_sumreduce(mybnd->size,comm);  
    if (myid==0) {
      refstat.npasses = pass;
      refstat.nbnd = tnbnd;
      refinement_stats_push(refstat,refstats);
    }
  }

  /* free combuffers */
  update_combuffer_free(combuffer);

  move_buffer_free(moves);

  dl_free(cdegs);

  /* free cluster information */
  dl_free(cnvtxs);
  dl_free(cviadjwgt);
  dl_free(ciadjwgt);
  dl_free(ceadjwgt);

  vtx_iset_free(mk);
  rv_pq_free(queue);

  if (empty) {
    __par_remove_empty_clusters(clustering,mgraph,comm);
  }

  if (myid == 0) {
    dprintf("Finished refinement of graph at level "PF_SIZE_T"\n",
        mgraph->level);
  }

  return clustering;
}


size_t __prr_done;
vtx_t ** __prr_cnvtxs;
wgt_t ** __prr_ceadjwgt;
wgt_t ** __prr_ciadjwgt;
wgt_t ** __prr_cviadjwgt;
/**
 * @brief Refine a clustering on a graph using random refinement to maximize
 * the modularity of the clustering -- to be called by each thread in a
 * parallel region
 *
 * @param objective The objective containing refinement parameters
 * @param clustering The clustering to refine
 * @param mgraph The clustered graph
 *
 * @return The updated clustering 
 */
static clustering_t * __par_refine_RANDOM(
    clustering_t * const clustering, 
    mgraph_t const * const mgraph,
    size_t const npasses,
    const int empty,
    unsigned int seed,
    refinement_stats_t * const refstats,
    dlthread_comm_t const comm)
{
  int dir;
  size_t pass;
  vtx_t v,i,k,nmoved,nunmoved,tnbnd, sq;
  adj_t j;
  cid_t from, to,c,tmp_to,e, ec, esize;
  wgt_t maxgain, gain,vd,da,db,aconn,bconn,ewgt,oiw,iwgt;
  tid_t o, t;
  move_t move;
  vtx_t * q;
  update_t * uptr;
  update_t up;
  cid_t * nbradj;
  wgt_t * nbrwgt;
  nbrinfo_t * myinfo;
  update_combuffer_t * combuffer;
  update_buffer_t * updates;
  move_buffer_t * moves;
  refinement_stat_t refstat;

  tid_t const myid = dlthread_get_id(comm);

  if (myid == 0) {
    dprintf("Using Modularity Refinement\n");
  }

  /* expose graph stuff */
  graph_t const * const graph = mgraph->graph;    
  tid_t const nthreads = graph->npar;
  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const * const xadj = graph->xadj[myid];
  vtx_t const * const adjncy = graph->adjncy[myid];
  wgt_t const * const adjwgt = (graph->adjwgt ? graph->adjwgt[myid] : NULL);
  wgt_t const * const iadjwgt = (graph->iadjwgt ? graph->iadjwgt[myid] : NULL);
  wgt_t const * const eadjwgt = (graph->eadjwgt ? graph->eadjwgt[myid] : NULL);
  wgt_t const gadjwgt = graph->gadjwgt;
  wgt_t const invgadjwgt = 1.0/gadjwgt;

  /* expose cluster stuff */
  cid_t const nclusters = clustering->nclusters;
  cid_t * const where = clustering->where[myid];
  cluster_t * clusters = clustering->clusters;

  if (refstats) {
    init_refinement_stat(&refstat);
  }

  /* bail out if there's nothing to refine */
  if (nclusters <= 1) {
    return clustering;
  }

  q = vtx_alloc(mynvtxs);

  /* expose mgraph stuff */
  ucinfo_t * const ucinfo = mgraph->ucinfo[myid];
  vtx_iset_t * const mybnd = ucinfo->bnd;

  if (myid==0) {
    __prr_cnvtxs = r_vtx_alloc(nthreads);
    __prr_ciadjwgt = r_wgt_alloc(nthreads);
    __prr_ceadjwgt = r_wgt_alloc(nthreads);
    __prr_cviadjwgt = r_wgt_alloc(nthreads);
  }
  dlthread_barrier(comm);

  /* private cluster vectors */
  wgt_t * cdegs = wgt_alloc(nclusters);
  vtx_t * cnvtxs = __prr_cnvtxs[myid] = vtx_calloc(nclusters);

  /* use one array for reduction purposes later */
  wgt_t * ceadjwgt = __prr_ceadjwgt[myid] = wgt_calloc(nclusters);
  wgt_t * ciadjwgt = __prr_ciadjwgt[myid] = wgt_calloc(nclusters);
  wgt_t * cviadjwgt = __prr_cviadjwgt[myid] = wgt_calloc(nclusters);

  /* create com buffers */
  combuffer = update_combuffer_create(comm);

  /* create move buffers */
  moves = move_buffer_create(dl_max(mynvtxs/4,16));

  /* out pass loop */
  for (pass=0;pass<npasses;++pass) {
    nmoved = 0;
    nunmoved = 0;
    /* clear mark */
    for (dir=0;dir<2;++dir) {
      sq = 0;
      if (myid==0) {
        __prr_done = 0;
      }
      for (c=0;c<nclusters;++c) {
        cdegs[c] = clusters[c].viadjwgt + clusters[c].eadjwgt 
            + clusters[c].iadjwgt;
      }
      /* Add boundary vertices to the queue */
      for (i=0;i<(vtx_t)mybnd->size;++i) {
        v = vtx_iset_get(i,mybnd);
        q[sq++] = v;
      }
      vtx_pseudo_shuffle_r(q,sq/8,sq,&seed);

      /* process vertices in my queue */
      while (sq > 0) {
        v = q[--sq];
        myinfo = ucinfo->nbrinfo+v;

        /* set all variables relevant to v */
        from = where[v];
        nbradj = ucinfo->adjncy + myinfo->nbrstart;
        nbrwgt = ucinfo->adjwgt + myinfo->nbrstart;

        /* don't empty a cluster in kway clustering */
        if (!empty && clusters[from].nvtxs + cnvtxs[from] <= 1) {
          continue;
        }

        to = from;
        maxgain = 0;

        aconn = myinfo->iw;

        if (iadjwgt) {
          vd = (eadjwgt[v]+iadjwgt[v]);
        } else {
          vd = eadjwgt[v];
        }
        da = (cviadjwgt[from] + ceadjwgt[from] +
            ciadjwgt[from]) + cdegs[from] - vd;
        vd *= invgadjwgt; /* scale m out ahead of time */
        /* decide who v is most connected to */
        for (c=0;c<myinfo->nnbrs;++c) {
          tmp_to = nbradj[c];
          if (__valid_move(from,tmp_to,dir)) {
            bconn = nbrwgt[c];
            db = cviadjwgt[tmp_to] + ceadjwgt[tmp_to] +
                ciadjwgt[tmp_to] + cdegs[tmp_to];
            gain = (bconn - aconn) + (vd*(da-db));
            if (maxgain < gain) { 
              maxgain = gain;
              to = nbradj[c];
            }
          }
        }

        /* if we v is in a good spot, leave it there */
        if (to == from) {
          continue;
        }

        /* record the move */
        ++nmoved;
        move.vid = v;
        move.to = to;
        move.from = from;
        move_buffer_add(move,moves);

        /* save info */
        oiw = myinfo->iw;

        if (iadjwgt) {
          iwgt = iadjwgt[v];
        } else {
          iwgt = 0;
        }

        /* actually add the vertex */
        __move_vertex(v,from,to,where,ucinfo,NULL,
            dl_min(nclusters,xadj[v+1]-xadj[v]),iwgt);

        /* vertex and edges */
        --cnvtxs[from];
        ++cnvtxs[to];

        /* vertex-internal edges only get counted by the owning thread */
        if (iwgt) {
          cviadjwgt[from] -= iwgt;
          cviadjwgt[to] += iwgt;
        }

        /* internal edges get counted twice -- once by each thread */
        ciadjwgt[from] -= oiw;
        ciadjwgt[to] += myinfo->iw;

        /* external edges get counted once -- half by each thread */
        ceadjwgt[from] -= ((eadjwgt[v]-oiw)-oiw)/2.0;
        ceadjwgt[to] += ((eadjwgt[v]-myinfo->iw)-myinfo->iw)/2.0;

        /* update neighbors I own */
        for (j=xadj[v];j<xadj[v+1];++j) {
          k = adjncy[j];
          if (k < mynvtxs) {
            if (adjwgt) {
              ewgt = adjwgt[j];
            } else {
              ewgt = 1.0;
            }
            c = where[k];
            if (iadjwgt) {
              iwgt = iadjwgt[k];
            } else {
              iwgt = 0;
            }
            __par_update_vertex(k,where[k],ewgt,from,to,invgadjwgt,
                ucinfo,ceadjwgt,ciadjwgt,NULL,NULL,
                dl_min(nclusters,xadj[k+1]-xadj[k]),iwgt);
          } else {
            /* tell the other thread about it once we decide to keep it */
          }
        }
      }
      /* clear out any extra vertices */
      sq = 0;

      dlthread_atomic_add(__prr_done,comm);

      dlthread_barrier(comm);
      for (e=myid;e<(nclusters/CLUSTER_CHUNK)+1;e+=nthreads) {
        ec = e*CLUSTER_CHUNK;
        esize = dl_min(CLUSTER_CHUNK,nclusters-ec);
        /* cnvtxs */
        for (t=(myid+1)%nthreads;t!=myid;t=(t+1)%nthreads) {
          for (c=0;c<esize;++c) {
            cnvtxs[c+ec] += __prr_cnvtxs[t][c+ec];
          }
        }
        for (c=0;c<esize;++c) {
          clusters[c+ec].nvtxs += cnvtxs[c+ec];
        }
        /* eadjwgt */
        for (t=(myid+1)%nthreads;t!=myid;t=(t+1)%nthreads) {
          for (c=0;c<esize;++c) {
            ceadjwgt[c+ec] += __prr_ceadjwgt[t][c+ec];
          }
        }
        for (c=0;c<esize;++c) {
          clusters[c+ec].eadjwgt += ceadjwgt[c+ec];
        }
        /* iadjwgt */
        for (t=(myid+1)%nthreads;t!=myid;t=(t+1)%nthreads) {
          for (c=0;c<esize;++c) {
            ciadjwgt[c+ec] += __prr_ciadjwgt[t][c+ec];
          }
        }
        for (c=0;c<esize;++c) {
          clusters[c+ec].iadjwgt += ciadjwgt[c+ec];
        }
        /* viadjwgt */
        for (t=(myid+1)%nthreads;t!=myid;t=(t+1)%nthreads) {
          for (c=0;c<esize;++c) {
            cviadjwgt[c+ec] += __prr_cviadjwgt[t][c+ec];
          }
        }
        for (c=0;c<esize;++c) {
          clusters[c+ec].viadjwgt += cviadjwgt[c+ec];
        }
      }

      dlthread_barrier(comm);

      vtx_set(cnvtxs,0,nclusters);
      wgt_set(ceadjwgt,0,nclusters);
      wgt_set(ciadjwgt,0,nclusters);
      wgt_set(cviadjwgt,0,nclusters);

      for (c=0;c<nclusters;++c) {
        cdegs[c] = clusters[c].viadjwgt + clusters[c].eadjwgt 
            + clusters[c].iadjwgt;
      }

      /* undo bad moves */
      for (i=moves->size;i>0;) {
        move = moves->elements[--i];
        v = move.vid;
        to = move.to;
        from = move.from;

        myinfo = ucinfo->nbrinfo+v;
        nbradj = ucinfo->adjncy + myinfo->nbrstart;
        nbrwgt = ucinfo->adjwgt + myinfo->nbrstart;

        /* calculate the gain of undoing the move */
        aconn = myinfo->iw;
        if (iadjwgt) {
          vd = (eadjwgt[v]+iadjwgt[v]);
        } else {
          vd = eadjwgt[v];
        }
        da = (cviadjwgt[to] + ceadjwgt[to] +
            ciadjwgt[to]) + cdegs[to] - vd;
        vd /= gadjwgt; /* scale m out ahead of time */
        for (c=0;c<myinfo->nnbrs;++c) {
          if (nbradj[c] == from) {
            break;
          }
        }
        if (c < myinfo->nnbrs) {
          bconn = nbrwgt[c];
        } else {
          bconn = 0;
        }
        db = cviadjwgt[from] + ceadjwgt[from] +
            ciadjwgt[from] + cdegs[from];
        gain = (bconn - aconn) + (vd*(da-db));

        if (gain > 0 || (!empty && clusters[from].nvtxs + cnvtxs[from] == 0)) {
          /* undo move */
          --nmoved;
          ++nunmoved;

          /* save info */
          oiw = myinfo->iw;

          if (iadjwgt) {
            iwgt = iadjwgt[v];
          } else {
            iwgt = 0;
          }

          /* actually add the vertex */
          __move_vertex(v,to,from,where,ucinfo,NULL,
              dl_min(nclusters,xadj[v+1]-xadj[v]),iwgt);

          /* vertex and edges */
          --cnvtxs[to];
          ++cnvtxs[from];

          /* vertex-internal edges only get counted by the owning thread */
          if (iwgt) {
            cviadjwgt[to] -= iwgt;
            cviadjwgt[from] += iwgt;
          }

          /* internal edges get counted twice -- once by each thread */
          ciadjwgt[to] -= oiw;
          ciadjwgt[from] += myinfo->iw;

          /* external edges get counted once -- half by each thread */
          ceadjwgt[to] -= (((eadjwgt[v]-oiw)-oiw)/2.0);
          ceadjwgt[from] += (((eadjwgt[v]-myinfo->iw)-myinfo->iw)/2.0);

          /* update neighbors I own */
          for (j=xadj[v];j<xadj[v+1];++j) {
            k = adjncy[j];
            if (k < mynvtxs) {
              if (adjwgt) {
                ewgt = adjwgt[j];
              } else {
                ewgt = 1.0;
              }

              if (iadjwgt) {
                iwgt = iadjwgt[k];
              } else {
                iwgt = 0;
              }

              __par_update_vertex(k,where[k],ewgt,to,from,
                  invgadjwgt,ucinfo,ceadjwgt,ciadjwgt,NULL,NULL,
                  dl_min(nclusters,xadj[k+1]-xadj[k]),iwgt);
            } else {
              /* The other thread never needs to know */
            }
          }
        } else {
          /* tell other threads since we're keeping the move */
          for (j=xadj[v];j<xadj[v+1];++j) {
            k = adjncy[j];
            if (k >= mynvtxs) {
              /* tell the other thread about it later */
              o = gvtx_to_tid(k,graph->dist);
              k = gvtx_to_lvtx(k,graph->dist);
              up.vid = k;
              if (adjwgt) {
                up.wgt = adjwgt[j];
              } else {
                up.wgt = 1.0f;
              }
              up.to = to;
              up.from = from;
              update_combuffer_add(o,up,combuffer);
            }
          }
        }
      }
      move_buffer_clear(moves);

      /* send updates */
      update_combuffer_send(combuffer);
        
      /* process updates */
      for (o=(myid+1)%nthreads;o!=myid;o=(o+1)%nthreads) {
        updates = update_combuffer_get(o,combuffer);
        uptr = updates->elements;
        for (i=0;i<(vtx_t)updates->size;++i) {
          k = uptr[i].vid;
          ewgt = uptr[i].wgt;
          to = uptr[i].to;
          from = uptr[i].from;

          if (iadjwgt) {
            iwgt = iadjwgt[k];
          } else {
            iwgt = 0;
          }

          __par_update_vertex(k,where[k],ewgt,from,to,invgadjwgt,ucinfo,
              ceadjwgt,ciadjwgt,NULL,NULL,
              dl_min(nclusters,xadj[k+1]-xadj[k]),iwgt);
        }
      }
      update_combuffer_clear(combuffer);

      dlthread_barrier(comm);
      for (e=myid;e<(nclusters/CLUSTER_CHUNK)+1;e+=nthreads) {
        ec = e*CLUSTER_CHUNK;
        esize = dl_min(CLUSTER_CHUNK,nclusters-ec);
        /* cnvtxs */
        for (t=(myid+1)%nthreads;t!=myid;t=(t+1)%nthreads) {
          for (c=0;c<esize;++c) {
            cnvtxs[c+ec] += __prr_cnvtxs[t][c+ec];
          }
        }
        for (c=0;c<esize;++c) {
          clusters[c+ec].nvtxs += cnvtxs[c+ec];
        }
        /* eadjwgt */
        for (t=(myid+1)%nthreads;t!=myid;t=(t+1)%nthreads) {
          for (c=0;c<esize;++c) {
            ceadjwgt[c+ec] += __prr_ceadjwgt[t][c+ec];
          }
        }
        for (c=0;c<esize;++c) {
          clusters[c+ec].eadjwgt += ceadjwgt[c+ec];
        }
        /* iadjwgt */
        for (t=(myid+1)%nthreads;t!=myid;t=(t+1)%nthreads) {
          for (c=0;c<esize;++c) {
            ciadjwgt[c+ec] += __prr_ciadjwgt[t][c+ec];
          }
        }
        for (c=0;c<esize;++c) {
          clusters[c+ec].iadjwgt += ciadjwgt[c+ec];
        }
        /* viadjwgt */
        for (t=(myid+1)%nthreads;t!=myid;t=(t+1)%nthreads) {
          for (c=0;c<esize;++c) {
            cviadjwgt[c+ec] += __prr_cviadjwgt[t][c+ec];
          }
        }
        for (c=0;c<esize;++c) {
          clusters[c+ec].viadjwgt += cviadjwgt[c+ec];
        }
      }

      dlthread_barrier(comm);

      vtx_set(cnvtxs,0,nclusters);
      wgt_set(ceadjwgt,0,nclusters);
      wgt_set(ciadjwgt,0,nclusters);
      wgt_set(cviadjwgt,0,nclusters);
    }

    /* sum moves */
    nmoved = vtx_dlthread_sumreduce(nmoved,comm);
    if (refstats) {
      nunmoved = vtx_dlthread_sumreduce(nunmoved,comm);
      if (myid==0) {
        refstat.nmoves += nmoved;
        refstat.nunmoves += nunmoved;
      }
    }

    if (myid == 0) {
      dprintf("Made "PF_VTX_T" moves in pass "PF_SIZE_T"\n",nmoved,pass);
    }

    #ifdef DEBUG
    DL_ASSERT(check_ucinfo(mgraph->graph,
          (ucinfo_t const * const *)mgraph->ucinfo,
          (cid_t const * const *)clustering->where,nclusters) == 1, 
        "Bad ucinfo after pass "PF_SIZE_T"\n",pass);
    DL_ASSERT(check_clustering(clustering,mgraph->graph) == 1,
        "Bad clustering after pass "PF_SIZE_T"\n",pass);
    dlthread_barrier(comm);
    #endif

    if (nmoved == 0) {
      ++pass;
      break;
    } 
  }

  if (myid == 0) {
    dl_free(__prr_cnvtxs);
    dl_free(__prr_ciadjwgt);
    dl_free(__prr_ceadjwgt);
    dl_free(__prr_cviadjwgt);
  }

  if (refstats) {
    tnbnd = vtx_dlthread_sumreduce(mybnd->size,comm);  
    if (myid==0) {
      refstat.npasses = pass;
      refstat.nbnd = tnbnd;
      refinement_stats_push(refstat,refstats);
    }
  }

  /* free combuffers */
  update_combuffer_free(combuffer);

  move_buffer_free(moves);

  dl_free(cdegs);

  /* free cluster information */
  dl_free(cnvtxs);
  dl_free(cviadjwgt);
  dl_free(ciadjwgt);
  dl_free(ceadjwgt);
  dl_free(q);


  if (empty) {
    __par_remove_empty_clusters(clustering,mgraph,comm);
  }

  if (myid == 0) {
    dprintf("Finished refinement of graph at level "PF_SIZE_T"\n",
        mgraph->level);
  }

  return clustering;
}



/******************************************************************************
* PUBLIC SERIAL FUNCTIONS *****************************************************
******************************************************************************/


clustering_t * refine_graph(
    objective_t * const objective, 
    clustering_t * const clustering, 
    mgraph_t const * const mgraph)
{
  clustering_t * fc;

  const int empty = objective->parttype != NERSTRAND_PARTITION_KWAY;

  dl_start_timer(&(objective->timers.refinement));
  fc = explicit_refine_graph(clustering,mgraph,objective->reftype,
      objective->nrefpass,empty,objective->seed++,objective->refstats);
  dl_stop_timer(&(objective->timers.refinement));

  return fc;
}


clustering_t * explicit_refine_graph(
    clustering_t * const clustering, 
    mgraph_t const * const mgraph,
    const nerstrand_reftype_t type,
    size_t const npasses,
    const int empty,
    unsigned int seed,
    refinement_stats_t * refstats)
{
  clustering_t * fc;
  dprintf("Performing serial refinement on level "PF_SIZE_T" with "PF_VTX_T
      " vertices\n",mgraph->level,mgraph->graph->nvtxs);

  switch(type) {
    case NERSTRAND_REFINEMENT_GREEDY :
      fc = __refine_GREEDY(clustering,mgraph,npasses,empty,refstats);
      break;
    case NERSTRAND_REFINEMENT_RANDOM:
      fc = __refine_RANDOM(clustering,mgraph,npasses,empty,seed,refstats);
      break;
    default :
      dl_error("Unknown refinement type\n");
      break;
  }

  DL_ASSERT(check_clustering(clustering,mgraph->graph) == 1,
      "Bad clustering after refinement\n");

  return fc;
}




/******************************************************************************
* PUBLIC PARALLEL FUNCTIONS ***************************************************
******************************************************************************/


clustering_t * par_refine_graph(
    objective_t * const objective, 
    clustering_t * const clustering, 
    mgraph_t const * const mgraph)
{
  clustering_t * fc;

  tid_t const myid = dlthread_get_id(objective->comm);
  int const empty = objective->parttype != NERSTRAND_PARTITION_KWAY;

  if (myid == 0) {
    dl_start_timer(&(objective->timers.refinement));
  }

  fc = par_explicit_refine_graph(clustering,mgraph,objective->reftype,
      objective->nrefpass,empty,objective->seed+myid,
      objective->refstats,objective->comm);

  if (myid == 0) {
    dl_stop_timer(&(objective->timers.refinement));
  }

  return fc;
}


clustering_t * par_explicit_refine_graph(
    clustering_t * const clustering, 
    mgraph_t const * const mgraph,
    const nerstrand_reftype_t type,
    size_t const npasses,
    const int empty,
    unsigned int seed,
    refinement_stats_t * refstats,
    dlthread_comm_t const comm)
{
  clustering_t * fc;

  const tid_t myid = dlthread_get_id(comm);

  dprintf("Performing parallel refinement on level "PF_SIZE_T" with "PF_VTX_T
      " vertices\n",mgraph->level,mgraph->graph->nvtxs);

  DL_ASSERT(check_clustering(clustering,mgraph->graph) == 1,
      "Bad clustering before refinement\n");

  switch(type) {
    case NERSTRAND_REFINEMENT_GREEDY :
      fc = __par_refine_GREEDY(clustering,mgraph,npasses,empty,refstats, \
          comm);
      break;
    case NERSTRAND_REFINEMENT_RANDOM:
      fc = __par_refine_RANDOM(clustering,mgraph,npasses,empty,seed+myid, \
          refstats,comm);
      break;
    default :
      dl_error("Unknown refinement type\n");
      break;
  }

  DL_ASSERT(check_clustering(clustering,mgraph->graph) == 1,
      "Bad clustering after refinement\n");

  return fc;
}




#endif
