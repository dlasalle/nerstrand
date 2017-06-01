/**
 * @file cluster.c
 * @brief Functions for cluster operations
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-10-08
 */




#ifndef NERSTRAND_CLUSTER_C
#define NERSTRAND_CLUSTER_C




#include "cluster.h"




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLMEM_PREFIX cluster
#define DLMEM_TYPE_T cluster_t
#define DLMEM_DLTYPE DLTYPE_STRUCT
#define DLMEM_INITFUNCTION init_cluster
#include "dlmem_funcs.h"
#undef DLMEM_INITFUNCTION
#undef DLMEM_DLTYPE
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX


#define DLMEM_PREFIX clustering
#define DLMEM_TYPE_T clustering_t
#define DLMEM_DLTYPE DLTYPE_STRUCT
#define DLMEM_INITFUNCTION init_clustering
#include "dlmem_funcs.h"
#undef DLMEM_INITFUNCTION
#undef DLMEM_DLTYPE
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX




/******************************************************************************
* PUBLIC SERIAL FUNCTIONS *****************************************************
******************************************************************************/


clustering_t * init_clustering(
    clustering_t * const clustering)
{
  /* zero */
  clustering->nclusters = 0;
  clustering->clusters = NULL;
  clustering->where = NULL;
  clustering->free_clusters = 0;
  clustering->free_where = 0;
  return clustering;
}


cluster_t * init_cluster(
    cluster_t * const cluster)
{
  /* zero */
  cluster->id = 0;
  cluster->nvtxs = 0;
  cluster->eadjwgt = 0;
  cluster->iadjwgt = 0;
  cluster->viadjwgt = 0;
  return cluster;
}


clustering_t * setup_clustering(
    cid_t const nclusters,  
    cid_t ** const where, 
    graph_t const * const graph)
{
  vtx_t i,k,mynvtxs;
  adj_t j;
  cid_t c,me;
  tid_t o,myid;
  wgt_t ewgt, tew;

  tid_t const nthreads = graph->npar;

  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy;
  wgt_t const * const * const adjwgt = (wgt_t const * const *)graph->adjwgt;
  wgt_t const * const * const iadjwgt = (wgt_t const * const *)graph->iadjwgt;
  wgt_t const * const * const eadjwgt = (wgt_t const * const *)graph->eadjwgt;

  clustering_t * clustering = clustering_calloc(1);
  cluster_t * cluster;

  clustering->nclusters = nclusters;
  clustering->npar = nthreads;
  clustering->where = where;

  clustering->clusters = cluster_calloc(nclusters);
  for (c=0;c<nclusters;++c) {
    clustering->clusters[c].id = c;
  }
  /* count vertices per cluster */
  if (graph->nvtxs > nclusters) {
    for (myid=0;myid<nthreads;++myid) {
      mynvtxs = graph->mynvtxs[myid];
      for (i=0;i<mynvtxs;++i) {
        me = where[myid][i];
        DL_ASSERT(me < nclusters,
            "Bad cluster number "PF_CID_T"/"PF_CID_T" for vertex "PF_VTX_T":"
            PF_TID_T"\n",me,nclusters,i,myid);
        cluster = clustering->clusters+me;
        if (iadjwgt) {
          cluster->viadjwgt += iadjwgt[myid][i];
        }
        ++cluster->nvtxs;
        tew = eadjwgt[myid][i];
        for (j=xadj[myid][i];j<xadj[myid][i+1];++j) {
          k = adjncy[myid][j]; 
          if (k < mynvtxs) {
            o = myid;
          } else {
            o = gvtx_to_tid(k,graph->dist);
            k = gvtx_to_lvtx(k,graph->dist);
          }
          if (adjwgt) {
            ewgt = adjwgt[myid][j];
          } else {
            ewgt = 1.0;
          }
          if (me == where[o][k]) {
            cluster->iadjwgt += ewgt;
          } else {
            cluster->eadjwgt += ewgt;
          }
          tew -= ewgt;
        }
        cluster->eadjwgt += tew;
      }
    }
  } else {
    for (myid=0;myid<nthreads;++myid) {
      mynvtxs = graph->mynvtxs[myid];
      for (i=0;i<mynvtxs;++i) {
        me = where[myid][i];
        cluster = clustering->clusters+me;
        if (iadjwgt) {
          cluster->viadjwgt = iadjwgt[myid][i];
        } else {
          cluster->viadjwgt = 0.0;
        }
        cluster->nvtxs = 1;
        cluster->iadjwgt = 0.0;
        cluster->eadjwgt = eadjwgt[myid][i];
      }
    }
  }

  clustering->free_clusters = 1;
  clustering->free_where = 1;

  return clustering;
}


clustering_t * clone_clustering(
    vtx_t const * const nvtxs, 
    clustering_t const * const clustering)
{
  tid_t t;

  clustering_t * nclustering = clustering_alloc(1);

  nclustering->nclusters = clustering->nclusters;
  nclustering->npar = clustering->npar;
  nclustering->where = r_cid_alloc(clustering->npar);
  for (t=0;t<clustering->npar;++t) {
    nclustering->where[t] = cid_duplicate(clustering->where[t],nvtxs[t]);
  }
  nclustering->clusters = cluster_duplicate(clustering->clusters,
      clustering->nclusters);

  nclustering->free_clusters = 1;
  nclustering->free_where = 1;

  return nclustering;
}


cluster_t * setup_cluster(
    cid_t const id, 
    vtx_t const nvtxs)
{
  cluster_t * cluster = cluster_calloc(1);

  cluster->id = id;
  cluster->nvtxs = nvtxs;

  return cluster;
}


void unify_where(
    cid_t const * const * const where, 
    graph_t const * const graph, 
    cid_t * const uwhere)
{
  vtx_t i,mynvtxs,k;
  tid_t myid;

  tid_t const nthreads = graph->npar;

  if (graph->alias) {
    for (myid=0;myid<nthreads;++myid) {
      mynvtxs = graph->mynvtxs[myid];
      for (i=0;i<mynvtxs;++i) {
        uwhere[graph->alias[myid][i]] = where[myid][i];
      }
    }
  } else {
    k = 0;
    for (myid=0;myid<nthreads;++myid) {
      mynvtxs = graph->mynvtxs[myid];
      for (i=0;i<mynvtxs;++i) {
        uwhere[k++] = where[myid][i];
      }
    }
  }
}


int calc_clustering(
    clustering_t * const clustering, 
    graph_t const * const graph)
{
  vtx_t i,k;
  adj_t j;
  tid_t t,o;
  cid_t c;

  DL_ASSERT_EQUALS(check_clustering(clustering,graph),1,"%d");

  tid_t const nthreads = clustering->npar;

  DL_ASSERT_EQUALS(nthreads,graph->npar,PF_TID_T);

  cid_t const nclusters = clustering->nclusters;

  cluster_t * const clusters = clustering->clusters;

  cid_t const * const * const where = (cid_t const*const*)clustering->where;

  /* initialize our clusters */
  for (c=0;c<nclusters;++c) {
    init_cluster(clusters+c);
    clustering->clusters[c].id = c;
    DL_ASSERT_EQUALS((vtx_t)0,clusters[c].nvtxs,PF_VTX_T);
  }

  /* find out private values for each cluster */
  for (t=0;t<nthreads;++t) {
    for (i=0;i<graph->mynvtxs[t];++i) {
      c = where[t][i];
      ++clusters[c].nvtxs;
      if (graph->iadjwgt) {
        clusters[c].viadjwgt += graph->iadjwgt[t][i];
      }
      for (j=graph->xadj[t][i];j<graph->xadj[t][i+1];++j) {
        k = graph->adjncy[t][j];
        if (k < graph->mynvtxs[t]) {
          o = t;
        } else {
          o = gvtx_to_tid(k,graph->dist);
          k = gvtx_to_lvtx(k,graph->dist);
        }
        if (c == where[o][k]) {
          if (graph->adjwgt) {
            clusters[c].iadjwgt += graph->adjwgt[t][j];
          } else {
            clusters[c].iadjwgt += 1.0;
          }
        } else {
          if (graph->adjwgt) {
            clusters[c].eadjwgt += graph->adjwgt[t][j];
          } else {
            clusters[c].eadjwgt += 1.0;
          }
        }
      }
    }
  }

  DL_ASSERT_EQUALS(check_clustering(clustering,graph),1,"%d");

  return 1;
}


int free_clustering(
    clustering_t * clustering)
{
  tid_t i;
  if (clustering->free_clusters) {
    dl_free(clustering->clusters);
  }
  if (clustering->free_where) {
    for (i=0;i<clustering->npar;++i) {
      dl_free(clustering->where[i]);
    }
    dl_free(clustering->where);
  }
  dl_free(clustering);

  return 1;
}


int free_cluster(
    cluster_t * cluster)
{
  dl_free(cluster);
  return 1;
}


mod_t calc_modularity(
    cluster_t const * const cluster, 
    graph_t const * const graph)
{
  mod_t q;

  twgt_t m = graph->gadjwgt;
  twgt_t id = cluster->iadjwgt + cluster->viadjwgt;
  twgt_t ed = cluster->eadjwgt;
  twgt_t td = id + ed;
  /*
   * Q = 1 / 2*edges (internal_edges - expected_edges) + voodoo math
   */
  q = (mod_t)((1.0 / m) * (id - ((td*td)/m))); 

  return q;
}


mod_t calc_total_modularity(
    clustering_t const * const clustering,
    graph_t const * const graph)
{
  cid_t i;
  long double q = 0;
  for (i=0;i<clustering->nclusters;++i) {
    q += calc_modularity(clustering->clusters+i,graph);
  }
  return (mod_t)q;
}


int check_clustering(
    clustering_t const * clustering, 
    graph_t const * graph)
{
  int rv = 1;
  vtx_t i,k,gnvtxs;
  adj_t j;
  cid_t c;
  tid_t myid, o;
  vtx_t mynvtxs;
  wgt_t tew;
  long double gadjwgt;

  tid_t const nthreads = graph->npar;
  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy;
  wgt_t const * const * const adjwgt = (wgt_t const * const *)graph->adjwgt;
  wgt_t const * const * const iadjwgt = (wgt_t const * const *)graph->iadjwgt;
  wgt_t const * const * const eadjwgt = (wgt_t const * const *)graph->eadjwgt;

  cid_t const * const * const where = (cid_t const * const *)clustering->where;
  cid_t const nclusters = clustering->nclusters;
  cluster_t const * const clusters = clustering->clusters;

  vtx_t * cnvtxs = vtx_calloc(nclusters);
  twgt_t * cwgts = twgt_calloc(nclusters);
  twgt_t * ewgts = twgt_calloc(nclusters);
  twgt_t * iwgts = twgt_calloc(nclusters);
  twgt_t * viwgts = twgt_calloc(nclusters);

  gnvtxs = 0;

  /* compute real cluster statistics */
  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    for (i=0;i<mynvtxs;++i) {
      c = where[myid][i];
      if (c >= nclusters) {
        eprintf("The vertex "PF_VTX_T" is in an invalid cluster "PF_CID_T
            "/"PF_CID_T"\n",i,c,nclusters);
        return 0;
      }
      ++cnvtxs[c];
      if (iadjwgt) {
        viwgts[c] += iadjwgt[myid][i];
      }
      tew = eadjwgt[myid][i];
      for (j=xadj[myid][i];j<xadj[myid][i+1];++j) {
        k = adjncy[myid][j];
        ++gnvtxs;
        if (k < mynvtxs) {
          o = myid;
        } else {
          o = gvtx_to_tid(k,graph->dist);
          k = gvtx_to_lvtx(k,graph->dist);
        }
        if (c == where[o][k]) {
          if (adjwgt) {
            iwgts[c] += adjwgt[myid][j];  
            tew -= adjwgt[myid][j];
          } else {
            iwgts[c] += 1.0;
            tew -= 1.0;
          } 
        }
      }
      ewgts[c] += tew;
    }
  }

  /* verify computed cluster statistics */
  for (c=0;c<nclusters;++c) {
    if (cnvtxs[c] != clusters[c].nvtxs) {
      eprintf("Cluster "PF_CID_T" thinks it has "PF_VTX_T" vertices but has "
          PF_VTX_T" vertices\n",c,clusters[c].nvtxs,cnvtxs[c]);
      rv = 0;
    }
    if (!dl_near_equal(viwgts[c],clusters[c].viadjwgt)) {
      eprintf("Cluster "PF_CID_T" thinks it has an vertex-internal edge weight" 
          " of "PF_TWGT_T" but has a weight of "PF_TWGT_T"\n",c,
          clusters[c].viadjwgt,viwgts[c]);
      rv = 0;
    }
    if (!dl_near_equal(iwgts[c],clusters[c].iadjwgt)) {
      eprintf("Cluster "PF_CID_T" thinks it has an internal edge weight of "
          PF_TWGT_T" but has a weight of "PF_TWGT_T"\n",c,clusters[c].iadjwgt,
          iwgts[c]);
      rv = 0;
    }
    if (!dl_near_equal(ewgts[c],clusters[c].eadjwgt)) {
      eprintf("Cluster "PF_CID_T" thinks it has an external edge weight of "
          PF_TWGT_T" but has a weight of "PF_TWGT_T"\n",c,clusters[c].eadjwgt,
          ewgts[c]);
      rv = 0;
    }
  }

  /* verify totals */
  gadjwgt = 0.0;
  for (c=0;c<nclusters;++c) {
    gadjwgt += clusters[c].viadjwgt + clusters[c].iadjwgt + 
        clusters[c].eadjwgt;
  }
  if (!dl_near_equal(gadjwgt,graph->gadjwgt)) {
    eprintf("Total cluster edge weights (%Lf) do not equal total "
      "graph edge weight ("PF_TWGT_T")\n",gadjwgt,graph->gadjwgt);
    rv = 0;
  }

  dl_free(cnvtxs);
  dl_free(cwgts); 
  dl_free(ewgts); 
  dl_free(iwgts); 
  dl_free(viwgts); 

  return rv;
}




/******************************************************************************
* PUBLIC PARALLEL FUNCTIONS ***************************************************
******************************************************************************/


clustering_t * par_setup_clustering(
    cid_t const nclusters,  
    vtx_t const nvtxs,
    const int clusterless_vertices, 
    cluster_t * const clusters,
    cid_t ** const where,
    dlthread_comm_t const comm)
{
  vtx_t i, tmpnvtxs;
  cid_t c, cblock, cstart, cend;
  vtx_t * cnvtxs;
  cluster_t * cluster;

  static clustering_t * clustering;

  tid_t const myid = dlthread_get_id(comm);
  tid_t const nthreads = dlthread_get_nthreads(comm);

  if (myid == 0) {
    clustering = clustering_calloc(1);
    clustering->nclusters = nclusters;
    if (where) {
      clustering->where = where;
    } else {
      clustering->where = r_cid_calloc(nthreads);
    }
    clustering->npar = nthreads;
    clustering->free_clusters = 1;
    clustering->free_where = 1;
  }
  dlthread_barrier(comm);

  if (clusters) {
    if (myid == 0) {
      clustering->clusters = clusters;
    }
  } else {
    if (myid == 0) {
      clustering->clusters = cluster_calloc(nclusters);
    }
    dlthread_barrier(comm);

    if (where) {
      cnvtxs = vtx_calloc(nclusters);
      /* count vertices per cluster */
      for (i=0;i<nvtxs;++i) {
        ++cnvtxs[where[myid][i]];
      }
      for (c=0;c<nclusters;++c) {
        cluster = clustering->clusters+c;
        tmpnvtxs = vtx_dlthread_sumreduce(cnvtxs[c],comm);
        if (myid == 0) {
          cluster->nvtxs = tmpnvtxs;
        }
      }
    } else {
      for (cblock=myid;cblock*64<nclusters;cblock+=nthreads) {
        cstart = cblock*64;
        cend = dl_min((cblock+1)*64,nclusters);
        for (c=cstart;c<cend;++c) {
          clustering->clusters[c].nvtxs = 0;
        }
      }
    }
  }

  return clustering;
}


void par_calc_clustering(
    clustering_t * const clustering, 
    graph_t const * const graph,
    dlthread_comm_t const comm)
{
  vtx_t i,k;
  adj_t j;
  cid_t c, d, cblock, cstart, cend;
  tid_t o;
  cluster_t * cluster;

  tid_t const myid = dlthread_get_id(comm);
  tid_t const nthreads = dlthread_get_nthreads(comm);

  DL_ASSERT_EQUALS(nthreads,clustering->npar,PF_TID_T);
  DL_ASSERT_EQUALS(nthreads,graph->npar,PF_TID_T);

  cid_t const nclusters = clustering->nclusters;

  cid_t const size = dl_max(nthreads,nclusters);
  cid_t const * const where = clustering->where[myid];

  /* initialize our clusters */
  for (cblock=myid;cblock*64<nclusters;cblock+=nthreads) {
    cstart = cblock*64;
    cend = dl_min((cblock+1)*64,nclusters);
    for (c=cstart;c<cend;++c) {
      init_cluster(clustering->clusters+c);
      clustering->clusters[c].id = c;
    }
  }

  /* create a private cluster arrays */
  adj_t * nedges = adj_calloc(clustering->nclusters);
  vtx_t * cnvtxs = vtx_calloc(clustering->nclusters);
  wgt_t * viadjwgt = wgt_calloc(clustering->nclusters);
  wgt_t * iadjwgt = wgt_calloc(clustering->nclusters);
  wgt_t * eadjwgt = wgt_calloc(clustering->nclusters);

  /* find out private values for each cluster */
  for (i=0;i<graph->mynvtxs[myid];++i) {
    c = where[i];
    ++cnvtxs[c];
    if (graph->iadjwgt) {
      viadjwgt[c] += graph->iadjwgt[myid][i];
    }
    for (j=graph->xadj[myid][i];j<graph->xadj[myid][i+1];++j) {
      k = graph->adjncy[myid][j];
      o = gvtx_to_tid(k,graph->dist);
      if (c == clustering->where[o][k]) {
        ++nedges[c];
        if (graph->adjwgt) {
          iadjwgt[c] += graph->adjwgt[myid][j];
        } else {
          iadjwgt[c] += 1.0;
        }
      } else {
        if (graph->adjwgt) {
          eadjwgt[c] += graph->adjwgt[myid][j];
        } else {
          eadjwgt[c] += 1.0;
        }
      }
    }
  }

  /* sum and store private values */
  for (c=0;c<size;++c) {
    d = (c + myid) % size;
    if (d < nclusters) {
      cluster = clustering->clusters+c;
      cluster->nvtxs += cnvtxs[c];
      cluster->viadjwgt += viadjwgt[c];
      cluster->iadjwgt += iadjwgt[c];
      cluster->eadjwgt += eadjwgt[c];
    }
    dlthread_barrier(comm);
  }

  dl_free(cnvtxs);
  dl_free(viadjwgt);
  dl_free(iadjwgt);
  dl_free(eadjwgt);
}


mod_t par_calc_total_modularity(
    clustering_t const * const clustering,
    graph_t const * const graph,
    dlthread_comm_t const comm)
{
  cid_t c, cblock, cstart, cend;
  static long double q;

  cid_t const nclusters = clustering->nclusters;

  tid_t const myid = dlthread_get_id(comm);
  tid_t const nthreads = dlthread_get_nthreads(comm);

  for (cblock=myid;cblock*64<nclusters;cblock+=nthreads) {
    cstart = cblock*64;
    cend = dl_min((cblock+1)*64,nclusters);
    for (c=cstart;c<cend;++c) {
      q += calc_modularity(clustering->clusters+c,graph);
    }
  }

  ld_dlthread_sumreduce(q,comm);

  return (mod_t)q;
}


void par_free_clustering(
    clustering_t * clustering,
    dlthread_comm_t const comm)
{
  tid_t const myid = dlthread_get_id(comm);

  dlthread_barrier(comm);
  if (myid == 0) {
    if (clustering->free_clusters) {
      dl_free(clustering->clusters);
    }
  }
  if (clustering->free_where) {
    dl_free(clustering->where[myid]);
    dlthread_barrier(comm);

    if (myid == 0) {
      dl_free(clustering->where);
    }
  }
  dlthread_barrier(comm);

  if (myid == 0) {
    dl_free(clustering);
  }
}




#endif
