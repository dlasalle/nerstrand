/**
 * @file uncoarsen.c
 * @brief Functions for uncoarsening a graph
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2014-03-03
 */




#ifndef NERSTRAND_UNCOARSEN_C
#define NERSTRAND_UNCOARSEN_C




#include "uncoarsen.h"
#include "project.h"
#include "refine.h"




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


static const adj_t MIN_NADJ = 0x4000;




/******************************************************************************
* PUBLIC SERIAL FUNCTIONS *****************************************************
******************************************************************************/


clustering_t * uncoarsen_graph(
    objective_t * objective, 
    clustering_t * clustering, 
    mgraph_t * const fmgraph, 
    const mgraph_t * const cmgraph)
{
  dl_start_timer(&(objective->timers.uncoarsening));

  DL_ASSERT(check_graph(cmgraph->graph) == 1, "Bad graph passed to "
      "uncoarsen_graph\n");
  DL_ASSERT(check_clustering(clustering,cmgraph->graph) == 1,
      "Bad clustering at start of uncoarsen_graph\n");
  DL_ASSERT(check_ucinfo(cmgraph->graph,
        (const ucinfo_t * const *)cmgraph->ucinfo,
        (const cid_t * const *)clustering->where,clustering->nclusters) == 1, 
      "Bad ucinfo passed to uncoarsen_graph\n");
  DL_ASSERT(check_bnd(cmgraph,(const cid_t * const *)clustering->where) == 1,
      "Bad boundary at start of uncoarsen_graph\n");

  clustering = project_clustering(objective,clustering,fmgraph,cmgraph);

  DL_ASSERT(check_clustering(clustering,fmgraph->graph) == 1,
      "Bad clustering after projection\n");
  DL_ASSERT(check_ucinfo(fmgraph->graph,
        (const ucinfo_t * const *)fmgraph->ucinfo,
        (const cid_t * const *)clustering->where,clustering->nclusters) == 1, 
      "Bad ucinfo after projection\n");
  DL_ASSERT(check_bnd(fmgraph,(const cid_t * const *)clustering->where) == 1,
      "Bad boundary after projection\n");

  clustering = refine_graph(objective,clustering,fmgraph);

  DL_ASSERT(check_clustering(clustering,fmgraph->graph) == 1,
      "Bad clustering after refinement\n");
  DL_ASSERT(check_ucinfo(fmgraph->graph,
        (const ucinfo_t * const *)fmgraph->ucinfo,
        (const cid_t * const *)clustering->where,clustering->nclusters) == 1, 
      "Bad ucinfo after refinement\n");
  DL_ASSERT(check_bnd(fmgraph,(const cid_t * const *)clustering->where) == 1,
      "Bad boundary after refinement\n");

  dl_stop_timer(&(objective->timers.uncoarsening));

  return clustering;
}


ucinfo_t ** derive_ucinfo(
    const clustering_t * const clustering, 
    const mgraph_t * const mgraph)
{
  vtx_t i,k,l,mynvtxs,tin;
  adj_t j,size;
  cid_t me,other,m;
  tid_t myid,o;
  wgt_t tiw, tew, tie;
  vtx_t * myadjnum;
  cid_t * myadjncy;
  wgt_t * myadjwgt;
  vtx_iset_t * mybnd;
  ucinfo_t * myucinfo;
  nbrinfo_t * myinfo;

  /* expose the graph parts */
  const graph_t * const graph = mgraph->graph;
  const tid_t nthreads = graph->npar;
  const adj_t * const * const xadj = (const adj_t * const *)graph->xadj;
  const vtx_t * const * const adjncy = (const vtx_t * const *)graph->adjncy;
  const wgt_t * const * const adjwgt = (const wgt_t * const *)graph->adjwgt;
  const wgt_t * const * const iadjwgt = (const wgt_t * const *)graph->iadjwgt;

  /* expose clustering parts */
  const cid_t * const * const where = (const cid_t * const *)clustering->where;
  const cid_t nclusters = clustering->nclusters;

  /* the ucinfo */
  ucinfo_t ** const ucinfo = r_ucinfo_sym_calloc(1,nthreads); 
  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    myucinfo = ucinfo[myid];
    myucinfo->nbrinfo = nbrinfo_alloc(mynvtxs);
    myucinfo->maxnadj = dl_max(graph->mynedges[myid]*4,MIN_NADJ);
    myucinfo->nadj = 0;
    myucinfo->adjncy = cid_alloc(myucinfo->maxnadj);
    myucinfo->adjwgt = wgt_alloc(myucinfo->maxnadj);
    myucinfo->adjnum = vtx_alloc(myucinfo->maxnadj);
    myucinfo->free_adjncy = 1;
    myucinfo->free_adjwgt = 1;
    myucinfo->free_adjnum = 1;
    myucinfo->bnd = vtx_iset_create(0,mynvtxs);
  }

  /* other stuff */
  cid_t * htable = cid_init_alloc(NULL_VTX,nclusters);

  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    myucinfo = ucinfo[myid];
    mybnd = myucinfo->bnd;
    /* parse the vertex list */
    for (i=0;i<mynvtxs;++i) {
      tin = 0;
      tiw = tew = 0.0;
      me = where[myid][i];
      size = dl_min(xadj[myid][i+1]-xadj[myid][i],nclusters);
      myinfo = myucinfo->nbrinfo + i;
      myinfo->nnbrs = 0;
      myinfo->nbrstart = next_nbr_adj(myucinfo,size);
      myadjncy = myucinfo->adjncy+myinfo->nbrstart;
      myadjwgt = myucinfo->adjwgt+myinfo->nbrstart;
      myadjnum = myucinfo->adjnum+myinfo->nbrstart;
      /* find external and internal vertices */
      for (j=xadj[myid][i];j<xadj[myid][i+1];++j) {
        k = adjncy[myid][j];
        if (k < mynvtxs) {
          o = myid;
          l = k;
        } else {
          o = gvtx_to_tid(k,graph->dist);
          l = gvtx_to_lvtx(k,graph->dist);
        }
        other = where[o][l];
        if (me != other) {
          tew += adjwgt[myid][j]; 
          if ((m = htable[other]) == NULL_VTX) {
            myadjncy[myinfo->nnbrs] = other;
            myadjwgt[myinfo->nnbrs] = adjwgt[myid][j];
            myadjnum[myinfo->nnbrs] = 1;
            htable[other] = myinfo->nnbrs++;
          } else {
            myadjwgt[m] += adjwgt[myid][j];
            ++myadjnum[m];
          }
        } else {
          tiw += adjwgt[myid][j];
          ++tin;
        }
      }
      myinfo->in = tin;
      myinfo->iw = tiw;
      myinfo->ew = tew;
      if (tew == 0) {
        release_nbr_adj(myucinfo,myinfo,size);
      } else {
        if (iadjwgt) {
          tie = iadjwgt[myid][i];
        } else {
          tie = 0;
        }
        if (is_bnd(tiw,tew,tie)) {
          vtx_iset_add(i,mybnd);
        }
        /* reset the htable */
        for (m=0;m<myinfo->nnbrs;++m) {
          DL_ASSERT(htable[myadjncy[m]]!=NULL_VTX,
              "Clearing a clear spot in the htable!\n");
          htable[myadjncy[m]] = NULL_VTX;
        }
      }
    }
  }

  dl_free(htable);
  return ucinfo;
}


int check_bnd(
    const mgraph_t * const mgraph, 
    const cid_t * const * const where)
{
  vtx_t i, v, u, n;
  adj_t j;
  wgt_t tiw, tew, tie, ewgt;
  tid_t myid, o;

  const vtx_iset_t * bnd;
  vtx_t mynvtxs;
  
  const ucinfo_t * const * const ucinfo = 
      (const ucinfo_t * const *)mgraph->ucinfo;

  const graph_t * const graph = mgraph->graph;
  const tid_t nthreads = graph->npar;
  const adj_t * const * const xadj = (const adj_t * const *)graph->xadj;
  const vtx_t * const * const adjncy = (const vtx_t * const *)graph->adjncy;
  const wgt_t * const * const adjwgt = (const wgt_t * const *)graph->adjwgt; 
  const wgt_t * const * const iadjwgt = (const wgt_t * const *)graph->iadjwgt; 

  for (myid=0;myid<nthreads;++myid) {
    bnd = ucinfo[myid]->bnd;
    mynvtxs = graph->mynvtxs[myid];
    n = 0;
    for (v=0;v<mynvtxs;++v) {
      i = vtx_iset_indexof(v,bnd);
      if (i != NULL_VTX) {
        ++n;
      }
      tiw = tew = 0.0;
      for (j=xadj[myid][v];j<xadj[myid][v+1];++j) {
        u = adjncy[myid][j];
        if (u < graph->mynvtxs[myid]) {
          o = myid;
        } else {
          o = gvtx_to_tid(u,graph->dist);
          u = gvtx_to_lvtx(u,graph->dist);
        }
        if (adjwgt) {
          ewgt = adjwgt[myid][j];
        } else {
          ewgt = 1.0;
        }
        if (where[myid][v] == where[o][u]) {
          tiw += ewgt;
        } else {
          tew += ewgt;
        }
      }
      if (iadjwgt) {
        tie = iadjwgt[myid][v];
      } else {
        tie = 0;
      }
      if (is_bnd(tiw,tew,tie)) {
        if (!vtx_iset_contains(v,bnd)) {
          eprintf("Vertex "PF_VTX_T" is not in the boundary but has a total "
              "interior edge weight of "PF_WGT_T" and an exterior weight of "
              PF_WGT_T"\n",v,tiw,tew);
          return 0;
        }
      } else {
        if (vtx_iset_contains(v,bnd)) {
          eprintf("Vertex "PF_VTX_T" is in the boundary but has a total "
              "interior edge weight of "PF_WGT_T" and an exterior wegiht of "
              PF_WGT_T"\n",v,tiw,tew);
          return 0;
        }
      }
    }

    if (n != (vtx_t)bnd->size) {
      eprintf("Counted "PF_VTX_T" bnd vertices, but nbnd = "PF_VTX_T"\n",n,
          (vtx_t)bnd->size);
      return 0;
    }
  }

  return 1;
}




/******************************************************************************
* PUBLIC PARALLEL FUNCTIONS ***************************************************
******************************************************************************/


clustering_t * par_uncoarsen_graph(
    objective_t * const objective, 
    clustering_t * clustering, 
    mgraph_t * const fmgraph, 
    const mgraph_t * const cmgraph)
{
  const tid_t myid = dlthread_get_id(objective->comm);

  if (myid == 0) {
    dl_start_timer(&(objective->timers.uncoarsening));
  }

  DL_ASSERT(check_graph(cmgraph->graph) == 1, "Bad graph passed to "
      "uncoarsen_graph\n");
  DL_ASSERT(check_clustering(clustering,cmgraph->graph) == 1,
      "Bad clustering at start of uncoarsen_graph\n");
  DL_ASSERT(check_ucinfo(cmgraph->graph,
        (const ucinfo_t * const *)cmgraph->ucinfo,
        (const cid_t * const *)clustering->where,clustering->nclusters) == 1, 
      "Bad ucinfo passed to uncoarsen_graph\n");
  DL_ASSERT(check_bnd(cmgraph,(const cid_t * const *)clustering->where) == 1,
      "Bad boundary at start of uncoarsen_graph\n");

  clustering = par_project_clustering(objective,clustering,fmgraph,cmgraph);

  DL_ASSERT(check_clustering(clustering,fmgraph->graph) == 1,
      "Bad clustering after projection\n");
  DL_ASSERT(check_ucinfo(fmgraph->graph,
        (const ucinfo_t * const *)fmgraph->ucinfo,
        (const cid_t * const *)clustering->where,clustering->nclusters) == 1, 
      "Bad ucinfo after projection\n");
  DL_ASSERT(check_bnd(fmgraph,(const cid_t * const *)clustering->where) == 1,
      "Bad boundary after projection\n");

  clustering = par_refine_graph(objective,clustering,fmgraph);

  DL_ASSERT(check_clustering(clustering,fmgraph->graph) == 1,
      "Bad clustering after refinement\n");
  DL_ASSERT(check_ucinfo(fmgraph->graph,
        (const ucinfo_t * const *)fmgraph->ucinfo,
        (const cid_t * const *)clustering->where,clustering->nclusters) == 1, 
      "Bad ucinfo after refinement\n");
  DL_ASSERT(check_bnd(fmgraph,(const cid_t * const *)clustering->where) == 1,
      "Bad boundary after refinement\n");

  if (myid == 0) {
    dl_stop_timer(&(objective->timers.uncoarsening));
  }

  return clustering;
}


ucinfo_t * par_derive_ucinfo(
    const clustering_t * const clustering, 
    const mgraph_t * const mgraph,
    dlthread_comm_t comm)
{
  vtx_t i,k,l,tin;
  adj_t j,size;
  cid_t me,other,m;
  tid_t o;
  wgt_t tiw, tew, tie;
  vtx_t * myadjnum;
  cid_t * myadjncy;
  wgt_t * myadjwgt;
  nbrinfo_t * myinfo;

  const tid_t myid = dlthread_get_id(comm);

  /* expose the graph parts */
  const graph_t * const graph = mgraph->graph;
  const vtx_t mynvtxs = graph->mynvtxs[myid];
  const adj_t mynedges = graph->mynedges[myid];
  const tid_t nthreads = graph->npar;
  const adj_t * const xadj = graph->xadj[myid];
  const vtx_t * const adjncy = graph->adjncy[myid];
  const wgt_t * const adjwgt = graph->adjwgt[myid];
  const wgt_t * const iadjwgt = graph->iadjwgt[myid];

  /* expose clustering parts */
  const cid_t * const * const where = (const cid_t * const *)clustering->where;
  const cid_t nclusters = clustering->nclusters;

  /* the ucinfo */
  static ucinfo_t ** ucinfos;
  if (myid == 0) {
    ucinfos = r_ucinfo_alloc(nthreads); 
  }
  dlthread_barrier(comm);

  ucinfo_t * const myucinfo = ucinfos[myid] = ucinfo_calloc(1);
  myucinfo->nbrinfo = nbrinfo_alloc(mynvtxs);
  myucinfo->maxnadj = dl_max(mynedges*2,MIN_NADJ);
  myucinfo->nadj = 0;
  myucinfo->adjncy = cid_alloc(myucinfo->maxnadj);
  myucinfo->adjwgt = wgt_alloc(myucinfo->maxnadj);
  myucinfo->adjnum = vtx_alloc(myucinfo->maxnadj);
  myucinfo->free_adjncy = 1;
  myucinfo->free_adjwgt = 1;
  myucinfo->free_adjnum = 1;
  vtx_iset_t * const mybnd = myucinfo->bnd = vtx_iset_create(0,mynvtxs);

  /* other stuff */
  cid_t * htable = cid_init_alloc(NULL_VTX,nclusters);

  /* parse the vertex list */
  for (i=0;i<mynvtxs;++i) {
    tin = 0;
    tiw = tew = 0.0;
    me = where[myid][i];
    size = dl_min(xadj[i+1]-xadj[i],nclusters);
    myinfo = myucinfo->nbrinfo + i;
    myinfo->nnbrs = 0;
    myinfo->nbrstart = next_nbr_adj(myucinfo,size);
    myadjncy = myucinfo->adjncy+myinfo->nbrstart;
    myadjwgt = myucinfo->adjwgt+myinfo->nbrstart;
    myadjnum = myucinfo->adjnum+myinfo->nbrstart;
    /* find external and internal vertices */
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (k < mynvtxs) {
        o = myid;
        l = k;
      } else {
        o = gvtx_to_tid(k,graph->dist);
        l = gvtx_to_lvtx(k,graph->dist);
      }
      other = where[o][l];
      if (me != other) {
        tew += adjwgt[j]; 
        if ((m = htable[other]) == NULL_VTX) {
          myadjncy[myinfo->nnbrs] = other;
          myadjwgt[myinfo->nnbrs] = adjwgt[j];
          myadjnum[myinfo->nnbrs] = 1;
          htable[other] = myinfo->nnbrs++;
        } else {
          myadjwgt[m] += adjwgt[j];
          ++myadjnum[m];
        }
      } else {
        tiw += adjwgt[j];
        ++tin;
      }
    }
    myinfo->in = tin;
    myinfo->iw = tiw;
    myinfo->ew = tew;
    if (tew == 0) {
      release_nbr_adj(myucinfo,myinfo,size);
    } else {
      if (iadjwgt) {
        tie = iadjwgt[i];
      } else {
        tie = 0;
      }
      if (is_bnd(tiw,tew,tie)) {
        vtx_iset_add(i,mybnd);
      }
      /* reset the htable */
      for (m=0;m<myinfo->nnbrs;++m) {
        DL_ASSERT(htable[myadjncy[m]]!=NULL_VTX,
            "Clearing a clear spot in the htable!\n");
        htable[myadjncy[m]] = NULL_VTX;
      }
    }
  }

  dl_free(htable);

  return myucinfo;
}




#endif
