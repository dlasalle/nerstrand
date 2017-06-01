/**
 * @file project.c
 * @brief Functions for projecting a clustering
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2014-03-03
 */




#ifndef NERSTRAND_PROJECT_C
#define NERSTRAND_PROJECT_C




#include "project.h"




/******************************************************************************
* PRIVATE SERIAL FUNCTIONS ****************************************************
******************************************************************************/


/**
 * @brief Project a clustering from coarse graph to the fine graph one level
 * down directly (using cmap) 
 *
 * @param clustering The clustering to project
 * @param fmgraph The fine graph
 * @param cmgraph The coarse graph
 *
 * @return The updated clustering
 */
static clustering_t * __project_DIRECT(
    clustering_t * clustering, 
    mgraph_t * const fmgraph, 
    mgraph_t const * const cmgraph)
{
  vtx_t i,v,k,l,mynvtxs, tin;
  adj_t j,size;
  cid_t c,me,other; 
  wgt_t tiw, tew, tie;
  tid_t myid, o;
  nbrinfo_t * myinfo;
  ucinfo_t * myucinfo;
  vtx_iset_t * mybnd;

  DL_ASSERT_EQUALS(clustering->npar,fmgraph->graph->npar,PF_TID_T);
  DL_ASSERT_EQUALS(clustering->npar,cmgraph->graph->npar,PF_TID_T);

  /* cluster stuff */
  cid_t const nclusters = clustering->nclusters;
  tid_t const nthreads = clustering->npar;
  cluster_t * const clusters = clustering->clusters;
  cid_t ** cwhere = clustering->where;
  vtx_t ** ccounts = r_vtx_sym_calloc(nclusters,nthreads);
  wgt_t ** viadjwgt = r_wgt_sym_calloc(nclusters,nthreads);
  /* eadjwgt and vwgt do not change in projection */

  /* project specific */
  cid_t * htable = cid_init_alloc(NULL_CID,nclusters);

  /* coarse graph stuff */
  ucinfo_t ** const cucinfo = cmgraph->ucinfo;
  graph_t const * const cgraph = cmgraph->graph;

  /* expose fine graph stuff */
  graph_t const * const graph = fmgraph->graph;
  adj_t const * const * const xadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const * const *)graph->adjncy;
  wgt_t const * const * const adjwgt = (wgt_t const * const *)graph->adjwgt;
  wgt_t const * const * const iadjwgt = (wgt_t const * const *)graph->iadjwgt;
  cid_t ** const where = r_cid_dalloc(graph->mynvtxs,sizeof(vtx_t),nthreads);
  vtx_t * const * const cmap = fmgraph->cmap; 

  /* allocate fine ucinfo and friends */
  ucinfo_t ** const ucinfo = fmgraph->ucinfo = r_ucinfo_sym_calloc(1,nthreads);
  fmgraph->free_ucinfo = 1;
  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    myucinfo = ucinfo[myid];
    myucinfo->bnd = vtx_iset_create(0,mynvtxs);
    myucinfo->nbrinfo = nbrinfo_alloc(mynvtxs);
    myucinfo->maxnadj = cucinfo[myid]->maxnadj;
    myucinfo->nadj = cucinfo[myid]->nadj;
    myucinfo->adjncy = cucinfo[myid]->adjncy;
    myucinfo->adjwgt = cucinfo[myid]->adjwgt; 
    myucinfo->adjnum = cucinfo[myid]->adjnum;
    cucinfo[myid]->free_adjncy = 0;
    cucinfo[myid]->free_adjwgt = 0;
    cucinfo[myid]->free_adjnum = 0;
    myucinfo->free_adjncy = 1;
    myucinfo->free_adjwgt = 1;
    myucinfo->free_adjnum = 1;
  }

  /* project where */
  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    for (i=0;i<mynvtxs;++i) {
      v = cmap[myid][i]; 
      if (v < cgraph->mynvtxs[myid]) {
        o = myid;
      } else {
        o = gvtx_to_tid(v,cgraph->dist);
        DL_ASSERT(o < nthreads,"Invalid thread id "PF_TID_T"/"PF_TID_T
            " derived from vertex "PF_VTX_T"/"PF_VTX_T"\n",o,nthreads,v,
            cgraph->nvtxs);
        v = gvtx_to_lvtx(v,cgraph->dist); 
      }
      DL_ASSERT(v < cgraph->mynvtxs[o],"Invalid local vertex id "PF_VTX_T"/"
          PF_VTX_T"\n",v,cgraph->mynvtxs[o]);
      me = where[myid][i] = cwhere[o][v];
      DL_ASSERT(me < nclusters,"Invalid cluster value found in cwhere["
          PF_TID_T"]["PF_VTX_T"] = "PF_CID_T"/"PF_CID_T"\n",o,v,me,nclusters);
      ++ccounts[myid][me];
    }
  }

  /* project borders */
  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid]; 
    myucinfo = ucinfo[myid];
    mybnd = myucinfo->bnd;
    for (i=0;i<mynvtxs;++i) {
      me = where[myid][i]; 
      myinfo = ucinfo[myid]->nbrinfo+i;
      myinfo->nnbrs = 0;
      v = cmap[myid][i]; 
      if (v < cgraph->mynvtxs[myid]) {
        o = myid;
      } else {
        o = gvtx_to_tid(v,cgraph->dist);
        v = gvtx_to_lvtx(v,cgraph->dist); 
      }
      if (cucinfo[o]->nbrinfo[v].nnbrs == 0) { /* non border node */
        myinfo->in = xadj[myid][i+1] - xadj[myid][i];
        myinfo->iw = graph->eadjwgt[myid][i];      
        myinfo->ew = 0;
        myinfo->nbrstart = NO_NBRS;
      } else { /* possibly a border node */
        /* using clustering->nclusters as a max could be a problem for variable
           number of clusters in the future */
        size = dl_min(xadj[myid][i+1]-xadj[myid][i],nclusters);
        myinfo->nbrstart = next_nbr_adj(myucinfo,size);

        tiw = 0;
        tew = 0;
        tin = 0;
        /* calculate neighbor information */
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
          if (me == other) {
            if (adjwgt) {
              tiw += adjwgt[myid][j];
            } else {
              tiw += 1.0;
            }
            ++tin;
          } else {
            if (adjwgt) {
              tew += adjwgt[myid][j];
            } else {
              tew += 1.0;
            }
            if ((c = htable[other]) == NULL_CID) {
              c = htable[other] = myinfo->nnbrs;
              myucinfo->adjncy[myinfo->nbrstart+c] = other;
              if (adjwgt) {
                myucinfo->adjwgt[myinfo->nbrstart+c] = adjwgt[myid][j];
              } else {
                myucinfo->adjwgt[myinfo->nbrstart+c] = 1;
              }
              myucinfo->adjnum[myinfo->nbrstart+c] = 1;
              ++myinfo->nnbrs;
            } else {
              if (adjwgt) {
                myucinfo->adjwgt[myinfo->nbrstart+c] += adjwgt[myid][j];
              } else {
                myucinfo->adjwgt[myinfo->nbrstart+c] += 1.0;
              }
              ++myucinfo->adjnum[myinfo->nbrstart+c];
            }
          }
        }
        myinfo->in = tin;
        myinfo->iw = tiw;
        myinfo->ew = tew;

        /* cleanup unused space */
        if (myinfo->nnbrs == 0) {
          release_nbr_adj(myucinfo,myinfo,size);
        } else {
          /* reset the htable */
          for (c=0;c<myinfo->nnbrs;++c) {
            DL_ASSERT(htable[myucinfo->adjncy[myinfo->nbrstart+c]]!=NULL_CID,
                "Clearing a clear spot in the htable!\n");
            htable[myucinfo->adjncy[myinfo->nbrstart+c]] = NULL_CID;
          }
          if (iadjwgt) {
            tie = iadjwgt[myid][i];
          } else {
            tie = 0;
          }
          if (is_bnd(myinfo->iw,myinfo->ew,tie)) {
            vtx_iset_add(i,mybnd);
          }
        }
      }
    }
    if (iadjwgt) {
      for (i=0;i<mynvtxs;++i) {
        viadjwgt[myid][where[myid][i]] += iadjwgt[myid][i];
      }
    }
  }
  dl_free(htable);

  /* out with the old, in with the new */
  r_cid_free(cwhere,nthreads);
  clustering->where = where;

  for (c=0;c<nclusters;++c) {
    for (myid=1;myid<nthreads;++myid) {
      ccounts[0][c] += ccounts[myid][c];
      viadjwgt[0][c] += viadjwgt[myid][c];
    }
    clusters[c].nvtxs = ccounts[0][c];
    /* all vertex internal edgeweights that disappear become clsuter internal
     * edge weights -- clsuter external edge weights are unchanged */
    clusters[c].iadjwgt = dl_max(0, \
        clusters[c].iadjwgt + clusters[c].viadjwgt - viadjwgt[0][c]);
    clusters[c].viadjwgt = viadjwgt[0][c];
  }

  r_vtx_free(ccounts,nthreads);
  r_wgt_free(viadjwgt,nthreads);

  DL_ASSERT(check_ucinfo(fmgraph->graph,
        (ucinfo_t const * const *)fmgraph->ucinfo,
        (cid_t const * const *)where,nclusters) == 1,
      "Bad ucinfo from projection\n");

  return clustering;
}




/******************************************************************************
* PRIVATE PARALLEL FUNCTIONS **************************************************
******************************************************************************/


/**
 * @brief Project a clustering from a coarse graph to the fine graph one level
 * down directly (using cmap) -- to be called by each thread from within a
 * parallel region
 *
 * @param clustering The clustering to project
 * @param fmgraph The fine graph
 * @param cmgraph The coarse graph
 *
 * @return The updated clustering
 */
static clustering_t * __par_project_DIRECT(
    clustering_t * clustering, 
    mgraph_t * const fmgraph, 
    mgraph_t const * const cmgraph,
    dlthread_comm_t const comm) 
{
  vtx_t i,v,k,l, tnvtxs, tin;
  adj_t j,size;
  cid_t c,me,other,cblock,cstart,cend; 
  wgt_t tiw, tew, tviadjwgt, tie;
  tid_t o;
  nbrinfo_t * myinfo;

  tid_t const nthreads = dlthread_get_nthreads(comm);
  tid_t const myid = dlthread_get_id(comm);

  DL_ASSERT_EQUALS(clustering->npar,fmgraph->graph->npar,PF_TID_T);
  DL_ASSERT_EQUALS(clustering->npar,cmgraph->graph->npar,PF_TID_T);

  /* cluster stuff */
  cid_t const nclusters = clustering->nclusters;
  cluster_t * const clusters = clustering->clusters;
  cid_t ** cwhere = clustering->where;
  vtx_t * ccounts = vtx_calloc(nclusters);
  wgt_t * viadjwgt = wgt_calloc(nclusters);
  /* eadjwgt and vwgt do not change in projection */

  /* project specific */
  cid_t * htable = cid_init_alloc(NULL_CID,nclusters);

  /* coarse graph stuff */
  ucinfo_t ** const cucinfo = cmgraph->ucinfo;
  graph_t const * const cgraph = cmgraph->graph;
  vtx_t const mycnvtxs = cgraph->mynvtxs[myid];

  /* expose fine graph stuff */
  graph_t const * const graph = fmgraph->graph;
  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const * const xadj = graph->xadj[myid];
  vtx_t const * const adjncy = graph->adjncy[myid];
  wgt_t const * const adjwgt = (graph->adjwgt ? graph->adjwgt[myid] : NULL);
  wgt_t const * const eadjwgt = (graph->eadjwgt ? graph->eadjwgt[myid] : NULL);
  wgt_t const * const iadjwgt = (graph->iadjwgt ? graph->iadjwgt[myid] : NULL);
  vtx_t * const cmap = fmgraph->cmap[myid]; 

  /* allocate fine ucinfo and friends */
  dlthread_barrier(comm);
  if (myid == 0) {
    fmgraph->ucinfo = r_ucinfo_sym_calloc(1,nthreads);
    fmgraph->free_ucinfo = 1;

    clustering->where = r_cid_alloc(nthreads);
  }
  dlthread_barrier(comm);

  cid_t ** const where = clustering->where;
  where[myid] = cid_alloc(mynvtxs);
  ucinfo_t * const myucinfo = fmgraph->ucinfo[myid];
  vtx_iset_t * const mybnd = myucinfo->bnd = vtx_iset_create(0,mynvtxs);

  myucinfo->nbrinfo = nbrinfo_alloc(mynvtxs);
  myucinfo->maxnadj = cucinfo[myid]->maxnadj;
  myucinfo->nadj = cucinfo[myid]->nadj;
  myucinfo->adjncy = cucinfo[myid]->adjncy;
  myucinfo->adjwgt = cucinfo[myid]->adjwgt; 
  myucinfo->adjnum = cucinfo[myid]->adjnum;
  cucinfo[myid]->free_adjncy = 0;
  cucinfo[myid]->free_adjwgt = 0;
  cucinfo[myid]->free_adjnum = 0;
  myucinfo->free_adjncy = 1;
  myucinfo->free_adjwgt = 1;
  myucinfo->free_adjnum = 1;

  /* project where */
  for (i=0;i<mynvtxs;++i) {
    v = cmap[i]; 
    if (v < cgraph->mynvtxs[myid]) {
      o = myid;
    } else {
      o = gvtx_to_tid(v,cgraph->dist);
      DL_ASSERT(o < nthreads,"Invalid thread id "PF_TID_T"/"PF_TID_T
          " derived from vertex "PF_VTX_T"/"PF_VTX_T"\n",o,nthreads,v,
          cgraph->nvtxs);
      v = gvtx_to_lvtx(v,cgraph->dist); 
    }
    DL_ASSERT(v < cgraph->mynvtxs[o],"Invalid local vertex id "PF_VTX_T"/"
        PF_VTX_T"\n",v,cgraph->mynvtxs[o]);
    me = where[myid][i] = cwhere[o][v];
    DL_ASSERT(me < nclusters,"Invalid cluster value found in cwhere["
        PF_TID_T"]["PF_VTX_T"] = "PF_CID_T"/"PF_CID_T"\n",o,v,me,nclusters);
    ++ccounts[me];
  }

  dlthread_barrier(comm);
  dl_free(cwhere[myid]);
  dlthread_barrier(comm);
  if (myid == 0) {
    dl_free(cwhere);
  }

  /* project borders */
  for (i=0;i<mynvtxs;++i) {
    me = where[myid][i]; 
    myinfo = myucinfo->nbrinfo+i;
    myinfo->nnbrs = 0;
    v = cmap[i]; 
    if (v < mycnvtxs) {
      o = myid;
    } else {
      o = gvtx_to_tid(v,cgraph->dist);
      v = gvtx_to_lvtx(v,cgraph->dist); 
    }
    if (cucinfo[o]->nbrinfo[v].nnbrs == 0) { /* non border node */
      myinfo->in = xadj[i+1] - xadj[i];
      myinfo->iw = eadjwgt[i];      
      myinfo->ew = 0;
      myinfo->nbrstart = NO_NBRS;
    } else { /* possibly a border node */
      /* using clustering->nclusters as a max could be a problem for variable
         number of clusters in the future */
      size = dl_min(xadj[i+1]-xadj[i],nclusters);
      myinfo->nbrstart = next_nbr_adj(myucinfo,size);

      tiw = 0;
      tew = 0;
      tin = 0;
      /* calculate neighbor information */
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
        if (me == other) {
          if (adjwgt) {
            tiw += adjwgt[j];
          } else {
            tiw += 1.0;
          }
          ++tin;
        } else {
          if (adjwgt) {
            tew += adjwgt[j];
          } else {
            tew += 1.0;
          }
          if ((c = htable[other]) == NULL_CID) {
            c = htable[other] = myinfo->nnbrs;
            myucinfo->adjncy[myinfo->nbrstart+c] = other;
            if (adjwgt) {
              myucinfo->adjwgt[myinfo->nbrstart+c] = adjwgt[j];
            } else {
              myucinfo->adjwgt[myinfo->nbrstart+c] = 1;
            }
            myucinfo->adjnum[myinfo->nbrstart+c] = 1;
            ++myinfo->nnbrs;
          } else {
            if (adjwgt) {
              myucinfo->adjwgt[myinfo->nbrstart+c] += adjwgt[j];
            } else {
              myucinfo->adjwgt[myinfo->nbrstart+c] += 1.0;
            }
            ++myucinfo->adjnum[myinfo->nbrstart+c];
          }
        }
      }
      myinfo->in = tin;
      myinfo->iw = tiw;
      myinfo->ew = tew;

      /* cleanup unused space */
      if (myinfo->nnbrs == 0) {
        release_nbr_adj(myucinfo,myinfo,size);
      } else {
        /* reset the htable */
        for (c=0;c<myinfo->nnbrs;++c) {
          DL_ASSERT(htable[myucinfo->adjncy[myinfo->nbrstart+c]]!=NULL_CID,
              "Clearing a clear spot in the htable!\n");
          htable[myucinfo->adjncy[myinfo->nbrstart+c]] = NULL_CID;
        }
        if (iadjwgt) {
          tie = iadjwgt[i];
        } else {
          tie = 0;
        }
        if (is_bnd(myinfo->iw,myinfo->ew,tie)) {
          vtx_iset_add(i,mybnd);
        }
      }
    }
  }
  if (iadjwgt) {
    for (i=0;i<mynvtxs;++i) {
      viadjwgt[where[myid][i]] += iadjwgt[i];
    }
  }
  dl_free(htable);

  /* this is poorly formulated for large numbers of clusters */
  vtx_dlthread_sumareduce(ccounts,nclusters,comm);
  wgt_dlthread_sumareduce(viadjwgt,nclusters,comm);

  /* have each thread process a chunk of 64 clusters */
  for (cblock=myid;cblock*64<nclusters;cblock+=nthreads) {
    cstart = cblock*64;
    cend = dl_min((cblock+1)*64,nclusters);
    for (c=cstart;c<cend;++c) {
      tnvtxs = ccounts[c];
      tviadjwgt = viadjwgt[c];
      clusters[c].nvtxs = tnvtxs;
      /* all vertex internal edgeweights that disappear become clsuter 
       * internal edge weights -- clsuter external edge weights are 
       * unchanged */
      clusters[c].iadjwgt = dl_max(0, \
          clusters[c].iadjwgt + clusters[c].viadjwgt - tviadjwgt);
      clusters[c].viadjwgt = tviadjwgt;
    }
  }

  dl_free(ccounts);
  dl_free(viadjwgt);

  dlthread_barrier(comm);

  DL_ASSERT(check_ucinfo(fmgraph->graph,
        (ucinfo_t const * const *)fmgraph->ucinfo,
        (cid_t const * const *)where,nclusters) == 1,
      "Bad ucinfo from projection\n");

  return clustering;
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


clustering_t * project_clustering(
    objective_t * const objective, 
    clustering_t * clustering, 
    mgraph_t * const fmgraph, 
    mgraph_t const * const cmgraph)
{
  dl_start_timer(&(objective->timers.projection));
  DL_ASSERT_EQUALS(clustering->npar,fmgraph->graph->npar,PF_TID_T);
  DL_ASSERT_EQUALS(clustering->npar,cmgraph->graph->npar,PF_TID_T);

  dprintf("Performing projection from level "PF_SIZE_T" to level "PF_SIZE_T
      "\n",cmgraph->level, fmgraph->level);

  switch(objective->projtype) {
    case NERSTRAND_PROJECT_DIRECT:
      clustering = __project_DIRECT(clustering,fmgraph,cmgraph);
      break;
    default:
      dl_error("Unknown projection type\n");
      break;
  }

  dl_stop_timer(&(objective->timers.projection));

  return clustering;
}




/******************************************************************************
* PUBLIC PARALELL FUNCTIONS ***************************************************
******************************************************************************/


clustering_t * par_project_clustering(
    objective_t * const objective, 
    clustering_t * clustering, 
    mgraph_t * const fmgraph, 
    mgraph_t const * const cmgraph)
{
  const tid_t myid = dlthread_get_id(objective->comm);

  if (myid == 0) {
    dl_start_timer(&(objective->timers.projection));
  }

  DL_ASSERT_EQUALS(clustering->npar,fmgraph->graph->npar,PF_TID_T);
  DL_ASSERT_EQUALS(clustering->npar,cmgraph->graph->npar,PF_TID_T);

  if (myid == 0) {
    dprintf("Performing projection from level "PF_SIZE_T" to level "PF_SIZE_T
       "\n",cmgraph->level, fmgraph->level);
  }

  switch(objective->projtype) {
    case NERSTRAND_PROJECT_DIRECT:
      clustering = __par_project_DIRECT(clustering,fmgraph,cmgraph, \
          objective->comm);
      break;
    default:
      dl_error("Unknown projection type\n");
      break;
  }

  if (myid == 0) {
    dl_stop_timer(&(objective->timers.projection));
  }

  return clustering;
}




#endif
