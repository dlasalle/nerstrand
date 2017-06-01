/**
 * @file multilevel.c
 * @brief Toplevel functions for multilevel clustering
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-10-08
 */




#ifndef NERSTRAND_MULTILEVEL_C
#define NERSTRAND_MULITLEVEL_C




#include "multilevel.h"




/******************************************************************************
* PRIVATE SERIAL FUNCTIONS ****************************************************
******************************************************************************/


static int __squash_graphs(
    mgraph_t * const cmgraph, 
    mgraph_t * mgraph,
    mgraph_t * const fmgraph) 
{
  tid_t myid, o;
  vtx_t i, k, c, mynvtxs;
  vtx_t ** fcmap, ** cmap;

  const tid_t nthreads = cmgraph->graph->npar;

  fcmap = fmgraph->cmap;
  cmap = mgraph->cmap;

  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = fmgraph->graph->mynvtxs[myid];
    for (i=0;i<mynvtxs;++i) {
      k = fcmap[myid][i];
      if (k < mgraph->graph->mynvtxs[myid]) {
        c = k;
        o = myid;
      } else {
        c = gvtx_to_lvtx(k,mgraph->graph->dist); 
        o = gvtx_to_tid(k,mgraph->graph->dist); 
      }
      fcmap[myid][i] = cmap[o][c];
    }
  }

  cmgraph->level = mgraph->level;
  cmgraph->finer = fmgraph;
  fmgraph->coarser = cmgraph;

  mgraph->free_cmap = 1;
  mgraph->free_graph = 1;
  free_mgraph(mgraph);

  return NERSTRAND_SUCCESS;
}


static clustering_t * __cluster_graph_kway(
    objective_t * const objective,
    const graph_t * const graph)
{
  int stop;
  size_t lb_offset;
  mgraph_t * cmgraph, * fmgraph, * mgraph;
  mod_t maxmod = 0;
  mod_t oldmod;
  clustering_t * clustering;
  loadbalance_stat_t lb;
  dl_timer_t tmr;

  const vtx_t tnvtxs = objective->cnvtxs_per_cluster*objective->nclusters;
  const real_t stopratio = objective->stopratio;
  mgraph_t * const ograph = setup_mgraph(0,graph,NULL,NULL,NULL);

  if (objective->lbstats) {
    lb_offset = objective->lbstats->size;
  } else {
    lb_offset = 0;
  }

  dl_init_timer(&tmr);

  /* coarsen the graph */
  cmgraph = ograph;
  do { 
    if (objective->lbstats) {
      lb.adj_imbalance = 1.0;
      lb.vtx_imbalance = 1.0;
      lb.nbr_imbalance = lb.bnd_imbalance = 1.0;
      loadbalance_stats_push(lb,objective->lbstats);
    }

    fmgraph = cmgraph;
    vprintf(objective->verbosity,NERSTRAND_VERBOSITY_HIGH,"Coarsening graph "
        "of "PF_VTX_T" vertices and "PF_ADJ_T" edges\n",fmgraph->graph->nvtxs,
        fmgraph->graph->nedges);

    dl_reset_timer(&tmr);
    dl_start_timer(&tmr);

    cmgraph = coarsen_graph(objective,fmgraph);

    if (cmgraph->graph->nedges > fmgraph->graph->nedges*objective->restep) {
      vprintf(objective->verbosity,NERSTRAND_VERBOSITY_HIGH,"Making re-step "
          "do to small edge shrinkage %0.03lf < %0.03lf\n",
          ((double)cmgraph->graph->nedges)/fmgraph->graph->nedges,
          objective->restep);

      mgraph = cmgraph;
      cmgraph = coarsen_graph(objective,mgraph);
      __squash_graphs(cmgraph,mgraph,fmgraph);
    }

    dl_stop_timer(&tmr);

    vprintf(objective->verbosity,NERSTRAND_VERBOSITY_HIGH,"Coarsening at "
        "level %zu took %0.05lfs\n",fmgraph->level,dl_poll_timer(&tmr));

    /* figure out if we should exit */
    stop = 0;
    switch (objective->stopcondition) {
      case NERSTRAND_STOPCONDITION_EDGES:
        stop = (cmgraph->graph->nedges >= fmgraph->graph->nedges*stopratio);
        break;
      case NERSTRAND_STOPCONDITION_VERTICES:
        stop = (cmgraph->graph->nvtxs >= fmgraph->graph->nvtxs*stopratio); 
        break;
      case NERSTRAND_STOPCONDITION_SIZE:
        stop = ((3.0*cmgraph->graph->nvtxs+2.0*cmgraph->graph->nedges) >= 
            (3.0*fmgraph->graph->nvtxs + 
                 2.0*fmgraph->graph->nedges)*stopratio); 
        break;
      default:
        dl_error("Unknown stop condition: %d\n",objective->stopcondition);
    }
    if (stop || cmgraph->graph->nvtxs <= tnvtxs || 
        cmgraph->graph->nvtxs > objective->nclusters+cmgraph->graph->nedges ||
        cmgraph->level >= objective->maxlevel) {
      break;
    }
  } while (1);

  vprintf(objective->verbosity,NERSTRAND_VERBOSITY_MEDIUM,"Initial clustering "
      "being made at level "PF_SIZE_T" on a graph with "PF_VTX_T" vertices "
      "and "PF_ADJ_T" edges.\n",cmgraph->level,cmgraph->graph->nvtxs,
      cmgraph->graph->nedges);

  clustering = generate_kclustering(objective,cmgraph);

  maxmod = calc_total_modularity(clustering,cmgraph->graph);

  DL_ASSERT(check_clustering(clustering,cmgraph->graph) == 1,
      "Bad clustering generated\n");
  DL_ASSERT(check_ucinfo(cmgraph->graph,
        (const ucinfo_t * const *)cmgraph->ucinfo,
        (const cid_t * const *)clustering->where,clustering->nclusters) == 1, 
      "Bad ucinfo derived\n");
  DL_ASSERT(check_bnd(cmgraph,(const cid_t * const *)clustering->where) == 1,
      "Bad boundary after derive_ucinfo\n");

  /* project the clustering */
  while (cmgraph != ograph) {
    dl_reset_timer(&tmr);
    dl_start_timer(&tmr);

    oldmod = maxmod;
    fmgraph = cmgraph->finer;
    clustering = uncoarsen_graph(objective,clustering,fmgraph,cmgraph);
    free_mgraph(cmgraph);
    cmgraph = fmgraph;
    maxmod = calc_total_modularity(clustering,cmgraph->graph);

    dl_stop_timer(&tmr);

    if (objective->lbstats) {
      lb = loadbalance_stats_get(lb_offset+cmgraph->level,
          objective->lbstats);
      lb.nbr_imbalance = 1.0;
      lb.bnd_imbalance = 1.0;
      loadbalance_stats_replace(lb,lb_offset+cmgraph->level,
          objective->lbstats);
    }

    vprintf(objective->verbosity,NERSTRAND_VERBOSITY_HIGH,"Clustering "
        "modularity went from "PF_MOD_T" at level "PF_SIZE_T" to "PF_MOD_T
        " at level "PF_SIZE_T".\n",oldmod,cmgraph->level+1,maxmod,
        cmgraph->level);

    vprintf(objective->verbosity,NERSTRAND_VERBOSITY_HIGH,"Uncoarsening at "
        "level %zu took %0.05lfs\n",fmgraph->level,dl_poll_timer(&tmr));
  }

  /* cleanup */
  free_mgraph(ograph);

  return clustering; 
}


static clustering_t * __cluster_graph_anyway(
    objective_t * const objective,
    const graph_t * const graph)
{
  int stop;
  size_t lb_offset;
  mgraph_t * cmgraph, * fmgraph, * mgraph;
  clustering_t * clustering;
  mod_t maxmod, oldmod;
  loadbalance_stat_t lb;
  dl_timer_t tmr;

  const real_t stopratio = objective->stopratio;

  mgraph_t * const ograph = setup_mgraph(0,graph,NULL,NULL,NULL);

  if (objective->lbstats) {
    lb_offset = objective->lbstats->size;
  } else {
    lb_offset = 0;
  }

  dl_init_timer(&tmr);

  /* coarsen the graph */
  cmgraph = ograph;
  do { 
    if (objective->lbstats) {
      lb.adj_imbalance = 1.0;
      lb.vtx_imbalance = 1.0;
      lb.nbr_imbalance = lb.bnd_imbalance = 1.0;
      loadbalance_stats_push(lb,objective->lbstats);
    }

    fmgraph = cmgraph;

    vprintf(objective->verbosity,NERSTRAND_VERBOSITY_HIGH,"Coarsening graph "
        "of "PF_VTX_T" vertices and "PF_ADJ_T" edges\n",fmgraph->graph->nvtxs,
        fmgraph->graph->nedges);

    dl_reset_timer(&tmr);
    dl_start_timer(&tmr);

    cmgraph = coarsen_graph(objective,fmgraph);

    if (cmgraph->graph->nedges > fmgraph->graph->nedges*objective->restep) {
      vprintf(objective->verbosity,NERSTRAND_VERBOSITY_HIGH,"Making re-step "
          "do to small edge shrinkage %0.03lf < %0.03lf\n",
          ((double)cmgraph->graph->nedges)/fmgraph->graph->nedges,
          objective->restep);

      mgraph = cmgraph;
      cmgraph = coarsen_graph(objective,mgraph);
      __squash_graphs(cmgraph,mgraph,fmgraph);
    }

    dl_stop_timer(&tmr);

    vprintf(objective->verbosity,NERSTRAND_VERBOSITY_HIGH,"Coarsening at "
        "level %zu took %0.05lfs\n",fmgraph->level,dl_poll_timer(&tmr));

    /* figure out if we should exit */
    stop = 0;
    switch (objective->stopcondition) {
      case NERSTRAND_STOPCONDITION_EDGES:
        stop = (cmgraph->graph->nedges >= fmgraph->graph->nedges*stopratio);
        break;
      case NERSTRAND_STOPCONDITION_VERTICES:
        stop = (cmgraph->graph->nvtxs >= fmgraph->graph->nvtxs*stopratio); 
        break;
      case NERSTRAND_STOPCONDITION_SIZE:
        stop = ((3.0*cmgraph->graph->nvtxs+2.0*cmgraph->graph->nedges) >= 
            (3.0*fmgraph->graph->nvtxs+2.0*fmgraph->graph->nedges) * 
                stopratio); 
        break;
      default:
        dl_error("Unknown stop condition: %d\n",objective->stopcondition);
    }
    if (stop || cmgraph->level >= objective->maxlevel) {
      break;
    }
  } while (1);

  vprintf(objective->verbosity,NERSTRAND_VERBOSITY_MEDIUM,"Initial clustering "
      "being made at level "PF_SIZE_T" on a graph with "PF_VTX_T" vertices "
      "and "PF_ADJ_T" edges.\n",cmgraph->level,cmgraph->graph->nvtxs,
      cmgraph->graph->nedges);

  clustering = generate_aclustering(objective,cmgraph);
  maxmod = calc_total_modularity(clustering,cmgraph->graph);

  DL_ASSERT(check_clustering(clustering,cmgraph->graph) == 1,
      "Bad clustering generated\n");
  DL_ASSERT(check_ucinfo(cmgraph->graph,
        (const ucinfo_t * const *)cmgraph->ucinfo,
        (const cid_t * const *)clustering->where,clustering->nclusters) == 1, 
      "Bad ucinfo derived\n");
  DL_ASSERT(check_bnd(cmgraph,(const cid_t * const *)clustering->where) == 1,
      "Bad boundary after derive_ucinfo\n");

  /* project the clustering */
  while (cmgraph != ograph) {
    dl_reset_timer(&tmr);
    dl_start_timer(&tmr);

    oldmod = maxmod;
    fmgraph = cmgraph->finer;
    clustering = uncoarsen_graph(objective,clustering,fmgraph,cmgraph);
    free_mgraph(cmgraph);
    cmgraph = fmgraph;
    maxmod = calc_total_modularity(clustering,cmgraph->graph);

    dl_stop_timer(&tmr);

    if (objective->lbstats) {
      lb = loadbalance_stats_get(lb_offset+cmgraph->level,
          objective->lbstats);
      lb.nbr_imbalance = 1.0;
      lb.bnd_imbalance = 1.0;
      loadbalance_stats_replace(lb,lb_offset+cmgraph->level,
          objective->lbstats);
    }
    vprintf(objective->verbosity,NERSTRAND_VERBOSITY_HIGH,"Clustering "
        "modularity went from "PF_MOD_T" at level "PF_SIZE_T" to "PF_MOD_T
        " at level "PF_SIZE_T".\n",oldmod,cmgraph->level+1,maxmod,
        cmgraph->level);

    vprintf(objective->verbosity,NERSTRAND_VERBOSITY_HIGH,"Uncoarsening at "
        "level %zu took %0.05lfs\n",fmgraph->level,dl_poll_timer(&tmr));
  }

  /* cleanup */
  free_mgraph(ograph);

  return clustering; 

}




/******************************************************************************
* PRIVATE PARALLEL FUNCTIONS **************************************************
******************************************************************************/


static int __par_squash_graphs(
    mgraph_t * const cmgraph, 
    mgraph_t * mgraph,
    mgraph_t * const fmgraph,
    dlthread_comm_t comm)
{
  tid_t o;
  vtx_t i, k, c;
  vtx_t * fcmap, ** cmap;

  const tid_t myid = dlthread_get_id(comm);
  const vtx_t mynvtxs = fmgraph->graph->mynvtxs[myid];

  fcmap = fmgraph->cmap[myid];
  cmap = mgraph->cmap;

  for (i=0;i<mynvtxs;++i) {
    k = fcmap[i];
    if (k < mgraph->graph->mynvtxs[myid]) {
      c = k;
      o = myid;
    } else {
      c = gvtx_to_lvtx(k,mgraph->graph->dist); 
      o = gvtx_to_tid(k,mgraph->graph->dist); 
    }
    fcmap[i] = cmap[o][c];
  }

  dlthread_barrier(comm);
  if (myid == 0) {
    cmgraph->level = mgraph->level;
    cmgraph->finer = fmgraph;
    fmgraph->coarser = cmgraph;

    mgraph->free_cmap = 1;
    mgraph->free_graph = 1;
  }
  dlthread_barrier(comm);

  par_free_mgraph(mgraph,comm);

  return NERSTRAND_SUCCESS;
}


static clustering_t * __par_cluster_graph_kway(
    objective_t * const objective,
    const graph_t * const graph)
{
  int stop;
  real_t badj, bvtx,bbnd,bnbr;
  size_t lb_offset;
  tid_t t;
  mgraph_t * cmgraph, * fmgraph, * mgraph;
  mod_t maxmod, oldmod;
  clustering_t * clustering;
  loadbalance_stat_t lb;
  dl_timer_t tmr;

  const tid_t myid = dlthread_get_id(objective->comm);
  const tid_t nthreads = dlthread_get_nthreads(objective->comm);

  const real_t stopratio = objective->stopratio;
  const vtx_t tnvtxs = objective->cnvtxs_per_cluster*objective->nclusters;
  mgraph_t * const ograph = par_setup_mgraph(0,graph,NULL,NULL,NULL, \
      objective->comm);

  if (objective->lbstats) {
    lb_offset = objective->lbstats->size;
  } else {
    lb_offset = 0;
  }

  dl_init_timer(&tmr);

  /* coarsen the graph */
  cmgraph = ograph;
  do { 
    if (myid == 0) {
      if (objective->lbstats) {
        badj = cmgraph->graph->nedges / (real_t)nthreads;
        bvtx = cmgraph->graph->nvtxs / (real_t)nthreads;
        lb.adj_imbalance = cmgraph->graph->mynedges[0] / badj;
        lb.vtx_imbalance = cmgraph->graph->mynvtxs[0] / badj;
        lb.nbr_imbalance = lb.bnd_imbalance = 0.0;
        for (t=1;t<nthreads;++t) {
          dl_storemax(lb.adj_imbalance,cmgraph->graph->mynedges[t]/badj);
          dl_storemax(lb.vtx_imbalance,cmgraph->graph->mynvtxs[t]/bvtx);
        }
        loadbalance_stats_push(lb,objective->lbstats);
      }
    }

    fmgraph = cmgraph;
    if (myid == 0) {
      vprintf(objective->verbosity,NERSTRAND_VERBOSITY_HIGH,"Coarsening graph "
          "of "PF_VTX_T" vertices and "PF_ADJ_T" edges\n",
          fmgraph->graph->nvtxs,fmgraph->graph->nedges);

      dl_reset_timer(&tmr);
      dl_start_timer(&tmr);
    }

    cmgraph = par_coarsen_graph(objective,fmgraph);

    if (cmgraph->graph->nedges > fmgraph->graph->nedges*objective->restep) {
      if (myid == 0) {
        vprintf(objective->verbosity,NERSTRAND_VERBOSITY_HIGH,"Making re-step "
            "do to small edge shrinkage %0.03lf < %0.03lf\n",
            ((double)cmgraph->graph->nedges)/fmgraph->graph->nedges,
            objective->restep);
      }

      mgraph = cmgraph;
      cmgraph = par_coarsen_graph(objective,mgraph);
      __par_squash_graphs(cmgraph,mgraph,fmgraph,objective->comm);
    }

    if (myid == 0) {
      dl_stop_timer(&tmr);
      vprintf(objective->verbosity,NERSTRAND_VERBOSITY_HIGH,"Coarsening at "
          "level %zu took %0.05lfs\n",fmgraph->level,dl_poll_timer(&tmr));
    }

    /* figure out if we should exit */
    stop = 0;
    switch (objective->stopcondition) {
      case NERSTRAND_STOPCONDITION_EDGES:
        stop = (cmgraph->graph->nedges >= fmgraph->graph->nedges*stopratio);
        break;
      case NERSTRAND_STOPCONDITION_VERTICES:
        stop = (cmgraph->graph->nvtxs >= fmgraph->graph->nvtxs*stopratio); 
        break;
      case NERSTRAND_STOPCONDITION_SIZE:
        stop = ((3.0*cmgraph->graph->nvtxs+2.0*cmgraph->graph->nedges) >= 
            (3.0*fmgraph->graph->nvtxs+2.0*fmgraph->graph->nedges)*stopratio); 
        break;
      default:
        dl_error("Unknown stop condition: %d\n",objective->stopcondition);
    }
    if (stop || cmgraph->graph->nvtxs <= tnvtxs || 
        cmgraph->graph->nvtxs > objective->nclusters+cmgraph->graph->nedges ||
        cmgraph->level >= objective->maxlevel) {
      break;

    }
  } while (1);

  if (myid == 0) {
    vprintf(objective->verbosity,NERSTRAND_VERBOSITY_MEDIUM,"Initial "
        "clustering being made at level "PF_SIZE_T" on a graph with "PF_VTX_T
        " vertices and "PF_ADJ_T" edges.\n",cmgraph->level,
        cmgraph->graph->nvtxs,cmgraph->graph->nedges);
  }

  clustering = par_generate_kclustering(objective,cmgraph);
  maxmod = calc_total_modularity(clustering,cmgraph->graph);

  DL_ASSERT(check_clustering(clustering,cmgraph->graph) == 1,
      "Bad clustering generated\n");
  DL_ASSERT(check_ucinfo(cmgraph->graph,
        (const ucinfo_t * const *)cmgraph->ucinfo,
        (const cid_t * const *)clustering->where,clustering->nclusters) == 1, 
      "Bad ucinfo derived\n");
  DL_ASSERT(check_bnd(cmgraph,(const cid_t * const *)clustering->where) == 1,
      "Bad boundary after derive_ucinfo\n");

  /* project the clustering */
  while (cmgraph != ograph) {
    dl_reset_timer(&tmr);
    dl_start_timer(&tmr);

    oldmod = maxmod;
    fmgraph = cmgraph->finer;
    clustering = par_uncoarsen_graph(objective,clustering,fmgraph,cmgraph);
    par_free_mgraph(cmgraph,objective->comm);
    cmgraph = fmgraph;
    maxmod = calc_total_modularity(clustering,cmgraph->graph);

    dl_stop_timer(&tmr);

    if (myid == 0) {
      if (objective->lbstats) {
        bbnd = 0;
        bnbr = 0;
        lb = loadbalance_stats_get(lb_offset+cmgraph->level,
            objective->lbstats);
        for (t=0;t<nthreads;++t) {
          bbnd += cmgraph->ucinfo[t]->bnd->size;
          bnbr += cmgraph->ucinfo[t]->nadj;
        }
        bbnd /= nthreads;
        bnbr /= nthreads;
        for (t=0;t<nthreads;++t) {
          dl_storemax(lb.bnd_imbalance,cmgraph->ucinfo[t]->bnd->size/bbnd);
          dl_storemax(lb.nbr_imbalance,cmgraph->ucinfo[t]->nadj/bnbr);
        }
        loadbalance_stats_replace(lb,lb_offset+cmgraph->level,
            objective->lbstats);
      }

      vprintf(objective->verbosity,NERSTRAND_VERBOSITY_HIGH,"Clustering "
          "modularity went from "PF_MOD_T" at level "PF_SIZE_T" to "PF_MOD_T
          " at level "PF_SIZE_T".\n",oldmod,cmgraph->level+1,maxmod,
          cmgraph->level);

      vprintf(objective->verbosity,NERSTRAND_VERBOSITY_HIGH,"Uncoarsening at "
          "level %zu took %0.05lfs\n",fmgraph->level,dl_poll_timer(&tmr));
    }
  }

  /* cleanup */
  par_free_mgraph(ograph,objective->comm);

  return clustering; 
}


static clustering_t * __par_cluster_graph_anyway(
    objective_t * const objective,
    const graph_t * const graph)
{
  int stop;
  real_t badj, bvtx,bbnd,bnbr;
  size_t lb_offset;
  tid_t t;
  mgraph_t * cmgraph, * fmgraph, * mgraph;
  clustering_t * clustering;
  mod_t maxmod, oldmod;
  loadbalance_stat_t lb;
  dl_timer_t tmr;

  const tid_t myid = dlthread_get_id(objective->comm);
  const tid_t nthreads = dlthread_get_nthreads(objective->comm);

  const real_t stopratio = objective->stopratio;
  mgraph_t * const ograph = par_setup_mgraph(0,graph,NULL,NULL,NULL, \
      objective->comm);

  if (objective->lbstats) {
    lb_offset = objective->lbstats->size;
  } else {
    lb_offset = 0;
  }

  dl_init_timer(&tmr);

  /* coarsen the graph */
  cmgraph = ograph;
  do { 
    if (myid == 0) {
      if (objective->lbstats) {
        badj = cmgraph->graph->nedges / (real_t)nthreads;
        bvtx = cmgraph->graph->nvtxs / (real_t)nthreads;
        lb.adj_imbalance = cmgraph->graph->mynedges[0] / badj;
        lb.vtx_imbalance = cmgraph->graph->mynvtxs[0] / badj;
        lb.nbr_imbalance = lb.bnd_imbalance = 0.0;
        for (t=1;t<nthreads;++t) {
          dl_storemax(lb.adj_imbalance,cmgraph->graph->mynedges[t]/badj);
          dl_storemax(lb.vtx_imbalance,cmgraph->graph->mynvtxs[t]/bvtx);
        }
        loadbalance_stats_push(lb,objective->lbstats);
      }
    }

    fmgraph = cmgraph;
    if (myid == 0) {
      vprintf(objective->verbosity,NERSTRAND_VERBOSITY_HIGH,"Coarsening graph "
          "of "PF_VTX_T" vertices and "PF_ADJ_T" edges\n",
          fmgraph->graph->nvtxs,fmgraph->graph->nedges);

      dl_reset_timer(&tmr);
      dl_start_timer(&tmr);
    }

    cmgraph = par_coarsen_graph(objective,fmgraph);

    if (cmgraph->graph->nedges > fmgraph->graph->nedges*objective->restep) {
      if (myid == 0) {
        vprintf(objective->verbosity,NERSTRAND_VERBOSITY_HIGH,"Making re-step "
            "do to small edge shrinkage %0.03lf < %0.03lf\n",
            ((double)cmgraph->graph->nedges)/fmgraph->graph->nedges,
            objective->restep);
      }

      mgraph = cmgraph;
      cmgraph = par_coarsen_graph(objective,mgraph);
      __par_squash_graphs(cmgraph,mgraph,fmgraph,objective->comm);
    }

    if (myid == 0) {
      dl_stop_timer(&tmr);
      vprintf(objective->verbosity,NERSTRAND_VERBOSITY_HIGH,"Coarsening at "
          "level %zu took %0.05lfs\n",fmgraph->level,dl_poll_timer(&tmr));
    }

    /* figure out if we should exit */
    stop = 0;
    switch (objective->stopcondition) {
      case NERSTRAND_STOPCONDITION_EDGES:
        stop = (cmgraph->graph->nedges >= fmgraph->graph->nedges*stopratio);
        break;
      case NERSTRAND_STOPCONDITION_VERTICES:
        stop = (cmgraph->graph->nvtxs >= fmgraph->graph->nvtxs*stopratio); 
        break;
      case NERSTRAND_STOPCONDITION_SIZE:
        stop = ((3.0*cmgraph->graph->nvtxs+2.0*cmgraph->graph->nedges) >= 
            (3.0*fmgraph->graph->nvtxs+2.0*fmgraph->graph->nedges)*stopratio); 
        break;
      default:
        dl_error("Unknown stop condition: %d\n",objective->stopcondition);
    }
    if (stop || cmgraph->level >= objective->maxlevel) {
      break;
    }
  } while (1);

  if (myid == 0) {
    vprintf(objective->verbosity,NERSTRAND_VERBOSITY_MEDIUM,"Initial "
        "clustering being made at level "PF_SIZE_T" on a graph with "PF_VTX_T
        " vertices and "PF_ADJ_T" edges.\n",cmgraph->level,
        cmgraph->graph->nvtxs,cmgraph->graph->nedges);
  }

  clustering = par_generate_aclustering(objective,cmgraph);
  maxmod = calc_total_modularity(clustering,cmgraph->graph);

  DL_ASSERT(check_clustering(clustering,cmgraph->graph) == 1,
      "Bad clustering generated\n");
  DL_ASSERT(check_ucinfo(cmgraph->graph,
        (const ucinfo_t * const *)cmgraph->ucinfo,
        (const cid_t * const *)clustering->where,clustering->nclusters) == 1, 
      "Bad ucinfo derived\n");
  DL_ASSERT(check_bnd(cmgraph,(const cid_t * const *)clustering->where) == 1,
      "Bad boundary after derive_ucinfo\n");

  /* project the clustering */
  while (cmgraph != ograph) {
    dl_reset_timer(&tmr);
    dl_start_timer(&tmr);

    oldmod = maxmod;
    fmgraph = cmgraph->finer;
    clustering = par_uncoarsen_graph(objective,clustering,fmgraph,cmgraph);
    par_free_mgraph(cmgraph,objective->comm);
    cmgraph = fmgraph;
    maxmod = calc_total_modularity(clustering,cmgraph->graph);

    dl_stop_timer(&tmr);

    if (myid == 0) {
      if (objective->lbstats) {
        bbnd = 0;
        bnbr = 0;
        lb = loadbalance_stats_get(lb_offset+cmgraph->level,
            objective->lbstats);
        for (t=0;t<nthreads;++t) {
          bbnd += cmgraph->ucinfo[t]->bnd->size;
          bnbr += cmgraph->ucinfo[t]->nadj;
        }
        bbnd /= nthreads;
        bnbr /= nthreads;
        for (t=0;t<nthreads;++t) {
          dl_storemax(lb.bnd_imbalance,cmgraph->ucinfo[t]->bnd->size/bbnd);
          dl_storemax(lb.nbr_imbalance,cmgraph->ucinfo[t]->nadj/bnbr);
        }
        loadbalance_stats_replace(lb,lb_offset+cmgraph->level,
            objective->lbstats);
      }
      vprintf(objective->verbosity,NERSTRAND_VERBOSITY_HIGH,"Clustering "
          "modularity went from "PF_MOD_T" at level "PF_SIZE_T" to "PF_MOD_T
          " at level "PF_SIZE_T".\n",oldmod,cmgraph->level+1,maxmod,
          cmgraph->level);

      vprintf(objective->verbosity,NERSTRAND_VERBOSITY_HIGH,"Uncoarsening at "
          "level %zu took %0.05lfs\n",fmgraph->level,dl_poll_timer(&tmr));
    }
  }

  /* cleanup */
  par_free_mgraph(ograph,objective->comm);

  return clustering; 
}




/******************************************************************************
* PUBLIC SERIAL FUNCTIONS *****************************************************
******************************************************************************/


clustering_t * cluster_graph(
    objective_t * const objective, 
    const graph_t * const graph)
{
  clustering_t * clustering;

  switch (objective->parttype) {
    case NERSTRAND_PARTITION_KWAY:
      clustering = __cluster_graph_kway(objective,graph);
      break;
    case NERSTRAND_PARTITION_ANYWAY:
      clustering = __cluster_graph_anyway(objective,graph);
      break;
    default :
      dl_error("Unknown partition type '%d'\n",objective->parttype);
  }

  return clustering; 
}




/******************************************************************************
* PUBLIC PARELLEL FUNCTIONS ***************************************************
******************************************************************************/


clustering_t * par_cluster_graph(
    objective_t * const objective, 
    const graph_t * const graph)
{
  clustering_t * clustering;

  switch (objective->parttype) {
    case NERSTRAND_PARTITION_KWAY:
      clustering = __par_cluster_graph_kway(objective,graph);
      break;
    case NERSTRAND_PARTITION_ANYWAY:
      clustering = __par_cluster_graph_anyway(objective,graph);
      break;
    default :
      dl_error("Unknown partition type '%d'\n",objective->parttype);
  }

  return clustering; 
}




#endif
