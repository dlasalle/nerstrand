/**
 * @file io.c
 * @brief Functions for input and output 
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2014-03-03
 */




#ifndef NERSTRAND_IO_C
#define NERSTRAND_IO_C




#include "io.h"




/******************************************************************************
* MACROS **********************************************************************
******************************************************************************/


#define __SETSTATS(stat,stat_max,stat_avg,level,level_max) \
  do { \
    if ((stat) > (stat_max)) { \
      (stat_max) = (stat); \
      (level_max) = (level); \
    } \
    (stat_avg) += (stat); \
  } while (0)




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void print_clustering(
    clustering_t * const clustering, 
    graph_t const * const graph, 
    int const update, 
    int const verbosity)
{
  cid_t i, nclusters;
  vtx_t * cnvtxs;
  cluster_t * cluster;

  nclusters = clustering->nclusters;

  if (update) {
    calc_clustering(clustering,graph);
  }
  if (verbosity > NERSTRAND_VERBOSITY_MINIMUM) {
    printf("Clustering of "PF_CID_T" clusters:\n",nclusters);
    printf("\tTotal number of vertices "PF_VTX_T" and edge weight of "PF_TWGT_T
        "\n",graph->nvtxs,graph->tadjwgt);
    printf("\tTotal modularity of "PF_MOD_T"\n",
        calc_total_modularity(clustering,graph));

    /* calculate and print median cluster size */
    cnvtxs = vtx_alloc(nclusters);
    for (i=0;i<nclusters;++i) {
      cnvtxs[i] = clustering->clusters[i].nvtxs;
    }
    printf("\tCluster sizes median: "PF_VTX_T", mean: %0.01lf, min.: " \
        PF_VTX_T", max.: "PF_VTX_T", std. dev.: %0.01lf\n", \
        vtx_median(cnvtxs,nclusters),vtx_arithmetic_mean(cnvtxs,nclusters), \
        vtx_min_value(cnvtxs,nclusters),vtx_max_value(cnvtxs,nclusters), \
        vtx_stddev(cnvtxs,nclusters));
    dl_free(cnvtxs);

    if (verbosity > NERSTRAND_VERBOSITY_HIGH) {
      for (i=0;i<nclusters;++i) {
        cluster = clustering->clusters+i;
        printf("\tCluster "PF_CID_T" has "PF_VTX_T" vertices\n",i, \
            cluster->nvtxs);
        if (verbosity > NERSTRAND_VERBOSITY_MEDIUM) {
          printf("\t\ttotal number of vertices of "PF_VTX_T" and total edge " \
              "weight of "PF_TWGT_T"\n",cluster->nvtxs, \
              cluster->eadjwgt+cluster->iadjwgt);
          printf("\t\tinternal degree of "PF_TWGT_T" and external degree of " \
              PF_TWGT_T"\n",cluster->iadjwgt,cluster->eadjwgt);
        }
        printf("\t\tmodularity of "PF_MOD_T"\n", \
            calc_modularity(cluster,graph));
        printf("\n");
      }
    }
  }
}


void print_timers(
    objective_t const * const objective) 
{
  timers_t const * const timers = &(objective->timers);

  dl_print_header("TIMING",'$');
  printf("Total Time: %.05fs\n",dl_poll_timer(&(timers->total)));
  printf("\tIO: %.05fs\n",dl_poll_timer(&(timers->io)));
  printf("\tClustering: %.05fs\n",dl_poll_timer(&(timers->clustering)));
  printf("\t\tCoarsening: %.05fs\n",dl_poll_timer(&(timers->coarsening)));
  printf("\t\t\tAggregation: %.05fs\n",dl_poll_timer(&(timers->aggregation)));
  printf("\t\t\tContraction: %.05fs\n",dl_poll_timer(&(timers->contraction)));
  if (objective->sparsifytype != NERSTRAND_SPARSIFY_NONE) {
    printf("\t\t\tSparsification: %.05fs\n",
        dl_poll_timer(&(timers->sparsification)));
    printf("\t\t\t\tEdge Selection: %.05fs\n",
        dl_poll_timer(&(timers->edgeselection)));
    printf("\t\t\t\tEdge Removal: %.05fs\n",
        dl_poll_timer(&(timers->edgeremoval)));
  }
  printf("\t\tInitial Clustering: %.05fs\n",
      dl_poll_timer(&(timers->initcluster)));
  printf("\t\tUncoarsening: %.05fs\n",dl_poll_timer(&(timers->uncoarsening)));
  printf("\t\t\tProjection: %.05fs\n",dl_poll_timer(&(timers->projection)));
  printf("\t\t\tRefinement: %.05fs\n",dl_poll_timer(&(timers->refinement)));
  dl_print_footer('$');
}


void print_objective(
    objective_t const * const objective)
{
  dl_print_header("PARAMETERS",'%');
  printf("Runtime Parameters:\n");
  printf("Number of Threads: "PF_TID_T" | Verbosity: %s\n",objective->nthreads,
      trans_verbosity_string(objective->verbosity));
  if (objective->nthreads > 1) {
    printf("Distribution: %s", \
        trans_distribution_string(objective->distribution));
    if (objective->distribution == NERSTRAND_DISTRIBUTION_BLOCKCYCLIC) {
      printf(" | Blocksize: "PF_VTX_T,objective->blocksize);
    }
    printf("\n");
  }
  printf("Number of Runs: "PF_SIZE_T"\n",objective->nruns);
  printf("Clustering Parameters:\n");
  printf("Number of Clusters: "PF_CID_T" | Random Seed: %u\n",
      objective->nclusters,objective->seed);
  printf("Partition Type: %s\n",trans_parttype_string(objective->parttype)); 
  printf("Aggregation Type: %s | Contraction Type: %s\n",
      trans_aggtype_string(objective->aggtype),
      trans_contype_string(objective->contype));
  printf("Initial Clustering Type: %s\n",
      trans_ictype_string(objective->ictype));
  printf("Projection Type: %s | Refinement Type: %s\n",
      trans_projtype_string(objective->projtype),
      trans_reftype_string(objective->reftype));
  printf("Number of Initial Solutions: "PF_SIZE_T" | Number of Refinement "
      "Passes: "PF_SIZE_T"\n",objective->initsolutions,objective->nrefpass);
  printf("Stop Condition: %s | Stop Ratio: "PF_REAL_T"\n",
      trans_stopcondition_string(objective->stopcondition),
      objective->stopratio);
  if (objective->sparsifytype != NERSTRAND_SPARSIFY_NONE) {
    printf("Sparsification: %s | Edge Removal: %s\n",
        trans_spatype_string(objective->sparsifytype),
        trans_distype_string(objective->edgeremtype));
  }
  dl_print_footer('%');
}


void print_loadbalance_stats(
    loadbalance_stats_t * const stats)
{
  size_t nlvls;
  size_t adj_max_lvl, vtx_max_lvl, bnd_max_lvl, nbr_max_lvl;
  real_t adj_max, vtx_max, bnd_max, nbr_max, adj_avg, vtx_avg, bnd_avg, 
         nbr_avg; 
  loadbalance_stat_t stat;

  nlvls = stats->size;
  adj_max_lvl = vtx_max_lvl = bnd_max_lvl = nbr_max_lvl = 0;
  adj_max = vtx_max = bnd_max = nbr_max = 0.0;
  adj_avg = vtx_avg = bnd_avg = nbr_avg = 0.0;
  while (stats->size > 0) {
    stat = loadbalance_stats_pop(stats);
    __SETSTATS(stat.adj_imbalance,adj_max,adj_avg,stats->size,adj_max_lvl);
    __SETSTATS(stat.vtx_imbalance,vtx_max,vtx_avg,stats->size,vtx_max_lvl);
    __SETSTATS(stat.bnd_imbalance,bnd_max,bnd_avg,stats->size,bnd_max_lvl);
    __SETSTATS(stat.nbr_imbalance,nbr_max,nbr_avg,stats->size,nbr_max_lvl);
  }
  adj_avg /= nlvls;
  vtx_avg /= nlvls;
  bnd_avg /= nlvls;
  nbr_avg /= nlvls;

  dl_print_header("LOAD BALANCING STATISTICS",'=');
  printf("Edge Load Balancing: Avg = "PF_REAL_T", Max = "PF_REAL_T
      " (level = "PF_SIZE_T")\n",adj_avg,adj_max,adj_max_lvl);
  printf("Vertex Load Balancing: Avg = "PF_REAL_T", Max = "PF_REAL_T
      " (level = "PF_SIZE_T")\n",vtx_avg,vtx_max,vtx_max_lvl);
  printf("Boundary Load Balancing: Avg = "PF_REAL_T", Max = "PF_REAL_T
      " (level = "PF_SIZE_T")\n",bnd_avg,bnd_max,bnd_max_lvl);
  printf("Neighbor Load Balancing: Avg = "PF_REAL_T", Max = "PF_REAL_T
      " (level = "PF_SIZE_T")\n",nbr_avg,nbr_max,nbr_max_lvl);
  dl_print_footer('=');
}


void print_aggregation_stats(
    aggregation_stats_t * const stats)
{
  size_t nlvls;
  size_t nsn_max_lvl, nln_max_lvl, cr_max_lvl, mns_max_lvl, ans_max_lvl;
  vtx_t nsn_max, nln_max, mns_max, ans_max;
  real_t cr_max, nsn_avg, nln_avg, cr_avg, mns_avg, ans_avg; 
  aggregation_stat_t stat;

  nlvls = stats->size;
  nsn_max_lvl = nln_max_lvl = cr_max_lvl = mns_max_lvl = ans_max_lvl = 0;
  nsn_max = nln_max = mns_max = ans_max = 0;
  cr_max  = 0.0;
  nsn_avg = nln_avg = cr_avg = mns_avg = ans_avg = 0.0;
  while (stats->size > 0) {
    stat = aggregation_stats_pop(stats);
    __SETSTATS(stat.nsupernodes,nsn_max,nsn_avg,stats->size,nsn_max_lvl);
    __SETSTATS(stat.nlonenodes,nln_max,nln_avg,stats->size,nln_max_lvl);
    __SETSTATS(stat.coarsenrate,cr_max,cr_avg,stats->size,cr_max_lvl);
    __SETSTATS(stat.maxnodesize,mns_max,mns_avg,stats->size,mns_max_lvl);
    __SETSTATS(stat.avgnodesize,ans_max,ans_avg,stats->size,ans_max_lvl);
  }
  nsn_avg /= nlvls;
  nln_avg /= nlvls;
  cr_avg /= nlvls;
  mns_avg /= nlvls;
  ans_avg /= nlvls;

  dl_print_header("AGGREGATION STATISTICS",'=');
  printf("Number of Supernodes: Avg = "PF_REAL_T", Max = "PF_VTX_T
      " (level = "PF_SIZE_T")\n",nsn_avg,nsn_max,nsn_max_lvl);
  printf("Number of Lone-Nodes: Avg = "PF_REAL_T", Max = "PF_VTX_T
      " (level = "PF_SIZE_T")\n",nln_avg,nln_max,nln_max_lvl);
  printf("Coarsening Rate: Avg = "PF_REAL_T", Max = "PF_REAL_T
      " (level = "PF_SIZE_T")\n",cr_avg,cr_max,cr_max_lvl);
  printf("Maximum Node Size: Avg = "PF_REAL_T", Max = "PF_VTX_T
      " (level = "PF_SIZE_T")\n",mns_avg,mns_max,mns_max_lvl);
  printf("Average Node Size: Avg = "PF_REAL_T", Max = "PF_VTX_T
      " (level = "PF_SIZE_T")\n",ans_avg,ans_max,ans_max_lvl);
  dl_print_footer('=');
}


void print_refinement_stats(
    refinement_stats_t * const stats)
{
  size_t nlvls, np_max;
  size_t nm_max_lvl, um_max_lvl, nb_max_lvl, ag_max_lvl, pg_max_lvl, 
         np_max_lvl;
  vtx_t nm_max, um_max, nb_max;
  real_t ag_max, pg_max, nm_avg, um_avg, nb_avg, ag_avg, pg_avg, 
         np_avg; 
  refinement_stat_t stat;

  nlvls = stats->size;
  nm_max_lvl = um_max_lvl = nb_max_lvl = 0;
  ag_max_lvl = pg_max_lvl = np_max_lvl = 0;
  nm_max = um_max = nb_max = np_max = 0;
  ag_max = pg_max = 0.0;
  nb_avg = nm_avg = um_avg = ag_avg = pg_avg = np_avg = 0.0;
  while (stats->size > 0) {
    stat = refinement_stats_pop(stats);
    __SETSTATS(stat.nmoves,nm_max,nm_avg,stats->size,nm_max_lvl);
    __SETSTATS(stat.nunmoves,um_max,um_avg,stats->size,um_max_lvl);
    __SETSTATS(stat.nbnd,nb_max,nb_avg,stats->size,nb_max_lvl);
    __SETSTATS(stat.actgain,ag_max,ag_avg,stats->size,ag_max_lvl);
    __SETSTATS(stat.pergain,pg_max,pg_avg,stats->size,pg_max_lvl);
    __SETSTATS(stat.npasses,np_max,np_avg,stats->size,np_max_lvl);
  }

  /* average unmoves is the averaged percent of total moves, not levels */
  if (nm_avg+um_avg > 0) {
    um_avg *= 100.0 / (real_t)(nm_avg+um_avg);
  } else {
    um_avg = 0;
  }

  nm_avg /= nlvls;
  nb_avg /= nlvls;
  ag_avg /= nlvls;
  pg_avg /= nlvls;
  np_avg /= nlvls;

  dl_print_header("REFINEMENT STATISTICS",'=');
  printf("Number of Moves: Avg = "PF_REAL_T", Max = "PF_VTX_T
      " (level = "PF_SIZE_T")\n",nm_avg,nm_max,nm_max_lvl);
  printf("Undone Moves: Avg = "PF_REAL_T"%%, Max = "PF_VTX_T
      " (level = "PF_SIZE_T")\n",um_avg,um_max,um_max_lvl);
  printf("Actual Gain: Avg = "PF_REAL_T", Max = "PF_REAL_T
      " (level = "PF_SIZE_T")\n",ag_avg,ag_max,ag_max_lvl);
  printf("Percieved Gain: Avg = "PF_REAL_T", Max = "PF_REAL_T
      " (level = "PF_SIZE_T")\n",pg_avg,pg_max,pg_max_lvl);
  printf("Number of Passes: Avg = "PF_REAL_T", Max = "PF_SIZE_T
      " (level = "PF_SIZE_T")\n",np_avg,np_max,np_max_lvl);
  printf("Number of Boundary Vertices: Avg = "PF_REAL_T", Max = "PF_VTX_T
      " (level = "PF_SIZE_T")\n",nb_avg,nb_max,nb_max_lvl);
  dl_print_footer('=');
}


void print_modularity_stats(
    mod_t const * const mods, 
    size_t const nmods)
{
  dl_print_header("MODULARITY STATISTICS",'=');
  printf("Number of runs: "PF_SIZE_T"\n",nmods);
  printf("Arithmetic Mean Modularity: "PF_MOD_T"\n",
      mod_arithmetic_mean(mods,nmods));
  printf("Geometric Mean Modularity: "PF_MOD_T"\n",
      mod_geometric_mean(mods,nmods));
  printf("Std. Deviation Modularity: "PF_MOD_T"\n",
      mod_stddev(mods,nmods));
  printf("Max Modularity: "PF_MOD_T"\n",
      mod_max_value(mods,nmods));
  printf("Min Modularity: "PF_MOD_T"\n",
      mod_min_value(mods,nmods));
  dl_print_footer('=');
}


void print_initialclustering_stats(
    initialclustering_stats_t const * const stats)
{ 
  size_t const n = stats->nclusterings;
  size_t const m = stats->nruns;
  size_t const * const lvls = stats->lvls; 
  mod_t const * const mods = stats->mods; 
  vtx_t const * const nvtxs = stats->nvtxs;
  adj_t const * const nedges = stats->nedges;
  wgt_t const * const tadjwgt = stats->tadjwgt;
  cid_t const * const nclusters = stats->nclusters;
  cid_t const * const nbclusters = stats->nbclusters;

  dl_print_header("INITIAL CLUSTERING STATISTICS",'=');
  printf("Number of clusterings: "PF_SIZE_T"\n",n);
  printf("Number of clusters: min "PF_CID_T", max "PF_CID_T", mean %0.0lf\n",
      cid_min_value(nclusters,n),cid_max_value(nclusters,n),
      cid_arithmetic_mean(nclusters,n));
  printf("Number of clusters (chosen): min "PF_CID_T", max "PF_CID_T
      ", mean %0.0lf\n",cid_min_value(nbclusters,m),
      cid_max_value(nbclusters,m),cid_arithmetic_mean(nbclusters,m));
  printf("Levels: Min "PF_SIZE_T" Avg %0.01lf Max "PF_SIZE_T"\n",
      size_min_value(lvls,n),size_arithmetic_mean(lvls,n),
      size_max_value(lvls,n));
  printf("Mean # of Vertices: %0.01lf # of Edges: %0.01lf Edge Weight: "
      "%0.02lf\n",vtx_arithmetic_mean(nvtxs,n),adj_arithmetic_mean(nedges,n),
      wgt_arithmetic_mean(tadjwgt,n));
  printf("Maximum # of Vertices: "PF_VTX_T" # of Edges: "PF_ADJ_T
      " Edge Weight: "PF_WGT_T"\n",vtx_max_value(nvtxs,n),
      adj_max_value(nedges,n),wgt_max_value(tadjwgt,n)); 
  printf("Minimum # of Vertices: "PF_VTX_T" # of Edges: "PF_ADJ_T
      " Edge Weight: "PF_WGT_T"\n",vtx_min_value(nvtxs,n),
      adj_min_value(nedges,n),wgt_min_value(tadjwgt,n)); 
  printf("Arithmetic Mean Modularity: "PF_MOD_T"\n",
      mod_arithmetic_mean(mods,n));
  printf("Geometric Mean Modularity: "PF_MOD_T"\n",
      mod_geometric_mean(mods,n));
  printf("Std. Deviation Modularity: "PF_MOD_T"\n",
      mod_stddev(mods,n));
  printf("Max Modularity: "PF_MOD_T"\n",
      mod_max_value(mods,n));
  printf("Min Modularity: "PF_MOD_T"\n",
      mod_min_value(mods,n));
  dl_print_footer('=');
}




#endif

