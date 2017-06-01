/**
 * @file objective.c
 * @brief Functions for manipulating and creating objectives
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-07-23
 */




#ifndef NERSTRAND_OBJECTIVE_C
#define NERSTRAND_OBJECTIVE_C




#include <omp.h>
#include "objective.h"




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLMEM_PREFIX objective
#define DLMEM_TYPE_T objective_t
#define DLMEM_DLTYPE DLTYPE_STRUCT
#define DLMEM_INITFUNCTION init_objective
#define DLMEM_STATIC
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_INITFUNCTION
#undef DLMEM_DLTYPE
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


static const cid_t DEFAULT_NCLUSTERS = 0;
static const double DEFAULT_STOP_RATIO = 0.95;
static const double DEFAULT_SUPERNODE_RATIO = 1.0;
static const nerstrand_parttype_t DEFAULT_PARTTYPE = \
    NERSTRAND_PARTITION_ANYWAY;
static const nerstrand_aggtype_t DEFAULT_KWAY_AGGTYPE = \
    NERSTRAND_AGGREGATE_AGC;
static const nerstrand_aggtype_t DEFAULT_ANYWAY_AGGTYPE = \
    NERSTRAND_AGGREGATE_AGC;
static const nerstrand_aggtype_t DEFAULT_AGGREGATE = NERSTRAND_AGGREGATE_AGC;
static const nerstrand_contype_t DEFAULT_CONTYPE = NERSTRAND_CONTRACT_SUM;
static const nerstrand_ictype_t DEFAULT_ICTYPE = \
    NERSTRAND_INITIAL_CLUSTERING_RVTX;
static const nerstrand_projtype_t DEFAULT_PROJTYPE = NERSTRAND_PROJECT_DIRECT; 
static const nerstrand_reftype_t DEFAULT_REFTYPE = \
    NERSTRAND_REFINEMENT_RANDOM;
static const size_t DEFAULT_NRUNS = 1; 
static const size_t DEFAULT_NREFPASS = 8;
static const size_t DEFAULT_MAXREFMOVES = 128;
static const size_t DEFAULT_INITSOLUTIONS = 16;
static const size_t DEFAULT_MAXLEVEL = 256;
static const vtx_t DEFAULT_AWAYOFFSET = 10000;
static const nerstrand_verbosity_t DEFAULT_VERBOSITY = \
   NERSTRAND_VERBOSITY_MINIMUM;
static const int DEFAULT_DO_TIMING = 0;
static const double DEFAULT_MAX_AGG_RATE = 0.5;
static const double DEFAULT_MAX_VWGT_RATIO = 0.1;
static const vtx_t DEFAULT_CNVTXS_PER_CLUSTER = 10; 
static const vtx_t DEFAULT_FM_VTX_LIMIT = 16;
static const vtx_t DEFAULT_BLOCKSIZE = 4096;
static const vtx_t DEFAULT_DISTRIBUTION = NERSTRAND_DISTRIBUTION_BLOCKCYCLIC;
static const size_t DEFAULT_RUNNUM = 0;
static const wgt_t DEFAULT_DEGREE_WEIGHT = 1.0;
static const nerstrand_stopcondition_t DEFAULT_STOPCONDITION = \
    NERSTRAND_STOPCONDITION_SIZE;
static const nerstrand_sparsifytype_t DEFAULT_SPATYPE = \
    NERSTRAND_SPARSIFY_NONE;
static const nerstrand_edgeremovaltype_t DEFAULT_EDGEREMTYPE = \
    NERSTRAND_EDGEREMOVAL_DROP;
static const vtx_t DEFAULT_CNVTXS = 0;
static const double DEFAULT_RESTEP = 1.0;


static const char * trans_table_part[] = {
  [NERSTRAND_PARTITION_KWAY] = NERSTRAND_STR_PARTITION_KWAY,
  [NERSTRAND_PARTITION_ANYWAY] = NERSTRAND_STR_PARTITION_ANYWAY
};


static const char * trans_table_agg[] = {
  [NERSTRAND_AGGREGATE_RM] = NERSTRAND_STR_AGGREGATE_RM,
  [NERSTRAND_AGGREGATE_SHEM] = NERSTRAND_STR_AGGREGATE_SHEM,
  [NERSTRAND_AGGREGATE_AGH] = NERSTRAND_STR_AGGREGATE_AGH,
  [NERSTRAND_AGGREGATE_RC] = NERSTRAND_STR_AGGREGATE_RC,
  [NERSTRAND_AGGREGATE_FC] = NERSTRAND_STR_AGGREGATE_FC,
  [NERSTRAND_AGGREGATE_AGM] = NERSTRAND_STR_AGGREGATE_AGM,
  [NERSTRAND_AGGREGATE_AGC] = NERSTRAND_STR_AGGREGATE_AGC
};


static const char * trans_table_spa[] = {
  [NERSTRAND_SPARSIFY_NONE] = NERSTRAND_STR_SPARSIFY_NONE,
  [NERSTRAND_SPARSIFY_RANDOM] = NERSTRAND_STR_SPARSIFY_RANDOM,
  [NERSTRAND_SPARSIFY_LIGHT] = NERSTRAND_STR_SPARSIFY_LIGHT,
  [NERSTRAND_SPARSIFY_HEAVY] = NERSTRAND_STR_SPARSIFY_HEAVY,
  [NERSTRAND_SPARSIFY_DEGREE] = NERSTRAND_STR_SPARSIFY_DEGREE
};


static const char * trans_table_dis[] = {
  [NERSTRAND_EDGEREMOVAL_DROP] = NERSTRAND_STR_EDGEREMOVAL_DROP,
  [NERSTRAND_EDGEREMOVAL_LOOP] = NERSTRAND_STR_EDGEREMOVAL_LOOP,
  [NERSTRAND_EDGEREMOVAL_DISTRIBUTE] = NERSTRAND_STR_EDGEREMOVAL_DISTRIBUTE,
  [NERSTRAND_EDGEREMOVAL_PHANTOM] = NERSTRAND_STR_EDGEREMOVAL_PHANTOM
};


static const char * trans_table_con[] = {
  [NERSTRAND_CONTRACT_SUM] = NERSTRAND_STR_CONTRACT_SUM
};


static const char * trans_table_ic[] = {
  [NERSTRAND_INITIAL_CLUSTERING_BFS] = NERSTRAND_STR_INITIAL_CLUSTERING_BFS,
  [NERSTRAND_INITIAL_CLUSTERING_VTX] = NERSTRAND_STR_INITIAL_CLUSTERING_VTX,
  [NERSTRAND_INITIAL_CLUSTERING_RVTX] = NERSTRAND_STR_INITIAL_CLUSTERING_RVTX,
  [NERSTRAND_INITIAL_CLUSTERING_RANDOM] = 
      NERSTRAND_STR_INITIAL_CLUSTERING_RANDOM,
  [NERSTRAND_INITIAL_CLUSTERING_SEED] = NERSTRAND_STR_INITIAL_CLUSTERING_SEED,
  [NERSTRAND_INITIAL_CLUSTERING_LP] = NERSTRAND_STR_INITIAL_CLUSTERING_LP,
  [NERSTRAND_INITIAL_CLUSTERING_NEWMAN] = 
      NERSTRAND_STR_INITIAL_CLUSTERING_NEWMAN,
  [NERSTRAND_INITIAL_CLUSTERING_GROW] = NERSTRAND_STR_INITIAL_CLUSTERING_GROW,
  [NERSTRAND_INITIAL_CLUSTERING_GROWKL] = 
      NERSTRAND_STR_INITIAL_CLUSTERING_GROWKL
};


static const char * trans_table_proj[] = {
  [NERSTRAND_PROJECT_DIRECT] = NERSTRAND_STR_PROJECT_DIRECT,
  [NERSTRAND_PROJECT_SPARSE] = NERSTRAND_STR_PROJECT_SPARSE
};


static const char * trans_table_ref[] = {
  [NERSTRAND_REFINEMENT_GREEDY] = NERSTRAND_STR_REFINEMENT_GREEDY,
  [NERSTRAND_REFINEMENT_RANDOM] = NERSTRAND_STR_REFINEMENT_RANDOM,
};


static const char * trans_table_verbosity[] = {
  [NERSTRAND_VERBOSITY_MINIMUM] = NERSTRAND_STR_VERBOSITY_MINIMUM,
  [NERSTRAND_VERBOSITY_LOW] = NERSTRAND_STR_VERBOSITY_LOW,
  [NERSTRAND_VERBOSITY_MEDIUM] = NERSTRAND_STR_VERBOSITY_MEDIUM,
  [NERSTRAND_VERBOSITY_HIGH] = NERSTRAND_STR_VERBOSITY_HIGH,
  [NERSTRAND_VERBOSITY_MAXIMUM] = NERSTRAND_STR_VERBOSITY_MAXIMUM,
};


static const char * trans_table_stopcondition[] = {
  [NERSTRAND_STOPCONDITION_VERTICES] = NERSTRAND_STR_STOPCONDITION_VERTICES,
  [NERSTRAND_STOPCONDITION_EDGES] = NERSTRAND_STR_STOPCONDITION_EDGES,
  [NERSTRAND_STOPCONDITION_SIZE] = NERSTRAND_STR_STOPCONDITION_SIZE
};


static const char * trans_table_distribution[] = {
  [NERSTRAND_DISTRIBUTION_BLOCK] = NERSTRAND_STR_DISTRIBUTION_BLOCK,
  [NERSTRAND_DISTRIBUTION_CYCLIC] = NERSTRAND_STR_DISTRIBUTION_CYCLIC,
  [NERSTRAND_DISTRIBUTION_BLOCKCYCLIC] = NERSTRAND_STR_DISTRIBUTION_BLOCKCYCLIC
};




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


#define __find_string(table,str) \
  __find_string_in(table,sizeof(table)/sizeof(*table),str)

static inline int __find_string_in(
    const char * const table[], 
    const int n,
    const char * const str)
{
  int i;
  if (str == NULL) {
    return -1;
  }
  for (i=0;i<n;++i) {
    if (strcmp(table[i],str) == 0) {
      break;
    }
  }
  if (i == n) {
    return -1;
  } else {
    return i;
  }
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


objective_t * create_objective(void)
{
  objective_t * objective;
  objective = objective_alloc(1);
  return init_objective(objective);
}


objective_t * init_objective(
    objective_t * const objective)
{
  objective->seed = (unsigned int)time(NULL);

  // TODO: remove this dependency on openmp
  if (omp_in_parallel()) {
    /* if we're called from inside of a parallel region, run serially */
    objective->nthreads = 1;
  } else {
    /* otherwise use as many threads as possible */
    objective->nthreads = omp_get_max_threads();
  }

  objective->nclusters = DEFAULT_NCLUSTERS;
  objective->parttype = DEFAULT_PARTTYPE;
  objective->stopratio = DEFAULT_STOP_RATIO;
  objective->snratio = DEFAULT_SUPERNODE_RATIO;
  objective->aggtype = DEFAULT_AGGREGATE;
  objective->contype = DEFAULT_CONTYPE;
  objective->sparsifytype = DEFAULT_SPATYPE;
  objective->edgeremtype = DEFAULT_EDGEREMTYPE;
  objective->ictype = DEFAULT_ICTYPE;
  objective->projtype = DEFAULT_PROJTYPE;
  objective->reftype = DEFAULT_REFTYPE;
  objective->nruns = DEFAULT_NRUNS;
  objective->runnum = DEFAULT_RUNNUM;
  objective->nrefpass = DEFAULT_NREFPASS;
  objective->maxrefmoves = DEFAULT_MAXREFMOVES;
  objective->verbosity = DEFAULT_VERBOSITY;
  objective->initsolutions = DEFAULT_INITSOLUTIONS;
  objective->max_agg_rate = DEFAULT_MAX_AGG_RATE;
  objective->max_vwgt_ratio = DEFAULT_MAX_VWGT_RATIO;
  objective->cnvtxs_per_cluster = DEFAULT_CNVTXS_PER_CLUSTER;
  objective->fm_vtx_limit = DEFAULT_FM_VTX_LIMIT;
  objective->time = DEFAULT_DO_TIMING;
  objective->maxlevel = DEFAULT_MAXLEVEL;
  objective->degwgt = DEFAULT_DEGREE_WEIGHT;
  objective->stopcondition = DEFAULT_STOPCONDITION;
  objective->awayoffset = DEFAULT_AWAYOFFSET;
  objective->cnvtxs = DEFAULT_CNVTXS;
  objective->blocksize = DEFAULT_BLOCKSIZE;
  objective->distribution = DEFAULT_DISTRIBUTION;
  objective->restep = DEFAULT_RESTEP;
  objective->modstats = NULL;
  objective->icstats = NULL;
  objective->lbstats = NULL;
  objective->aggstats = NULL;
  objective->refstats = NULL;

  /* timers */
  dl_init_timer(&(objective->timers.total));
  dl_init_timer(&(objective->timers.io));
  dl_init_timer(&(objective->timers.clustering));
  dl_init_timer(&(objective->timers.coarsening));
  dl_init_timer(&(objective->timers.aggregation));
  dl_init_timer(&(objective->timers.contraction));
  dl_init_timer(&(objective->timers.sparsification));
  dl_init_timer(&(objective->timers.edgeselection));
  dl_init_timer(&(objective->timers.edgeremoval));
  dl_init_timer(&(objective->timers.initcluster));
  dl_init_timer(&(objective->timers.uncoarsening));
  dl_init_timer(&(objective->timers.projection));
  dl_init_timer(&(objective->timers.refinement));

  /* derived fields */
  setup_objective(objective);

  return objective;
}


int setup_objective(
    objective_t * objective)
{
  int rv = NERSTRAND_SUCCESS;

  /* set parttype */
  if (objective->nclusters != 0) {
    objective->ictype = NERSTRAND_INITIAL_CLUSTERING_GROW;
    objective->parttype = NERSTRAND_PARTITION_KWAY;
    objective->cnvtxs = objective->nclusters*objective->cnvtxs_per_cluster;
  } else { 
    objective->parttype = NERSTRAND_PARTITION_ANYWAY;
  }

  /* fix initial solutions */
  if (objective->ictype == NERSTRAND_INITIAL_CLUSTERING_NEWMAN ||
      objective->ictype == NERSTRAND_INITIAL_CLUSTERING_VTX) {
    objective->initsolutions = 1;
  }
  
  /* set projection type */
  if (objective->sparsifytype != NERSTRAND_SPARSIFY_NONE) {
    objective->projtype = NERSTRAND_PROJECT_SPARSE;
  }

  /* allocate derived things */
  if (objective->modstats == YES_POINTER) {
    objective->modstats = mod_alloc(objective->nruns);
  }
  if (objective->icstats == YES_POINTER) {
    objective->icstats = initialclustering_stats_create(objective->nruns,
        objective->initsolutions);
  }

  /* start checks */
  if (objective->nruns < 1) {
    eprintf("Invalid number of runs specificied: %zu\n",objective->nruns);
    rv = NERSTRAND_ERROR_INVALIDOPTIONS;
    goto CLEANUP;
  }
  if (objective->parttype == NERSTRAND_PARTITION_KWAY) {
    /* insufficient number of cnvtxs */
    if (objective->cnvtxs < objective->nclusters) {
      eprintf("Invalid number of coarse vertices: "PF_VTX_T"\n", \
          objective->cnvtxs);
      rv = NERSTRAND_ERROR_INVALIDOPTIONS;
      goto CLEANUP;
    }
  } else if (objective->parttype == NERSTRAND_PARTITION_ANYWAY) {
    /* no anyway specific checks yet */
  }

  /* check distribution */
  if (objective->distribution == NERSTRAND_DISTRIBUTION_BLOCKCYCLIC && \
      objective->blocksize == 0) {
    eprintf("Invalid block size for block-cyclic distribution: "PF_VTX_T"\n", \
        objective->blocksize);
    rv = NERSTRAND_ERROR_INVALIDOPTIONS;
    goto CLEANUP;
  }

  if (objective->stopratio > 1.0) {
    eprintf("Invalid stop ratio: %lf\n",(double)objective->stopratio);
    rv = NERSTRAND_ERROR_INVALIDOPTIONS; 
    goto CLEANUP;
  }
  if (objective->initsolutions == 0) {
    eprintf("Invalid number of initial solutions: %zu\n", \
        objective->initsolutions);
    rv = NERSTRAND_ERROR_INVALIDOPTIONS;
    goto CLEANUP;
  }
  if (objective->nthreads == 0) {
    eprintf("Invalid number of threads: "PF_TID_T"\n",objective->nthreads);
    rv = NERSTRAND_ERROR_INVALIDOPTIONS;
    goto CLEANUP;
  }
  if (objective->max_agg_rate == 1.0) {
    eprintf("Invalid aggregation rate: %lf\n",(double)objective->max_agg_rate);
    rv = NERSTRAND_ERROR_INVALIDOPTIONS;
    goto CLEANUP;
  }

  return rv;

  CLEANUP:
  return rv;
}


int parse_objective(
    const double * const options, 
    objective_t ** const r_objective)
{
  int rv = NERSTRAND_SUCCESS;
  objective_t * objective = NULL;

  objective = objective_calloc(1);
  if (objective == NULL) {
    rv = NERSTRAND_ERROR_NOTENOUGHMEMORY;
    goto CLEANUP;
  }

  if (options[NERSTRAND_OPTION_SEED] != NERSTRAND_VAL_OFF) {
    objective->seed = (unsigned int)options[NERSTRAND_OPTION_SEED];
  }
  if (options[NERSTRAND_OPTION_NCLUSTERS] != NERSTRAND_VAL_OFF) {
    objective->nclusters = (cid_t)options[NERSTRAND_OPTION_NCLUSTERS];
    /* set agg type if they didn't specify */
    if (options[NERSTRAND_OPTION_AGGTYPE] == NERSTRAND_VAL_OFF) {
      if (objective->nclusters != 0) {
        objective->aggtype = DEFAULT_KWAY_AGGTYPE;
      } else {
        objective->aggtype = DEFAULT_ANYWAY_AGGTYPE;
      }
    }
  }
  if (options[NERSTRAND_OPTION_NRUNS] != NERSTRAND_VAL_OFF) {
    objective->nruns = (size_t)options[NERSTRAND_OPTION_NRUNS];
  }
  if (options[NERSTRAND_OPTION_NINITSOLUTIONS] != NERSTRAND_VAL_OFF) {
    objective->initsolutions = 
        (size_t)options[NERSTRAND_OPTION_NINITSOLUTIONS];
  }
  if (options[NERSTRAND_OPTION_AGGTYPE] != NERSTRAND_VAL_OFF) {
    objective->aggtype = 
        (nerstrand_aggtype_t)options[NERSTRAND_OPTION_AGGTYPE];
  }
  if (options[NERSTRAND_OPTION_CONTYPE] != NERSTRAND_VAL_OFF) {
    objective->contype = 
        (nerstrand_contype_t)options[NERSTRAND_OPTION_CONTYPE];
  }
  if (options[NERSTRAND_OPTION_SPATYPE] != NERSTRAND_VAL_OFF) {
    objective->sparsifytype = 
        (nerstrand_sparsifytype_t)options[NERSTRAND_OPTION_SPATYPE];
  }
  if (options[NERSTRAND_OPTION_DISTYPE] != NERSTRAND_VAL_OFF) {
    objective->edgeremtype = 
        (nerstrand_edgeremovaltype_t)options[NERSTRAND_OPTION_DISTYPE];
  }
  if (options[NERSTRAND_OPTION_REFTYPE] != NERSTRAND_VAL_OFF) {
    objective->reftype = 
        (nerstrand_reftype_t)options[NERSTRAND_OPTION_REFTYPE];
  }
  if (options[NERSTRAND_OPTION_INITYPE] != NERSTRAND_VAL_OFF) {
    objective->ictype = (nerstrand_ictype_t)options[NERSTRAND_OPTION_INITYPE];
  }
  if (options[NERSTRAND_OPTION_PARTYPE] != NERSTRAND_VAL_OFF) {
    objective->parttype = 
        (nerstrand_parttype_t)options[NERSTRAND_OPTION_PARTYPE];
  }
  if (options[NERSTRAND_OPTION_DEGREE_WEIGHT] != NERSTRAND_VAL_OFF) {
    objective->degwgt = (wgt_t)options[NERSTRAND_OPTION_DEGREE_WEIGHT];
  }
  if (options[NERSTRAND_OPTION_NREFPASS] != NERSTRAND_VAL_OFF) {
    objective->nrefpass = (size_t)options[NERSTRAND_OPTION_NREFPASS];
  }
  if (options[NERSTRAND_OPTION_MAXREFMOVES] != NERSTRAND_VAL_OFF) {
    objective->maxrefmoves = (size_t)options[NERSTRAND_OPTION_MAXREFMOVES];
  }
  if (options[NERSTRAND_OPTION_VERBOSITY] != NERSTRAND_VAL_OFF) {
    objective->verbosity = 
        (nerstrand_verbosity_t)options[NERSTRAND_OPTION_VERBOSITY];
  }
  if (options[NERSTRAND_OPTION_MODSTATS] != NERSTRAND_VAL_OFF) {
    objective->modstats = YES_POINTER;
  }
  if (options[NERSTRAND_OPTION_ICSTATS] != NERSTRAND_VAL_OFF) {
    objective->icstats = YES_POINTER;
  }
  if (options[NERSTRAND_OPTION_LBSTATS] != NERSTRAND_VAL_OFF) {
    objective->lbstats = loadbalance_stats_create();
  }
  if (options[NERSTRAND_OPTION_AGGSTATS] != NERSTRAND_VAL_OFF) {
    objective->aggstats = aggregation_stats_create();
  }
  if (options[NERSTRAND_OPTION_REFSTATS] != NERSTRAND_VAL_OFF) {
    objective->refstats = refinement_stats_create();
  }
  if (options[NERSTRAND_OPTION_TIME] != NERSTRAND_VAL_OFF) {
    objective->time = 1;
  }
  if (options[NERSTRAND_OPTION_AGG_RATE] != NERSTRAND_VAL_OFF) {
    objective->max_agg_rate = (double)options[NERSTRAND_OPTION_AGG_RATE];
  }
  if (options[NERSTRAND_OPTION_SUPERNODE_RATIO] != NERSTRAND_VAL_OFF) {
    objective->snratio = (double)options[NERSTRAND_OPTION_SUPERNODE_RATIO];
  }
  if (options[NERSTRAND_OPTION_CNVTXS_PER_CLUSTER] != NERSTRAND_VAL_OFF) {
    objective->cnvtxs_per_cluster = \
        (vtx_t)options[NERSTRAND_OPTION_CNVTXS_PER_CLUSTER];
  }
  if (options[NERSTRAND_OPTION_NTHREADS] != NERSTRAND_VAL_OFF) {
    objective->nthreads = (tid_t)options[NERSTRAND_OPTION_NTHREADS];
  }
  if (options[NERSTRAND_OPTION_STOPRATIO] != NERSTRAND_VAL_OFF) {
    objective->stopratio = (double)options[NERSTRAND_OPTION_STOPRATIO];
  }
  if (options[NERSTRAND_OPTION_STOPCONDITION] != NERSTRAND_VAL_OFF) {
    objective->stopcondition = \
        (nerstrand_stopcondition_t)options[NERSTRAND_OPTION_STOPCONDITION];
  }
  if (options[NERSTRAND_OPTION_BLOCKSIZE] != NERSTRAND_VAL_OFF) {
    objective->blocksize = (vtx_t)options[NERSTRAND_OPTION_BLOCKSIZE];
  }
  if (options[NERSTRAND_OPTION_DISTRIBUTION] != NERSTRAND_VAL_OFF) {
    objective->distribution = \
      (nerstrand_distribution_t)options[NERSTRAND_OPTION_DISTRIBUTION];
  }
  if (options[NERSTRAND_OPTION_RESTEP] != NERSTRAND_VAL_OFF) {
    objective->restep = (double)options[NERSTRAND_OPTION_RESTEP];
  }

  if ((rv = setup_objective(objective)) != NERSTRAND_SUCCESS) {
    goto CLEANUP;
  }

  *r_objective = objective;

  return rv; 

  CLEANUP:

  if (objective) {
    free_objective(objective);
    objective = NULL;
  }
  *r_objective = NULL;

  return rv;
}


int free_objective(
    objective_t * objective)
{
  if (objective->modstats) {
    dl_free(objective->modstats);
  }
  if (objective->icstats) {
    initialclustering_stats_free(objective->icstats);
  }
  if (objective->refstats) {
    refinement_stats_free(objective->refstats);
  }
  if (objective->aggstats) {
    aggregation_stats_free(objective->aggstats);
  }
  if (objective->lbstats) {
    loadbalance_stats_free(objective->lbstats);
  }
  dl_free(objective);

  return NERSTRAND_SUCCESS;
}


const char * trans_parttype_string(
    const nerstrand_parttype_t type)
{
  return trans_table_part[type];
}


const char * trans_aggtype_string(
    const nerstrand_aggtype_t type)
{
  return trans_table_agg[type];
}


const char * trans_contype_string(
    const nerstrand_contype_t type)
{
  return trans_table_con[type];
}


const char * trans_spatype_string(
    const nerstrand_sparsifytype_t type)
{
  return trans_table_spa[type];
}


const char * trans_distype_string(
    const nerstrand_edgeremovaltype_t type)
{
  return trans_table_dis[type];
}


const char * trans_ictype_string(
    const nerstrand_ictype_t type)
{
  return trans_table_ic[type];
}


const char * trans_projtype_string(
    const nerstrand_projtype_t type)
{
  return trans_table_proj[type];
}


const char * trans_reftype_string(
    const nerstrand_reftype_t type)
{
  return trans_table_ref[type];
}


const char * trans_verbosity_string(
    const nerstrand_verbosity_t type)
{
  return trans_table_verbosity[type];
}


const char * trans_stopcondition_string(
    const nerstrand_stopcondition_t type)
{
  return trans_table_stopcondition[type];
}


const char * trans_distribution_string(
    const nerstrand_distribution_t type)
{
  return trans_table_distribution[type];
}


nerstrand_parttype_t trans_string_parttype(
    const char * const str)
{
  int i = __find_string(trans_table_part,str);
  if (i < 0) {
    dl_error("Unknown Partition Type '%s'\n",str);
  } else {
    return (nerstrand_parttype_t)i;
  }
}


nerstrand_aggtype_t trans_string_aggtype(
    const char * const str)
{
  int i = __find_string(trans_table_agg,str);
  if (i < 0) {
    dl_error("Unknown Aggregation Type '%s'\n",str);
  } else {
    return (nerstrand_aggtype_t)i;
  }
}


nerstrand_contype_t trans_string_contype(
    const char * const str)
{
  int i = __find_string(trans_table_con,str);
  if (i < 0) {
    dl_error("Unknown Contraction Type '%s'\n",str);
  } else {
    return (nerstrand_contype_t)i;
  }
}


nerstrand_sparsifytype_t trans_string_spatype(
    const char * const str)
{
  int i = __find_string(trans_table_spa,str);
  if (i < 0) {
    dl_error("Unknown Sparsify Type '%s'\n",str);
  } else {
    return (nerstrand_sparsifytype_t)i;
  }
}


nerstrand_edgeremovaltype_t trans_string_distype(
    const char * const str)
{
  int i = __find_string(trans_table_dis,str);
  if (i < 0) {
    dl_error("Unknown Edge Removal Type '%s'\n",str);
  } else {
    return (nerstrand_edgeremovaltype_t)i;
  }
}


nerstrand_ictype_t trans_string_ictype(
    const char * const str)
{
  int i = __find_string(trans_table_ic,str);
  if (i < 0) {
    dl_error("Unknown Initial Clustering Type '%s'\n",str);
  } else {
    return (nerstrand_ictype_t)i;
  }
}


nerstrand_projtype_t trans_string_projtype(
    const char * const str)
{
  int i = __find_string(trans_table_proj,str);
  if (i < 0) {
    dl_error("Unknown Projection Type '%s'\n",str);
  } else {
    return (nerstrand_projtype_t)i;
  }
}


nerstrand_reftype_t trans_string_reftype(
    const char * const str)
{
  int i = __find_string(trans_table_ref,str);
  if (i < 0) {
    dl_error("Unknown Refinement Type '%s'\n",str);
  } else {
    return (nerstrand_reftype_t)i;
  }
}


nerstrand_verbosity_t trans_string_verbosity(
    const char * const str)
{
  int i = __find_string(trans_table_verbosity,str);
  if (i < 0) {
    dl_error("Unknown Verbosity Level '%s'\n",str);
  } else {
    return (nerstrand_verbosity_t)i;
  }
}


nerstrand_stopcondition_t trans_string_stopcondition(
    const char * const str)
{
  int i = __find_string(trans_table_stopcondition,str);
  if (i < 0) {
    dl_error("Unknown Stop Condition '%s'\n",str);
  } else {
    return (nerstrand_stopcondition_t)i;
  }
}


nerstrand_distribution_t trans_string_distribution(
    const char * const str)
{
  int i = __find_string(trans_table_distribution,str);
  if (i < 0) {
    dl_error("Unknown Distribution '%s'\n",str);
  } else {
    return (nerstrand_distribution_t)i;
  }
}


#endif
