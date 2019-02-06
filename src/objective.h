/**
 * @file objective.h
 * @brief Types and function prototypes for manipulating and maintaining 
 * objectives
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-10-08
 */




#ifndef NERSTRAND_OBJECTIVE_H
#define NERSTRAND_OBJECTIVE_H




#include "base.h"
#include "stats.h"




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef struct timers_t {
  dl_timer_t total;
  dl_timer_t io;
  dl_timer_t clustering;
  dl_timer_t coarsening;
  dl_timer_t aggregation; 
  dl_timer_t contraction; 
  dl_timer_t sparsification;
  dl_timer_t edgeselection;
  dl_timer_t edgeremoval;
  dl_timer_t initcluster;
  dl_timer_t uncoarsening;
  dl_timer_t projection;
  dl_timer_t refinement;
} timers_t;


typedef struct objective_t {
  unsigned int seed;
  vtx_t cnvtxs;
  int ignore_weight;
  /* clustering parameters */
  vtx_t cnvtxs_per_cluster;
  vtx_t fm_vtx_limit;
  cid_t nclusters;
  double stopratio;
  double snratio;
  double restep;
  wgt_t degwgt;
  nerstrand_stopcondition_t stopcondition;
  nerstrand_aggtype_t aggtype;
  nerstrand_contype_t contype;
  nerstrand_sparsifytype_t sparsifytype;
  nerstrand_edgeremovaltype_t edgeremtype;
  nerstrand_ictype_t ictype;
  nerstrand_projtype_t projtype;
  nerstrand_reftype_t reftype;
  nerstrand_parttype_t parttype;
  nerstrand_verbosity_t verbosity;
  size_t nruns;
  size_t runnum;
  size_t nrefpass;
  size_t maxrefmoves;
  size_t initsolutions;
  size_t maxlevel;
  real_t max_agg_rate;
  real_t max_vwgt_ratio;
  vtx_t awayoffset;
  /* runtime parameters */
  int distribution;
  vtx_t blocksize;
  int time;
  timers_t timers;
  /* stats objects */
  mod_t * modstats;
  initialclustering_stats_t * icstats;
  loadbalance_stats_t * lbstats;
  aggregation_stats_t * aggstats;
  refinement_stats_t * refstats;
  /* threading */
  tid_t nthreads;
  dlthread_comm_t comm;
} objective_t;




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define create_objective __nerstrand_create_objective
/**
 * @brief Allocate and initialize an objective structure.
 *
 * @return The objective structure.
 */
objective_t * create_objective(void);


#define init_objective __nerstrand_init_objective
/**
 * @brief Initialize an objective structure
 *
 * @param objective The objective to initalize
 *
 * @return The initialized structure
 */
objective_t * init_objective(
    objective_t * objective);


#define setup_objective __nerstrand_setup_objective
/**
 * @brief Calculate derived fields of an objective structure 
 *
 * @param objective The objective to setup
 *
 * @return NERSTRAND_SUCCESS if successful, or an error code otherwise
 */
int setup_objective(
    objective_t * objective);


#define parse_objective __nerstrand_parse_objective
/**
 * @brief Create an objective struct from a list of options
 *
 * @param options An array indexed by options.
 * @param r_objective A reference to the objective pointer that will be
 * allocated.
 *
 * @return NERSTRAND_SUCCESS if successful, or an error if invalid values are
 * encountered.
 */
int parse_objective(
    const double * options, 
    objective_t ** r_objective);


#define free_objective __nerstrand_free_objective
/**
 * @brief Free an objective structure and associated memory
 *
 * @param objective The objective structure to free
 *
 * @return !0 on success
 */
int free_objective(
    objective_t * objective);


/* translation functions */
#define trans_parttype_string __nerstrand_trans_parttype_string
const char * trans_parttype_string(
    nerstrand_parttype_t type);


#define trans_aggtype_string __nerstrand_trans_aggtype_string
const char * trans_aggtype_string(
    nerstrand_aggtype_t type);


#define trans_contype_string __nerstrand_trans_contype_string
const char * trans_contype_string(
    nerstrand_contype_t type);


#define trans_spatype_string __nerstrand_trans_spatype_string
const char * trans_spatype_string(
    nerstrand_sparsifytype_t type);


#define trans_distype_string __nerstrand_trans_distype_string
const char * trans_distype_string(
    nerstrand_edgeremovaltype_t type);


#define trans_ictype_string __nerstrand_trans_ictype_string
const char * trans_ictype_string(
    nerstrand_ictype_t type);


#define trans_projtype_string __nerstrand_trans_projtype_string
const char * trans_projtype_string(
    nerstrand_projtype_t type);


#define trans_reftype_string __nerstrand_trans_reftype_string
const char * trans_reftype_string(
    nerstrand_reftype_t type);


#define trans_verbosity_string __nerstrand_trans_verbosity_string
const char * trans_verbosity_string(
    nerstrand_verbosity_t type);


#define trans_stopcondition_string __nerstrand_trans_stopcondition_string
const char * trans_stopcondition_string(
    nerstrand_stopcondition_t type);


#define trans_distribution_string __nerstrand_trans_distribution_string
const char * trans_distribution_string(
    nerstrand_distribution_t type);


#define trans_string_parttype __nerstrand_trans_string_parttype
nerstrand_parttype_t trans_string_parttype(
    const char * str);


#define trans_string_aggtype __nerstrand_trans_string_aggtype
nerstrand_aggtype_t trans_string_aggtype(
    const char * str);


#define trans_string_contype __nerstrand_trans_string_contype
nerstrand_contype_t trans_string_contype(
    const char * str);


#define trans_string_spatype __nerstrand_trans_string_spatype
nerstrand_sparsifytype_t trans_string_spatype(
    const char * str);


#define trans_string_distype __nerstrand_trans_string_distype
nerstrand_edgeremovaltype_t trans_string_distype(
    const char * str);


#define trans_string_ictype __nerstrand_trans_string_ictype
nerstrand_ictype_t trans_string_ictype(
    const char * str);


#define trans_string_projtype __nerstrand_trans_string_projtype
nerstrand_projtype_t trans_string_projtype(
    const char * str);


#define trans_string_reftype __nerstrand_trans_string_reftype
nerstrand_reftype_t trans_string_reftype(
    const char * str);


#define trans_string_verbosity __nerstrand_trans_string_verbosity
nerstrand_verbosity_t trans_string_verbosity(
    const char * str);


#define trans_string_stopcondition __nerstrand_trans_string_stopcondition
nerstrand_stopcondition_t trans_string_stopcondition(
    const char * str);


#define trans_string_distribution __nerstrand_trans_string_distribution
nerstrand_distribution_t trans_string_distribution(
    const char * str);




#endif
