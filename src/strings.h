/**
 * @file strings.h
 * @brief String constants
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2012, Regents of the University of Minnesota
 * @version 1
 * @date 2013-04-02
 */




#ifndef STRINGS_H
#define STRINGS_H




/******************************************************************************
* PRINT FORMAT STRINGS ********************************************************
******************************************************************************/


#define PF_REAL_T "%0.5f"
#define RF_REAL_T "%f"
#define PF_TWGT_T "%0.07Lf"
#define RF_TWGT_T "%Lf"
#define PF_UINT32_T "%"PRIu32
#define RF_UINT32_T "%"PRIu32
#define PF_UINT64_T "%"PRIu64
#define RF_UINT64_T "%"PRIu64
#define PF_FLOAT_T "%0.05f"
#define RF_FLOAT_T "%f"
#define PF_DOUBLE_T "%0.05lf"
#define RF_DOUBLE_T "%lf"


#ifdef NERSTRAND_64BIT_CLUSTERS
#define PF_CID_T PF_UINT64_T
#define RF_CID_T RF_UINT64_T
#else
#define PF_CID_T PF_UINT32_T
#define RF_CID_T RF_UINT32_T
#endif
#ifdef NERSTRAND_64BIT_VERTICES
#define PF_VTX_T PF_UINT64_T
#define RF_VTX_T RF_UINT64_T
#else
#define PF_VTX_T PF_UINT32_T
#define RF_VTX_T RF_UINT32_T
#endif
#ifdef NERSTRAND_64BIT_EDGES
#define PF_ADJ_T PF_UINT64_T
#define RF_ADJ_T RF_UINT64_T
#else
#define PF_ADJ_T PF_UINT32_T
#define RF_ADJ_T RF_UINT32_T
#endif
#ifdef NERSTRAND_64BIT_THREADS
#define PF_TID_T PF_UINT64_T
#define RF_TID_T RF_UINT64_T
#else
#define PF_TID_T PF_UINT32_T
#define RF_TID_T RF_UINT32_T
#endif
#ifdef NERSTRAND_DOUBLE_WEIGHTS
#define PF_WGT_T PF_DOUBLE_T
#define RF_WGT_T RF_DOUBLE_T
#else
#define PF_WGT_T PF_FLOAT_T
#define RF_WGT_T RF_FLOAT_T
#endif
#ifdef NERSTRAND_DOUBLE_MODULARITY
#define PF_MOD_T PF_DOUBLE_T
#define RF_MOD_T RF_DOUBLE_T
#else
#define PF_MOD_T PF_FLOAT_T
#define RF_MOD_T RF_FLOAT_T
#endif




/******************************************************************************
* OPTION STRINGS **************************************************************
******************************************************************************/


/* partition types */
#define NERSTRAND_STR_PARTITION_KWAY "kway"
#define NERSTRAND_STR_PARTITION_RB "rb"
#define NERSTRAND_STR_PARTITION_MRB "mrb"
#define NERSTRAND_STR_PARTITION_ANYWAY "anyway"


/* aggregation types */
#define NERSTRAND_STR_AGGREGATE_RM "rm"
#define NERSTRAND_STR_AGGREGATE_SHEM "shem"
#define NERSTRAND_STR_AGGREGATE_AGM "agm"
#define NERSTRAND_STR_AGGREGATE_AGH "agh"
#define NERSTRAND_STR_AGGREGATE_RC "rc"
#define NERSTRAND_STR_AGGREGATE_FC "fc"
#define NERSTRAND_STR_AGGREGATE_AGC "agc"
#define NERSTRAND_STR_AGGREGATE_AGW "agw"


/* contraction types */
#define NERSTRAND_STR_CONTRACT_SUM "sum"


/* sparsify types */
#define NERSTRAND_STR_SPARSIFY_NONE "none"
#define NERSTRAND_STR_SPARSIFY_RANDOM "random"
#define NERSTRAND_STR_SPARSIFY_LIGHT "light"
#define NERSTRAND_STR_SPARSIFY_HEAVY "heavy"
#define NERSTRAND_STR_SPARSIFY_DEGREE "degree"
#define NERSTRAND_STR_EDGEREMOVAL_DROP "drop"
#define NERSTRAND_STR_EDGEREMOVAL_LOOP "loop"
#define NERSTRAND_STR_EDGEREMOVAL_DISTRIBUTE "distribute"
#define NERSTRAND_STR_EDGEREMOVAL_PHANTOM "phantom"


/* initical clustering types */
#define NERSTRAND_STR_INITIAL_CLUSTERING_VTX "vtx"
#define NERSTRAND_STR_INITIAL_CLUSTERING_RVTX "rvtx"
#define NERSTRAND_STR_INITIAL_CLUSTERING_BFS "bfs"
#define NERSTRAND_STR_INITIAL_CLUSTERING_RANDOM "random"
#define NERSTRAND_STR_INITIAL_CLUSTERING_LP "lp"
#define NERSTRAND_STR_INITIAL_CLUSTERING_NEWMAN "newman"
#define NERSTRAND_STR_INITIAL_CLUSTERING_SEED "seed"
#define NERSTRAND_STR_INITIAL_CLUSTERING_GROW "grow"
#define NERSTRAND_STR_INITIAL_CLUSTERING_GROWKL "grow"


/* projection types */
#define NERSTRAND_STR_PROJECT_DIRECT "direct"
#define NERSTRAND_STR_PROJECT_SPARSE "sparse"


/* refinement types */
#define NERSTRAND_STR_REFINEMENT_EDGECUT "edgecut"
#define NERSTRAND_STR_REFINEMENT_GREEDY "greedy"
#define NERSTRAND_STR_REFINEMENT_RANDOM "random"
#define NERSTRAND_STR_REFINEMENT_KWAYMFM "kwaymfm"


/* verbosity strings */
#define NERSTRAND_STR_VERBOSITY_MINIMUM "minimum"
#define NERSTRAND_STR_VERBOSITY_LOW "low"
#define NERSTRAND_STR_VERBOSITY_MEDIUM "medium"
#define NERSTRAND_STR_VERBOSITY_HIGH "high"
#define NERSTRAND_STR_VERBOSITY_MAXIMUM "maximum"


/* stop condition strings */
#define NERSTRAND_STR_STOPCONDITION_VERTICES "vertices"
#define NERSTRAND_STR_STOPCONDITION_EDGES "edges"
#define NERSTRAND_STR_STOPCONDITION_SIZE "size"


/* distribution strings */
#define NERSTRAND_STR_DISTRIBUTION_BLOCK "block"
#define NERSTRAND_STR_DISTRIBUTION_CYCLIC "cyclic"
#define NERSTRAND_STR_DISTRIBUTION_BLOCKCYCLIC "blockcyclic"



#endif
