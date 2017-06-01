/**
 * @file ucinfo.h
 * @brief Types and function prototypes for uncoarsening bookkeeping
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-10-08
 */



#ifndef NERSTRAND_UCINFO_H
#define NERSTRAND_UCINFO_H




#include "base.h"
#include "graph.h"




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef struct nbrinfo_t {
  wgt_t iw;
  vtx_t in;
  wgt_t ew;
  cid_t nnbrs;
  adj_t nbrstart;
} nbrinfo_t;


typedef struct ucinfo_t {
  tid_t npar;
  vtx_iset_t * bnd;
  nbrinfo_t * nbrinfo;
  adj_t nadj;
  adj_t maxnadj;
  cid_t * adjncy;
  wgt_t * adjwgt;
  vtx_t * adjnum;
  int free_adjncy;
  int free_adjwgt;
  int free_adjnum;
} ucinfo_t;




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


static adj_t const NO_NBRS = ((adj_t)-1);




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLMEM_PREFIX nbrinfo
#define DLMEM_TYPE_T nbrinfo_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMEM_PREFIX ucinfo
#define DLMEM_TYPE_T ucinfo_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLPQ_PREFIX rv
#define DLPQ_KEY_T real_t
#define DLPQ_VAL_T vtx_t
#include "dlpq_headers.h"
#undef DLPQ_PREFIX
#undef DLPQ_KEY_T
#undef DLPQ_VAL_T




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define init_nbrinfo __nerstrand_init_nbrinfo
/**
 * @brief Initialize a nbrinfo structure
 *
 * @param nbrinfo A pointer to the structure to initialize
 *
 * @return The initialized structure 
 */
nbrinfo_t * init_nbrinfo(
    nbrinfo_t * nbrinfo);


#define free_nbrinfo __nerstrand_free_nbrinfo
/**
 * @brief Free a nbrinfo structure and associated memory
 *
 * @param nbrinfo The nbrinfo to free 
 *
 * @return !0 on success
 */
int free_nbrinfo(
    nbrinfo_t * nbrinfo);


#define init_ucinfo __nerstrand_init_ucinfo
/**
 * @brief Initialize a ucinfo structure
 *
 * @param ucinfo A pointer to the structure to initialize
 *
 * @return The initialized structure
 */
ucinfo_t * init_ucinfo(
    ucinfo_t * ucinfo);


#define setup_ucinfo __nerstrand_setup_ucinfo
/**
 * @brief Setup a ucinfo structure given the supplied parameters
 *
 * @param nvtxs Vertex counts for each part of the graph 
 * @param nbrinfo The nbrinfo arrays (can be NULL)
 * @param adjncy The adjancent cluster arrays (can be NULL)
 * @param adjwgt The cluster connection weight arrays (can be NULL)
 * @param adjnum The cluster connection edge number arrays (can be NULL)
 *
 * @return The setup ucinfo
 */
ucinfo_t * setup_ucinfo(
    vtx_t nvtxs, 
    nbrinfo_t * nbrinfo, 
    cid_t * adjncy, 
    wgt_t * adjwgt,
    vtx_t * adjnum);


#define clone_ucinfo __nerstrand_clone_ucinfo
/**
 * @brief Create a copy of a uc info
 *
 * @param nvtxs The number of vertices contained within the ucinfo
 * @param ucinfo The ucinfo to copy
 *
 * @return The newly duplicated ucinfo
 */
ucinfo_t * clone_ucinfo(
    vtx_t nvtxs, 
    ucinfo_t const * ucinfo);


#define clone_ucinfos __nerstrand_clone_ucinfos
/**
 * @brief Create a copy of a set of ucinfos
 *
 * @param nvtxs The number of vertices for each ucinfo
 * @param ucinfos The ucinfos to copy
 * @param nthreads The number of ucinfos/threads
 *
 * @return  The newly allocated ucinfos
 */
ucinfo_t ** clone_ucinfos(
    vtx_t const * nvtxs, 
    ucinfo_t const * const * ucinfos,
    tid_t nthreads);


#define free_ucinfo __nerstrand_free_ucinfo
/**
 * @brief Free a ucinfo structure and associate memory
 *
 * @param ucinfo The structure to free
 *
 * @return !0 on success
 */
int free_ucinfo(
    ucinfo_t * ucinfo);


#define free_ucinfos __nerstrand_free_ucinfos
/**
 * @brief Free ucinfo structures and associate memory
 *
 * @param ucinfo The structures to free
 * @param nthreads The number of structures/threads to free
 *
 * @return !0 on success
 */
int free_ucinfos(
    ucinfo_t ** ucinfo, 
    tid_t  nthreads);


#define next_nbr_adj __nerstrand_next_nbr_adj
/**
 * @brief Acquire the offset for the next available chunk nbr adj space 
 *
 * @param ucinfo The ucinfo to get the offset from
 * @param size The amount of space required
 *
 * @return The index of the allocated space
 */
adj_t next_nbr_adj(
    ucinfo_t * ucinfo, 
    adj_t size);


#define release_nbr_adj __nerstrand_release_nbr_adj
/**
 * @brief Release the most recently allocated nbr adj space 
 *
 * @param ucinfo The ucinfo to release the space from
 * @param nbrinfo The nbrinfo who's space is to be released
 * @param size The size of the space
 *
 * @return !0 on success 
 */
int release_nbr_adj(
    ucinfo_t * ucinfo, 
    nbrinfo_t * nbrinfo, 
    adj_t size);


#define check_ucinfo __nerstrand_check_ucinfo
/**
 * @brief Verify the sanity of a graphs ucinfo 
 *
 * @param graph The clustered graph 
 * @param ucinfo The ucinfo to check
 * @param where The cluster information for each vertex 
 * @param nclusters The number of clusters
 *
 * @return 1 on success, 0 for invalid ucinfo 
 */
int check_ucinfo(
    graph_t const * graph, 
    ucinfo_t const * const * ucinfo,
    cid_t const * const * where, 
    cid_t nclusters);




#endif
