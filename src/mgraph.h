/**
 * @file mgraph.h
 * @brief Types and function prototypes for multi-level graphs
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-10-08
 */




#ifndef NERSTRAND_MGRAPH_H
#define NERSTRAND_MGRAPH_H




#include "base.h"
#include "graph.h"
#include "ucinfo.h"




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef struct mgraph_t {
  size_t level;
  const graph_t * graph;
  vtx_t ** cmap;
  struct mgraph_t * coarser;
  struct mgraph_t * finer;
  ucinfo_t ** ucinfo;
  int free_cmap;
  int free_graph;
  int free_ucinfo;
} mgraph_t;




/******************************************************************************
* SERIAL FUNCTION PROTOTYPES **************************************************
******************************************************************************/


#define init_mgraph __nerstrand_init_mgraph
/**
 * @brief Initialize an mgraph structure
 *
 * @param mgraph The structure to initialize
 *
 * @return The initialized structure 
 */
mgraph_t * init_mgraph(
    mgraph_t * mgraph);


#define setup_mgraph __nerstrand_setup_mgraph
/**
 * @brief Setup an mgraph structure with the given parameters
 *
 * @param level The level of the mgraph (0 is finest)
 * @param graph The underlying graph structure
 * @param cmap The fine vertex to coarse vertex mapping
 * @param coarser The coarser mgraph
 * @param finer The finer mgraph
 *
 * @return A pointer to the setup mgraph structure
 */
mgraph_t * setup_mgraph(
    size_t level,
    graph_t const * graph,
    vtx_t ** cmap, 
    mgraph_t * coarser,
    mgraph_t * finer);


#define free_mgraph __nerstrand_free_mgraph
/**
 * @brief Free an mgraph structure and associated memory
 *
 * @param mgraph The mgraph to free
 */
void free_mgraph(
    mgraph_t * mgraph);


#define adjust_cmap __nerstrand_adjust_cmap
/**
 * @brief Adjust a cmap generated with a the old offset to use the new offset 
 *
 * @param cmap The cmap to adjust
 * @param nvtxs The size of the cmap
 * @param nthreads The number of parts to this cmap
 * @param olddist The graph distribution from which the cmap was created
 * @param newdist The graph distribution that the cmap will be used with
 */
void adjust_cmap(
    vtx_t * const * cmap,
    vtx_t const * nvtxs,
    tid_t nthreads,
    graphdist_t const olddist,
    graphdist_t const newdist);




/******************************************************************************
* PARALLEL FUNCTION PROTOTYPES ************************************************
******************************************************************************/


#define par_setup_mgraph __nerstrand_par_setup_mgraph
/**
 * @brief Setup an mgraph structure in parallel
 *
 * @param level The level of the new mgraph
 * @param graph The underlying graph
 * @param cmap The fine to coarse vertex mapping
 * @param coarser The coarser mgraph
 * @param finer The finer mgraph
 * @param comm The thread communicator
 *
 * @return The setup mgraph
 */
mgraph_t * par_setup_mgraph(
    size_t level,
    graph_t const * graph,
    vtx_t * cmap, 
    mgraph_t * coarser,
    mgraph_t * finer,
    dlthread_comm_t comm);


#define par_free_mgraph __nerstrand_par_free_mgraph
/**
 * @brief Free the mgraph structure in parallel
 *
 * @param mgraph The mgraph to free.
 * @param comm The thread communicator.
 */
void par_free_mgraph(
    mgraph_t * mgraph,
    dlthread_comm_t comm);


#define par_adjust_cmap __nerstrand_par_adjust_cmap
/**
 * @brief Adjust a cmap generated with a the old offset to use the new offset 
 *
 * @param cmap The cmap to adjust
 * @param mynvtxs The size of the cmap
 * @param olddist The graph distribution from which the cmap was created
 * @param newdist The graph distribution that the cmap will be used with
 */
void par_adjust_cmap(
    vtx_t * cmap, 
    vtx_t mynvtxs,
    graphdist_t olddist, 
    graphdist_t newdist);




#endif
