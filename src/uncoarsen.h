/**
 * @file uncoarsen.h
 * @brief Function prototypes for uncoarsening a graph
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2014-03-03
 */




#ifndef NERSTRAND_UNCOARSEN_H
#define NERSTRAND_UNCOARSEN_H




#include "base.h"
#include "cluster.h"
#include "ucinfo.h"
#include "mgraph.h"
#include "objective.h"




/******************************************************************************
* SERIAL FUNCTION PROTOYPES ***************************************************
******************************************************************************/


#define derive_ucinfo __nerstrand_derive_ucinfo
/**
 * @brief Create and return the ucinfo for the supplied graph and clustering 
 *
 * @param clustering The clustering to use while deriving the ucinfo
 * @param mgraph The mgraph to use while deriving the ucinfo
 *
 * @return The derived ucinfo
 */
ucinfo_t ** derive_ucinfo(
    clustering_t const * clustering, 
    mgraph_t const * mgraph);



#define uncoarsen_graph __nerstrand_uncoarsen_graph
/**
 * @brief Project and refine a clustering from the coarser graph cmgraph to the
 * finer graph fmgraph
 *
 * @param objective The objective to use for projection and refinement
 * @param clustering The clustering to project and refine
 * @param fmgraph The fine graph
 * @param cmgraph The coarse graph
 *
 * @return The updated clustering 
 */
clustering_t * uncoarsen_graph(
    objective_t * objective, 
    clustering_t * clustering, 
    mgraph_t * fmgraph, 
    mgraph_t const * cmgraph);


#define check_bnd __nerstrand_check_bnd
/**
 * @brief Verify the sanity of the boundary
 *
 * @param mgraph The mgraph with the boundary info
 * @param where The where vector describing the clustering
 *
 * @return 1 on success, 0 for invalid boundary info
 */
int check_bnd(
    mgraph_t const * mgraph, 
    cid_t const * const * where);




/******************************************************************************
* PARALLEL FUNCTION PROTOYPES *************************************************
******************************************************************************/


#define par_derive_ucinfo __nerstrand_par_derive_ucinfo
/**
 * @brief Allocate and calculate a new ucinfo based on the graph and clustering
 * -- to be called by each thread in a parallel region
 *
 * @param clustering The clustering to calculate the ucinfo from
 * @param mgraph The graph to which the ucinfo belongs
 * @param comm The thread communicator.
 *
 * @return The newly created ucinfo 
 */
ucinfo_t * par_derive_ucinfo(
    clustering_t const * clustering, 
    mgraph_t const * mgraph,
    dlthread_comm_t comm);


#define par_uncoarsen_graph __nerstrand_par_uncoarsen_graph
/**
 * @brief Project and refine a clustering from the coarser graph cmgraph to the
 * finer graph fmgraph -- to be called by each thread in a parallel region
 *
 * @param objective The objective to use for projection and refinement
 * @param clustering The clustering to project and refine
 * @param fmgraph The fine graph
 * @param cmgraph The coarse graph
 *
 * @return The updated clustering 
 */
clustering_t * par_uncoarsen_graph(
    objective_t * objective, 
    clustering_t * clustering, 
    mgraph_t * fmgraph, 
    mgraph_t const * cmgraph);




#endif
