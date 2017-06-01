/**
 * @file io.h
 * @brief Function prototypes for input/output
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2014-03-03
 */




#ifndef NERSTRAND_IO_H
#define NERSTRAND_IO_H




#include "base.h"
#include "objective.h"
#include "graph.h"
#include "cluster.h"




/******************************************************************************
* DEFINES *********************************************************************
******************************************************************************/


#define NERSTRAND_BINARY_CLUSTER_HEADER 0x646f6d63
#define NERSTRAND_BINARY_CLUSTER_VERSION 1




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define print_clustering __nerstrand_print_clustering
/**
 * @brief Prints the clustering information for the graph
 *
 * @param clustering The clustering to print
 * @param graph The graph that has been clustered (maybe NULL if update = 0)
 * @param update Whether or not to update the clustering based on the graph
 * @param verbosity How verbose the print out should be
 */
void print_clustering(
    clustering_t * clustering, 
    graph_t const * graph,
    int update, 
    int verbosity);


#define print_timers __nerstrand_print_timers
/**
 * @brief Prints out the timers.
 *
 * @param objective Objective containing the timers to be print.
 */
void print_timers(
    objective_t const * objective);


#define print_objective __nerstrand_print_objective
/**
 * @brief Prints out the the different paramters of the objective struct.
 *
 * @param objective Objective of which to print the parameters.
 */
void print_objective(
    objective_t const * objective);


#define print_loadbalance_stats __nerstrand_print_loadbalance_stats
/**
 * @brief Print statistic loadbalancing from the various stages of clustering.
 *
 * @param stat The loadbalancing statistics
 */
void print_loadbalance_stats(
    loadbalance_stat_list_t * stat);


#define print_aggregation_stats __nerstrand_print_aggregation_stats
/**
 * @brief Print information about aggregation process. 
 *
 * @param stat The aggregation stats.
 */
void print_aggregation_stats(
    aggregation_stat_list_t * stat);


#define print_refinement_stats __nerstrand_print_refinement_stats
/**
 * @brief Print information about the refinement process. 
 *
 * @param stat The refinement stats.
 */
void print_refinement_stats(
    refinement_stat_list_t * stat);


#define print_initialclustering_stats __nerstrand_print_initialclustering_stats
/**
 * @brief Print information about the initial clustering process.
 *
 * @param stats The initial clustering stats.
 */
void print_initialclustering_stats(
    initialclustering_stats_t const * stats);


#define print_modularity_stats __nerstrand_print_modularity_stats
/**
 * @brief Print statistics on the generated modularities (mean, max, etc.)
 *
 * @param mods The series of modularities found.
 * @param nmods The number of modularities found.
 */
void print_modularity_stats(
    mod_t const * mods, 
    size_t nmods);





#endif
