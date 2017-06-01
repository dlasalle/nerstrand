/**
 * @file contract.h
 * @brief Function prototypes for contracting graphs
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2014-03-03
 */




#ifndef NERSTRAND_CONTRACT_H
#define NERSTRAND_CONTRACT_H




#include "base.h"
#include "aggregate.h"




/******************************************************************************
* SERIAL FUNCTION PROTOTYPES **************************************************
******************************************************************************/


#define contract_graph __nerstrand_contract_graph
/**
 * @brief Contract a graph given an aggregation 
 *
 * @param objective The objective specifying how to contract the graph
 * @param mgraph The graph to contract
 * @param aggregate The aggeregation to contract based on
 *
 */
graph_t * contract_graph(
    objective_t * objective, 
    mgraph_t * mgraph, 
    aggregate_t const * aggregate);




/******************************************************************************
* PARALLEL FUNCTION PROTOTYPES ************************************************
******************************************************************************/


#define par_contract_graph __nerstrand_par_contract_graph
/**
 * @brief Contract a graph given an aggregation -- to be called by each thread
 * from within a parallel region
 *
 * @param objective The objective specifying how to contract the graph 
 * @param mgraph The graph to contract
 * @param aggregate The aggregation to contract based on
 *
 * @return The contracted (coarser) graph
 */
graph_t * par_contract_graph(
    objective_t * objective,
    mgraph_t * mgraph,
    aggregate_t const * aggregate);




#endif
