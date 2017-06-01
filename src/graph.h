/**
 * @file graph.h
 * @brief Types and function prototypes for graphs
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-10-08
 */




#ifndef NERSTRAND_GRAPH_H
#define NERSTRAND_GRAPH_H




#include "base.h"




/******************************************************************************
* MACROS **********************************************************************
******************************************************************************/


#define is_supernode(i,myid,graph,alpha) \
  ((graph)->eadjwgt[myid][i]*(alpha) > (graph)->snthresh)





/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


/**
 * @brief The structure holding information for distribing vertices among
 * threads.
 */
typedef struct graphdist_t {
  vtx_t mask;
  int shift;
  vtx_t offset;
} graphdist_t;


/**
 * @brief The structure to hold the graph information.
 */
typedef struct graph_t {
  vtx_t nvtxs; /* number of vertices */
  adj_t nedges; /* number of edges */
  twgt_t tadjwgt; /* sum(adjwgt) */
  twgt_t gadjwgt; /* sum(adjwgt) + sum(iadjwgt) */
  graphdist_t dist;
  /* thread specific stuff */
  adj_t * maxdeg; /* maximum degree */
  vtx_t * mynvtxs;
  adj_t * mynedges;
  adj_t ** xadj;
  vtx_t ** adjncy;
  wgt_t ** adjwgt, **iadjwgt, **eadjwgt;
  vtx_t ** alias;
  int phantom_edges;
  int free_xadj, free_adjncy, free_adjwgt, free_iadjwgt;
  tid_t npar; /* number of thread parts */
  wgt_t snthresh;
} graph_t;




/******************************************************************************
* SERIAL FUNCTION PROTOTYPES **************************************************
******************************************************************************/


#define init_graph __nerstrand_init_graph
/**
 * @brief Initialize (zero) the graph
 *
 * @param graph The graph to initialize
 *
 * @return The initialized graph
 */
graph_t * init_graph(
    graph_t * graph);


#define setup_graph __nerstrand_setup_graph
/**
 * @brief Setup a graph with passed in edge and vertex arrays 
 *
 * @param nvtxs The number of vertices owned by each thread
 * @param xadj The adjncy array indicies for each vertex
 * @param adjncy The list of edges for each vertex as indexed by xadj
 * @param adjwgt The weights of the edges (indexed the same as adjancy)
 * @param iadjwgt The weights of vertex internal edges
 * @param eadjwgt The sum of edge weights per vertex (weighted degree)
 * @param maxdeg The max vertex degree for each thread
 * @param nthreads The number of threads that will be used with this graph
 *
 * @return The setup graph structure.
 */
graph_t * setup_graph(
    vtx_t * nvtxs, 
    adj_t ** xadj, 
    vtx_t ** adjncy, 
    wgt_t ** adjwgt, 
    wgt_t ** iadjwgt, 
    wgt_t ** eadjwgt, 
    adj_t * maxdeg, 
    tid_t nthreads);


#define distribue_graph __nerstrand_distribute_graph 
/**
 * @brief Creates a distributed graph structure using the specified
 * distribution strategy. 
 *
 * @param nvtxs The number of vertices in the graph (|xadj| = nvtxs + 1)
 * @param xadj The adjncy array indices for each vertex
 * @param adjncy The list of edges for each vertex as indexed by xadj
 * @param adjwgt The weights of the edges (indexed the same as adjancy)
 * @param iadjwgt The weights of vertex internal edges
 * @param eadjwgt The sum of edge weights per vertex (weighted degree)
 * @param nthreads The number of threads that will be used with this graph
 * @param distribution The distribution method to use.
 * @param block The size of the block for block cyclic distributions
 *
 * @return The initialized graph
 */
graph_t * distribute_graph(
    vtx_t nvtxs, 
    adj_t const * xadj, 
    vtx_t const * adjncy, 
    wgt_t const * adjwgt, 
    wgt_t const * iadjwgt, 
    wgt_t const * eadjwgt, 
    tid_t nthreads, 
    int distribution,
    vtx_t block);


#define free_graph __nerstrand_free_graph
/**
 * @brief Frees a graph and it's memory. It does not free arrays
 * (xadj,adjncy,adjwgt) unless the corresponding flag is set, free_*
 *
 * @param graph The graph to be free'd
 *
 * @return !0 if successful 
 */
int free_graph(
    graph_t * graph);


#define calc_graph_dist __nerstrand_calc_graph_dist
/**
 * @brief Populate a graphdist_t struct 
 *
 * @param maxnvtxs Maximum number of vertices for any thread
 * @param dist The distribution struct
 *
 * @return !0 for success
 */
int calc_graph_dist(
    vtx_t maxnvtxs, 
    graphdist_t * dist);


#define check_graph __nerstrand_check_graph
/**
 * @brief Verify the sanity of a graph structure
 *
 * @param graph The graph to verify
 *
 * @return 1 on success, 0 on invalid graph
 */
int check_graph(
    graph_t const * graph);


#define build_reverse_adjidx __nerstrand_build_reverse_adjidx
/**
 * @brief Builds a two dimensional reverse index of the edges
 *
 * @param graph The graph to build the reverse index from
 *
 * @return The two dimensional reverse index of edges
 */
adj_t ** build_reverse_adjidx(
    graph_t const * graph);




/******************************************************************************
* PARALLEL FUNCTION PROTOTYPES ************************************************
******************************************************************************/


#define par_setup_graph __nerstrand_par_setup_graph
/**
 * @brief Initializes a graph with the passed in edges and vertices
 *
 * @param nvtxs The number of vertices in the graph (|xadj| = nvtxs + 1)
 * @param xadj The adjncy array indices for each vertex
 * @param adjncy The list of edges for each vertex as indexed by xadj
 * @param adjwgt The weights of the edges (indexed the same as adjancy)
 * @param iadjwgt The weights of vertex internal edges
 * @param eadjwgt The sum of edge weights per vertex (weighted degree)
 * @param maxdeg The maximum vertex (NULL_VTX if unknown)
 * @param comm The thread communicator.
 *
 * @return The initialized graph
 */
graph_t * par_setup_graph(
    vtx_t nvtxs, 
    adj_t * xadj, 
    vtx_t * adjncy, 
    wgt_t * adjwgt, 
    wgt_t * iadjwgt, 
    wgt_t * eadjwgt, 
    adj_t maxdeg,
    dlthread_comm_t comm);


#define par_free_graph __nerstrand_par_free_graph
/**
 * @brief Frees a graph and it's memory. It does not free arrays
 * (xadj,adjncy,adjwgt) unless the corresponding flag is set, free_*
 *
 * @param graph The graph to be free'd
 * @param comm The thread communicator
 */
void par_free_graph(
    graph_t * graph,
    dlthread_comm_t comm);



#endif
