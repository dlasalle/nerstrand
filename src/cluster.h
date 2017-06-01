/**
 * @file cluster.h
 * @brief Types and function prototypes for mainting clusters
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-10-08
 */



#ifndef NERSTRAND_CLUSTER_H
#define NERSTRAND_CLUSTER_H




#include "base.h"
#include "graph.h"




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


/**
 * @brief The structure to hold information for each cluster
 */
typedef struct cluster_t {
  vtx_t nvtxs;
  cid_t id;
  twgt_t eadjwgt;  /* sum(adjwgt) for external facing edges */
  twgt_t iadjwgt;  /* sum(adjwgt) for internal edges */
  twgt_t viadjwgt; /* sum(iadjwgt) for vertices in cluster */
} cluster_t;


/**
 * @brief The structure to hold global clustering information
 */
typedef struct clustering_t {
  cid_t nclusters;
  cluster_t * clusters;
  cid_t ** where; 
  int free_clusters, free_where;
  tid_t npar;
} clustering_t;




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLMEM_PREFIX cluster
#define DLMEM_TYPE_T cluster_t
#include "dlmem_headers.h"
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX


#define DLMEM_PREFIX clustering
#define DLMEM_TYPE_T clustering_t
#include "dlmem_headers.h"
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX




/******************************************************************************
* SERIAL FUNCTION PROTOTYPES **************************************************
******************************************************************************/


#define init_clustering __nerstrand_init_clustering
/**
 * @brief Initialize the clustering
 *
 * @param clustering The clustering to initialize
 *
 * @return The initialized clustering
 */
clustering_t * init_clustering(
    clustering_t * clustering);


#define init_cluster __nerstrand_init_cluster
/**
 * @brief Initialize the cluster
 *
 * @param cluster The cluster to initialize
 *
 * @return The initialized cluster
 */
cluster_t * init_cluster(
    cluster_t * cluster);


#define setup_clustering __nerstrand_setup_clustering
/**
 * @brief Allocate and initialize a clustering based on supplied inputs
 *
 * @param nclusters The number of clusters in the clustering
 * @param where An two dimensional array of vertex cluster ids
 * @param graph The graph of which the clustering is on
 *
 * @return The initialized clustering 
 */
clustering_t * setup_clustering(
    cid_t nclusters, 
    cid_t ** where, 
    graph_t const * graph);


#define clone_clustering __nerstrand_clone_clustering
/**
 * @brief Create a copy of a clustering
 *
 * @param nvtxs The number of vertices per thread
 * @param clustering The clustering to copy
 *
 * @return The duplicate clustering
 */
clustering_t * clone_clustering(
    vtx_t const * nvtxs, 
    clustering_t const * clustering);


#define unify_where __nerstrand_unify_where
/**
 * @brief Create a unified where vector from a distributed one.
 *
 * @param where The distributed where vector.
 * @param graph The graph associated with the where vector.
 * @param uwhere The unified where vector.
 *
 * @return
 */
void unify_where(
    cid_t const * const * where, 
    graph_t const * graph, 
    cid_t * uwhere);


#define setup_cluster __nerstrand_setup_cluster
/**
 * @brief Initialize a cluster based on supplied inputs
 *
 * @param id The id of the cluster
 * @param nvtxs The number of vertices in the cluter
 *
 * @return The initialized cluster
 */
cluster_t * setup_cluster(
    cid_t id, 
    vtx_t nvtxs);


#define calc_clustering __nerstrand_calc_clustering
/**
 * @brief Update the calculated information on the clustering and contained
 * clusters.
 *
 * @param clustering The clustering to update.
 * @param graph The graph to update the clustering information from.
 *
 * @return !0 if successful 
 */
int calc_clustering(
    clustering_t * clustering, 
    graph_t const * graph);


#define free_clustering __nerstrand_free_clustering
/**
 * @brief Free clustering data
 *
 * @param clustering The clustering data to be free'd
 *
 * @return !0 if successful
 */
int free_clustering(
    clustering_t * clustering);


#define free_cluster __nerstrand_free_cluster
/**
 * @brief Free the cluster data
 *
 * @param cluster The cluster data to be free'd
 *
 * @return !0 if successful
 */
int free_cluster(
    cluster_t * cluster);


#define calc_modularity __nerstrand_calc_modularity
/**
 * @brief Calculates the modularity, Q(V_i), of the passed in graph and
 * corresponding cluster.
 *
 * @param cluster The cluster of which to calculate the modularity. 
 * @param graph The graph to the which the clustering has been applied.
 *
 * @return The modularity of the given cluster
 */
mod_t calc_modularity(
    cluster_t const * cluster, 
    graph_t const * graph);


#define calc_total_modularity __nerstrand_calc_total_modularity 
/**
 * @brief Calculates the modularity sum Q(V_i) of the whole clustering
 *
 * @param clustering The clustering of which to calculate the modularity.
 * @param graph The graph to which the clustering has been applied.
 *
 * @return 
 */
mod_t calc_total_modularity(
    clustering_t const * clustering,
    graph_t const * graph);


#define check_clustering __nerstrand_check_clustering
/**
 * @brief Verify the sanity of a clustering
 *
 * @param clustering The clustering to verify
 * @param graph The clustered graph
 *
 * @return 1 on success, 0 on invalid clustering
 */
int check_clustering(
    clustering_t const * clustering, 
    graph_t const * graph);




/******************************************************************************
* PARALLEL FUNCTION PROTOTYPES ************************************************
******************************************************************************/


#define par_setup_clustering __nerstrand_par_setup_clustering
/**
 * @brief Setup a clustering -- to be called by all threads in a parallel
 * section.
 *
 * @param nclusters The number of clustering in the clustering
 * @param nvtxs The number of vertices in the clustering
 * @param clusterless_vertices Whether or not there are clusterless vertices
 * @param clusters The array of clusters
 * @param where The where vector (split for each thread)
 * @param comm The thread communicator.
 *
 * @return The setup clustering
 */
clustering_t * par_setup_clustering(
    cid_t nclusters, 
    vtx_t nvtxs,
    int clusterless_vertices, 
    cluster_t * clusters, 
    cid_t ** where,
    dlthread_comm_t comm);


#define par_calc_clustering __nerstrand_par_calc_clustering
/**
 * @brief Update teh calculated infomration on the clustering and contained
 * clusters -- to be called from all threads inside of a parallel region.
 *
 * @param clustering The clustering to update>
 * @param graph The graph to update teh clustering from.
 * @param comm The thread communicator.
 */
void par_calc_clustering(
    clustering_t * clustering, 
    graph_t const * graph,
    dlthread_comm_t comm);


#define par_calc_total_modularity __nerstrand_par_calc_total_modularity
/**
 * @brief Calculate total modularity of a clustering in parallel
 *
 * @param clustering The clustering 
 * @param graph The clustered graph
 * @param comm The thread communicator.
 *
 * @return The modularity
 */
mod_t par_calc_total_modularity(
    clustering_t const * clustering,
    graph_t const * graph,
    dlthread_comm_t comm);


#define par_free_clustering __nerstrand_par_free_clustering
/**
 * @brief Free a clustering and associate memory -- to be called by all
 * threads from inside of a parallel region
 *
 * @param clustering The clustering structure.
 * @param comm The thread communicator.
 */
void par_free_clustering(
    clustering_t * clustering,
    dlthread_comm_t comm);




/******************************************************************************
* INLINE FUNCTIONS ************************************************************
******************************************************************************/


/**
 * @brief Determines if a vertex with internal degree id and external degree ed
 * is on the boundary of a cluster or not.
 *
 * @param id Internal cluster degree of the vertex
 * @param ed External cluster degree of the vertex
 * @param ie Internal degree of the vertex
 *
 * @return !0 if the vertex is on the boundary and 0 if it is not
 */
static inline int is_bnd(
    wgt_t const id, 
    wgt_t const ed, 
    wgt_t const ie)
{
  return ed > 0; //ed > dl_max(0,id - ie);
}


/**
 * @brief The difference in modularity of removing a vertex from a cluster
 *
 * @param vid The internal degree of the vertex
 * @param vcd The edgeweight between the vertex and the cluster
 * @param vtd The total degree of the vertex (excluding vid)
 * @param cd Degree of the cluster
 * @param gadjwgt Global edge weight
 *
 * @return Difference in modularity
 */
static inline mod_t calc_mod_difference(
    wgt_t const vid, 
    wgt_t const vcd, 
    wgt_t const vtd, 
    wgt_t const cd, 
    wgt_t const invgadjwgt)
{
  /* (Q_{A-v} Q_{B+v}) - (Q_{A} + Q_{B})
   * but we don't know anything about B.
   */
  wgt_t const va = vtd + vid;
  return ((cd - va)*va*invgadjwgt) - vcd;
}




#endif
