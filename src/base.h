/**
 * @file base.h
 * @brief Base types etc. for Nerstrand
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2014-02-13
 */




#ifndef NERSTRAND_BASE_H
#define NERSTRAND_BASE_H




/******************************************************************************
* EXTERNAL INCLUDES ***********************************************************
******************************************************************************/

/* external libraires */
#include <domlib.h>

/* c includes */
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <float.h>
#include <inttypes.h>




/******************************************************************************
* MACRO INCLUDES **************************************************************
******************************************************************************/
/* theses are needed before any real code gets included -- and should not 
 * refernce base.h */
#include "strings.h"




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


#ifdef NERSTRAND_GRAPH_TYPES_DEFINED
/* define these types for internal use */
#ifdef NERSTRAND_64BIT_VERTICES
typedef uint64_t vtx_t;
#else
typedef uint32_t vtx_t;
#endif
#ifdef NERSTRAND_64BIT_EDGES
typedef uint64_t adj_t;
#else
typedef uint32_t adj_t;
#endif
#ifdef NERSTRAND_DOUBLE_WEIGHTS
typedef double wgt_t;
#else
typedef float wgt_t;
#endif
#endif /* NERSTRAND_GRAPH_TYPES_DEFINED */


#ifdef NERSTRAND_64BIT_THREADS
typedef uint64_t tid_t;
#else
typedef uint32_t tid_t;
#endif


#ifdef NERSTRAND_DOUBLE_MODULARITY
typedef double mod_t;
#else
typedef float mod_t;
#endif



typedef uint16_t small_t;
typedef double real_t;
typedef long double twgt_t;


/**
 * @brief Flags specifying what weights are in a file
 */
typedef enum wgt_flag_t {
  ADJWGT_FLAG = 1,
  VWGT_FLAG = 1 << 1,
  IADJWGT_FLAG = 1 << 2
} wgt_flag_t;


/* the graph types need to be defined before this is included */
#include <nerstrand.h>



/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


/* macros */
#define YES_POINTER ((void*)-1)
#define DEF_NULL_VTX ((vtx_t)-1)
#define DEF_NULL_CID ((cid_t)-1)
#define DEF_NULL_WGT ((wgt_t)-1.0)
#define DEF_NULL_ADJ ((adj_t)-1)
#define DEF_NULL_TID ((tid_t)-1)


/* type null values */
static const vtx_t NULL_VTX = DEF_NULL_VTX;
static const cid_t NULL_CID = DEF_NULL_CID;
static const wgt_t NULL_WGT = DEF_NULL_WGT;
static const adj_t NULL_ADJ = DEF_NULL_ADJ;
static const tid_t NULL_TID = DEF_NULL_TID;


/* thread specific constants */
static const vtx_t BLOCKSIZE = 0x1000;
static const int BLOCKSHIFT = 12;
static const vtx_t BLOCKMASK = 0x0FFF;




/******************************************************************************
* DOMLIB FUNCTION PROTOTYPES **************************************************
******************************************************************************/


/* vtx_t */
#define DLMEM_PREFIX vtx
#define DLMEM_TYPE_T vtx_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX vtx
#define DLMATH_TYPE_T vtx_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#include "dlmath_headers.h"
#undef DLMATH_DLTYPE 
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


#define DLRAND_PREFIX vtx
#define DLRAND_TYPE_T vtx_t
#include "dlrand_headers.h"
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T


#define DLSTATS_PREFIX vtx
#define DLSTATS_TYPE_T vtx_t
#include "dlstats_headers.h"
#undef DLSTATS_PREFIX
#undef DLSTATS_TYPE_T


#define DLISET_PREFIX vtx
#define DLISET_TYPE_T vtx_t
#include "dliset_headers.h"
#undef DLISET_PREFIX
#undef DLISET_TYPE_T


#define DLTHREAD_PREFIX vtx
#define DLTHREAD_TYPE_T vtx_t
#define DLTHREAD_STATIC 1
#include "dlthread_reduction_headers.h"
#undef DLTHREAD_STATIC
#undef DLTHREAD_TYPE_T
#undef DLTHREAD_PREFIX


/* adj_t */
#define DLMEM_PREFIX adj
#define DLMEM_TYPE_T adj_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX adj
#define DLMATH_TYPE_T adj_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#include "dlmath_headers.h"
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


#define DLRAND_PREFIX adj
#define DLRAND_TYPE_T adj_t
#include "dlrand_headers.h"
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T


#define DLSTATS_PREFIX adj
#define DLSTATS_TYPE_T adj_t
#include "dlstats_headers.h"
#undef DLSTATS_PREFIX
#undef DLSTATS_TYPE_T


#define DLTHREAD_PREFIX adj
#define DLTHREAD_TYPE_T adj_t
#define DLTHREAD_STATIC 1
#include "dlthread_reduction_headers.h"
#undef DLTHREAD_STATIC
#undef DLTHREAD_TYPE_T
#undef DLTHREAD_PREFIX




/* cid_t */
#define DLMEM_PREFIX cid
#define DLMEM_TYPE_T cid_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX cid
#define DLMATH_TYPE_T cid_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#include "dlmath_headers.h"
#undef DLMATH_DLTYPE 
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


#define DLRAND_PREFIX cid
#define DLRAND_TYPE_T cid_t
#include "dlrand_headers.h"
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T


#define DLSTATS_PREFIX cid
#define DLSTATS_TYPE_T cid_t
#include "dlstats_headers.h"
#undef DLSTATS_PREFIX
#undef DLSTATS_TYPE_T


#define DLTHREAD_PREFIX cid
#define DLTHREAD_TYPE_T cid_t
#define DLTHREAD_STATIC 1
#include "dlthread_reduction_headers.h"
#undef DLTHREAD_STATIC
#undef DLTHREAD_TYPE_T
#undef DLTHREAD_PREFIX




/* tid_t */
#define DLMEM_PREFIX tid
#define DLMEM_TYPE_T tid_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX tid
#define DLMATH_TYPE_T tid_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#include "dlmath_headers.h"
#undef DLMATH_DLTYPE 
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


/* real_t */
#define DLMEM_PREFIX real
#define DLMEM_TYPE_T real_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX real
#define DLMATH_TYPE_T real_t
#define DLMATH_DLTYPE DLTYPE_FLOAT
#include "dlmath_headers.h"
#undef DLMATH_DLTYPE 
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


/* mod_t */
#define DLMEM_PREFIX mod
#define DLMEM_TYPE_T mod_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX mod
#define DLMATH_TYPE_T mod_t
#define DLMATH_DLTYPE DLTYPE_FLOAT
#include "dlmath_headers.h"
#undef DLMATH_DLTYPE 
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


#define DLSTATS_PREFIX mod
#define DLSTATS_TYPE_T mod_t
#include "dlstats_headers.h"
#undef DLSTATS_PREFIX
#undef DLSTATS_TYPE_T


#define DLTHREAD_PREFIX mod
#define DLTHREAD_TYPE_T mod_t
#define DLTHREAD_STATIC 1
#include "dlthread_reduction_headers.h"
#undef DLTHREAD_STATIC
#undef DLTHREAD_TYPE_T
#undef DLTHREAD_PREFIX




/* wgt_t */
#define DLMEM_PREFIX wgt
#define DLMEM_TYPE_T wgt_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX wgt
#define DLMATH_TYPE_T wgt_t
#define DLMATH_DLTYPE DLTYPE_FLOAT
#include "dlmath_headers.h"
#undef DLMATH_DLTYPE 
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


#define DLRAND_PREFIX wgt
#define DLRAND_TYPE_T wgt_t
#include "dlrand_headers.h"
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T


#define DLSTATS_PREFIX wgt
#define DLSTATS_TYPE_T wgt_t
#include "dlstats_headers.h"
#undef DLSTATS_PREFIX
#undef DLSTATS_TYPE_T


#define DLTHREAD_PREFIX wgt 
#define DLTHREAD_TYPE_T wgt_t
#define DLTHREAD_STATIC 1
#include "dlthread_reduction_headers.h"
#undef DLTHREAD_STATIC
#undef DLTHREAD_TYPE_T
#undef DLTHREAD_PREFIX




/* twgt_t */
#define DLMEM_PREFIX twgt
#define DLMEM_TYPE_T twgt_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX twgt
#define DLMATH_TYPE_T twgt_t
#define DLMATH_DLTYPE DLTYPE_FLOAT
#include "dlmath_headers.h"
#undef DLMATH_DLTYPE 
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


#define DLTHREAD_PREFIX twgt
#define DLTHREAD_TYPE_T twgt_t
#define DLTHREAD_STATIC 1
#include "dlthread_reduction_headers.h"
#undef DLTHREAD_STATIC
#undef DLTHREAD_TYPE_T
#undef DLTHREAD_PREFIX




/* double */
#define DLMEM_PREFIX double
#define DLMEM_TYPE_T double
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


/* small_t */
#define DLMEM_PREFIX small
#define DLMEM_TYPE_T small_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


/* long double */
#define DLTHREAD_PREFIX ld
#define DLTHREAD_TYPE_T long double
#define DLTHREAD_STATIC 1
#include "dlthread_reduction_headers.h"
#undef DLTHREAD_STATIC
#undef DLTHREAD_TYPE_T
#undef DLTHREAD_PREFIX



DL_MK_SORTKV_HEADERS(vtx,vtx_t,vtx_t)




/******************************************************************************
* INLINE FUNCTIONS ************************************************************
******************************************************************************/


/**
 * @brief Convert a global vertex number to a thread local vertex number.
 *
 * @param v The global vertex number.
 * @param mask The bitmask used for conversion.
 *
 * @return The local vertex number.
 */
static inline vtx_t gvtx_to_lvtx(
    vtx_t const v, 
    vtx_t const mask)
{
  DL_ASSERT(mask > 0,"The mask is set to 0!\n");
  DL_ASSERT(v > mask,"Global vertex number is smaller than mask (gvtx = "
      PF_VTX_T", mask = "PF_VTX_T")\n",v,mask);

  return v & mask;
}


/**
 * @brief Convert a thread local vertex number to a global vertex number.
 *
 * @param v The local vertex number.
 * @param t The owning thread.
 * @param shift The bitshift used for conversion.
 *
 * @return The gloabal vertex number.
 */
static inline vtx_t lvtx_to_gvtx(
    vtx_t const v, 
    tid_t const t, 
    int const shift)
{
  DL_ASSERT(shift > 0,"The mask size is set to 0!\n");
  DL_ASSERT(v < (vtx_t)(1 << shift),"Local vertex number is greater than "
      "shift (lvtx = "PF_VTX_T", shift = %d)\n",v,shift);

  return ((t+1) << shift) | v;
}


/**
 * @brief Find the thread ID given a global vertex number.
 *
 * @param v The global vertex.
 * @param shift THe bitshift used for conversion.
 *
 * @return The owning thread ID.
 */
static inline tid_t gvtx_to_tid(
    vtx_t const v, 
    int const shift)
{
  DL_ASSERT(shift > 0,"The shift size is set to %d!\n",shift);
  DL_ASSERT(v >= (vtx_t)(1 << shift),"Global vertex number is too small "
      "(gvtx = "PF_VTX_T", shift = %d)\n",v,shift);

  return (v >> shift)-1;
}


/**
 * @brief Determine the maximum global vertex number.
 *
 * @param shift The bitshift to use for conversion.
 * @param nthreads THe number of threads.
 *
 * @return The maximum global vertex number. 
 */
static inline vtx_t max_gvtx(
    int const shift, 
    tid_t const nthreads) 
{
  return (vtx_t)(1 << shift)*(nthreads+1);
}


/* avoid having to pass each element */
#define gvtx_to_lvtx(v,dist) gvtx_to_lvtx(v,(dist).mask)
#define lvtx_to_gvtx(v,t,dist) lvtx_to_gvtx(v,t,(dist).shift)
#define gvtx_to_tid(v,dist) gvtx_to_tid(v,(dist).shift)
#define max_gvtx(graph) max_gvtx((graph)->dist.shift,(graph)->npar)




#endif
