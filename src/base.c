/**
 * @file base.c
 * @brief Miscellaneous junk.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2014-08-06
 */



#ifndef NERSTRAND_TYPES_C
#define NERSTRAND_TYPES C




#include "base.h"




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


/* vtx_t */
#define DLMEM_PREFIX vtx
#define DLMEM_TYPE_T vtx_t
#define DLMEM_DLTYPE DLTYPE_INTEGRAL
#include "dlmem_funcs.h"
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX vtx
#define DLMATH_TYPE_T vtx_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#include "dlmath_funcs.h"
#undef DLMATH_DLTYPE 
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


#define DLRAND_PREFIX vtx
#define DLRAND_TYPE_T vtx_t
#define DLRAND_DLTYPE DLTYPE_INTEGRAL
#include "dlrand_funcs.h"
#undef DLRAND_DLTYPE
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T


#define DLSTATS_PREFIX vtx
#define DLSTATS_TYPE_T vtx_t
#include "dlstats_funcs.h"
#undef DLSTATS_PREFIX
#undef DLSTATS_TYPE_T


#define DLISET_PREFIX vtx
#define DLISET_TYPE_T vtx_t
#include "dliset_funcs.h"
#undef DLISET_PREFIX
#undef DLISET_TYPE_T


/* adj_t */
#define DLMEM_PREFIX adj
#define DLMEM_TYPE_T adj_t
#define DLMEM_DLTYPE DLTYPE_INTEGRAL
#include "dlmem_funcs.h"
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX adj
#define DLMATH_TYPE_T adj_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#include "dlmath_funcs.h"
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


#define DLRAND_PREFIX adj
#define DLRAND_TYPE_T adj_t
#define DLRAND_DLTYPE DLTYPE_INTEGRAL
#include "dlrand_funcs.h"
#undef DLRAND_DLTYPE
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T


#define DLSTATS_PREFIX adj
#define DLSTATS_TYPE_T adj_t
#include "dlstats_funcs.h"
#undef DLSTATS_PREFIX
#undef DLSTATS_TYPE_T


/* cid_t */
#define DLMEM_PREFIX cid
#define DLMEM_TYPE_T cid_t
#define DLMEM_DLTYPE DLTYPE_INTEGRAL
#include "dlmem_funcs.h"
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX cid
#define DLMATH_TYPE_T cid_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#include "dlmath_funcs.h"
#undef DLMATH_DLTYPE 
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


#define DLRAND_PREFIX cid
#define DLRAND_TYPE_T cid_t
#define DLRAND_DLTYPE DLTYPE_INTEGRAL
#include "dlrand_funcs.h"
#undef DLRAND_DLTYPE
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T


#define DLSTATS_PREFIX cid
#define DLSTATS_TYPE_T cid_t
#include "dlstats_funcs.h"
#undef DLSTATS_PREFIX
#undef DLSTATS_TYPE_T


/* tid_t */
#define DLMEM_PREFIX tid
#define DLMEM_TYPE_T tid_t
#define DLMEM_DLTYPE DLTYPE_INTEGRAL
#include "dlmem_funcs.h"
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX tid
#define DLMATH_TYPE_T tid_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#include "dlmath_funcs.h"
#undef DLMATH_DLTYPE 
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


/* real_t */
#define DLMEM_PREFIX real
#define DLMEM_TYPE_T real_t
#define DLMEM_DLTYPE DLTYPE_FLOAT
#include "dlmem_funcs.h"
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX real
#define DLMATH_TYPE_T real_t
#define DLMATH_DLTYPE DLTYPE_FLOAT
#include "dlmath_funcs.h"
#undef DLMATH_DLTYPE 
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


/* mod_t */
#define DLMEM_PREFIX mod
#define DLMEM_TYPE_T mod_t
#define DLMEM_DLTYPE DLTYPE_FLOAT
#include "dlmem_funcs.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX mod
#define DLMATH_TYPE_T mod_t
#define DLMATH_DLTYPE DLTYPE_FLOAT
#include "dlmath_funcs.h"
#undef DLMATH_DLTYPE 
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


#define DLSTATS_PREFIX mod
#define DLSTATS_TYPE_T mod_t
#include "dlstats_funcs.h"
#undef DLSTATS_PREFIX
#undef DLSTATS_TYPE_T


/* wgt_t */
#define DLMEM_PREFIX wgt
#define DLMEM_TYPE_T wgt_t
#define DLMEM_DLTYPE DLTYPE_FLOAT
#include "dlmem_funcs.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX wgt
#define DLMATH_TYPE_T wgt_t
#define DLMATH_DLTYPE DLTYPE_FLOAT
#include "dlmath_funcs.h"
#undef DLMATH_DLTYPE 
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


#define DLRAND_PREFIX wgt
#define DLRAND_TYPE_T wgt_t
#define DLRAND_DLTYPE DLTYPE_FLOAT
#include "dlrand_funcs.h"
#undef DLRAND_DLTYPE
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T


#define DLSTATS_PREFIX wgt
#define DLSTATS_TYPE_T wgt_t
#include "dlstats_funcs.h"
#undef DLSTATS_PREFIX
#undef DLSTATS_TYPE_T


/* twgt_t */
#define DLMEM_PREFIX twgt
#define DLMEM_TYPE_T twgt_t
#define DLMEM_DLTYPE DLTYPE_FLOAT
#include "dlmem_funcs.h"
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMATH_PREFIX twgt
#define DLMATH_TYPE_T twgt_t
#define DLMATH_DLTYPE DLTYPE_FLOAT
#include "dlmath_funcs.h"
#undef DLMATH_DLTYPE 
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T


/* small_t */
#define DLMEM_PREFIX small
#define DLMEM_TYPE_T small_t
#define DLMEM_DLTYPE DLTYPE_INTEGRAL
#include "dlmem_funcs.h"
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


DL_MK_SORTKV_FUNCS(vtx,vtx_t,vtx_t)




#endif
