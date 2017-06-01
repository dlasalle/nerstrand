/**
 * @file ucinfo.c
 * @brief Functions for uncoarsening information
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-10-08
 */





#ifndef NERSTRAND_UCINFO_C
#define NERSTRAND_UCINFO_C




#include "ucinfo.h"




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLMEM_PREFIX nbrinfo
#define DLMEM_TYPE_T nbrinfo_t
#define DLMEM_INITFUNCTION init_nbrinfo
#define DLMEM_DLTYPE DLTYPE_STRUCT
#include "dlmem_funcs.h"
#undef DLMEM_INITFUNCTION
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMEM_PREFIX ucinfo
#define DLMEM_TYPE_T ucinfo_t
#define DLMEM_INITFUNCTION init_ucinfo
#define DLMEM_DLTYPE DLTYPE_STRUCT
#include "dlmem_funcs.h"
#undef DLMEM_INITFUNCTION
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLPQ_PREFIX rv
#define DLPQ_KEY_T real_t
#define DLPQ_VAL_T vtx_t
#include "dlpq_funcs.h"
#undef DLPQ_PREFIX
#undef DLPQ_KEY_T
#undef DLPQ_VAL_T




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


nbrinfo_t * init_nbrinfo(
    nbrinfo_t * const nbrinfo) 
{
  nbrinfo->iw = 0;
  nbrinfo->in = 0;
  nbrinfo->ew = 0;
  nbrinfo->nnbrs = 0;
  nbrinfo->nbrstart = NO_NBRS;
  return nbrinfo;
}


int free_nbrinfo(
    nbrinfo_t * nbrinfo)
{
  dl_free(nbrinfo);
  return 1;
}


ucinfo_t * init_ucinfo(
    ucinfo_t * const ucinfo)
{
  ucinfo->nadj = 0;
  ucinfo->maxnadj = 0;
  ucinfo->bnd = NULL;
  ucinfo->adjncy = NULL;
  ucinfo->adjwgt = NULL;
  ucinfo->adjnum = NULL;
  ucinfo->free_adjncy = 0;
  ucinfo->free_adjwgt = 0;
  ucinfo->free_adjnum = 0;
  return ucinfo;
}


ucinfo_t * setup_ucinfo(
    const vtx_t nvtxs, 
    nbrinfo_t * const nbrinfo, 
    cid_t * const adjncy, 
    wgt_t * const adjwgt,
    vtx_t * const adjnum)
{
  ucinfo_t * const ucinfo = ucinfo_calloc(1);

  ucinfo->bnd = vtx_iset_create(0,nvtxs);
  ucinfo->nbrinfo = nbrinfo;
  ucinfo->adjncy = adjncy;
  ucinfo->adjwgt = adjwgt;
  ucinfo->adjnum = adjnum;

  return ucinfo;
}


ucinfo_t * clone_ucinfo(
    const vtx_t nvtxs, 
    const ucinfo_t * const ucinfo)
{
  ucinfo_t * const nucinfo = ucinfo_calloc(1);

  nucinfo->maxnadj = ucinfo->maxnadj;
  nucinfo->nadj = ucinfo->nadj;

  nucinfo->bnd = vtx_iset_clone(ucinfo->bnd);

  nucinfo->nbrinfo = nbrinfo_duplicate(ucinfo->nbrinfo,nvtxs);
  nucinfo->adjncy = cid_duplicate(ucinfo->adjncy,nucinfo->maxnadj);
  nucinfo->adjwgt = wgt_duplicate(ucinfo->adjwgt,nucinfo->maxnadj);
  nucinfo->adjnum = vtx_duplicate(ucinfo->adjnum,nucinfo->maxnadj);

  nucinfo->free_adjncy = 1;
  nucinfo->free_adjwgt = 1;

  return nucinfo;
}


ucinfo_t ** clone_ucinfos(
    const vtx_t * const nvtxs, 
    const ucinfo_t * const * const ucinfos, 
    const tid_t nthreads)
{
  tid_t t;

  ucinfo_t ** const nucinfos = r_ucinfo_alloc(nthreads);

  for (t=0;t<nthreads;++t) {
    nucinfos[t] = clone_ucinfo(nvtxs[t],ucinfos[t]);
  }

  return nucinfos;
}


int free_ucinfo(
    ucinfo_t * ucinfo)
{
  dl_free(ucinfo->nbrinfo);
  vtx_iset_free(ucinfo->bnd);

  if (ucinfo->free_adjncy) {
    dl_free(ucinfo->adjncy);
  }
  if (ucinfo->free_adjwgt) {
    dl_free(ucinfo->adjwgt);
  }
  if (ucinfo->free_adjnum) {
    dl_free(ucinfo->adjnum);
  }

  dl_free(ucinfo);

  return 1;
}


int free_ucinfos(
    ucinfo_t ** ucinfo, 
    const tid_t nthreads)
{
  tid_t myid;
  for (myid=0;myid<nthreads;++myid) {
    vtx_iset_free(ucinfo[myid]->bnd);
    dl_free(ucinfo[myid]->nbrinfo);
    if (ucinfo[myid]->free_adjncy) {
      dl_free(ucinfo[myid]->adjncy);
    }
    if (ucinfo[myid]->free_adjwgt) {
      dl_free(ucinfo[myid]->adjwgt);
    }
    if (ucinfo[myid]->free_adjnum) {
      dl_free(ucinfo[myid]->adjnum);
    }
    dl_free(ucinfo[myid]);
  }

  dl_free(ucinfo);

  return 1;
}


adj_t next_nbr_adj(
    ucinfo_t * const ucinfo, 
    const adj_t size)
{
  adj_t a = ucinfo->nadj;
  ucinfo->nadj += size;

  if (ucinfo->nadj > ucinfo->maxnadj) { /* expand */
    dprintf("Expanding ucinfo->adjncy, ucinfo->adjwgt, and ucinfo->adjnum " \
        "to "PF_ADJ_T" from "PF_ADJ_T"\n",ucinfo->maxnadj,ucinfo->maxnadj*2);
    ucinfo->maxnadj*=2;
    ucinfo->adjncy = cid_realloc(ucinfo->adjncy,ucinfo->maxnadj);
    ucinfo->adjwgt = wgt_realloc(ucinfo->adjwgt,ucinfo->maxnadj);
    ucinfo->adjnum = vtx_realloc(ucinfo->adjnum,ucinfo->maxnadj);
  }

  return a;
}


int release_nbr_adj(
    ucinfo_t * const ucinfo, 
    nbrinfo_t * const nbrinfo,
    const adj_t size)
{
  DL_ASSERT_EQUALS(nbrinfo->nbrstart,ucinfo->nadj-size,PF_ADJ_T);
  ucinfo->nadj -= size;
  nbrinfo->nbrstart = NULL_ADJ; 
  return 1;
}


int check_ucinfo(
    const graph_t * const graph, 
    const ucinfo_t * const * const ucinfo, 
    const cid_t * const * const where, 
    const cid_t nclusters)
{
  vtx_t i, k, nint;
  adj_t j;
  cid_t me, other, m, nnbrs;
  wgt_t tiw, tew, ewgt;
  tid_t myid, o;

  vtx_t * nbrlst;
  vtx_t mynvtxs;

  /* expose the graph parts */
  const tid_t nthreads = graph->npar;
  const adj_t * const * const xadj = (const adj_t * const *)graph->xadj;
  const vtx_t * const * const adjncy = (const vtx_t * const *)graph->adjncy;
  const wgt_t * const * const adjwgt = (const wgt_t * const *)graph->adjwgt;

  /* the ucinfo */
  const nbrinfo_t * nbrinfo;
  const cid_t * nadjncy;
  const wgt_t * nadjwgt;
  const vtx_t * nadjnum;

  /* other stuff */
  wgt_t * htwgt = wgt_alloc(nclusters);
  vtx_t * htnum = vtx_alloc(nclusters);

  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    nbrinfo = ucinfo[myid]->nbrinfo;
    nadjncy = ucinfo[myid]->adjncy;
    nadjwgt = ucinfo[myid]->adjwgt;
    nadjnum = ucinfo[myid]->adjnum;
    nbrlst = vtx_init_alloc(NULL_VTX,ucinfo[myid]->maxnadj);
    /* verify preliminary stuff */
    for (i=0;i<mynvtxs;++i) {
      for (j=0;j<nbrinfo[i].nnbrs;++j) {
        if (nbrinfo[i].nbrstart == NULL_CID) {
          eprintf("Vertex "PF_VTX_T" thinks it has "PF_CID_T" neighbors, but "
              "its nbrstart is unset\n",i,nbrinfo[i].nnbrs);
          return 0;
        } else if (nbrlst[nbrinfo[i].nbrstart+j] != NULL_VTX) {
          eprintf("Vertex "PF_VTX_T"'s nbrlist overlaps with vertex "PF_VTX_T
              "'s\n",i,nbrlst[nbrinfo[i].nbrstart+j]);
          return 0;
        } else {
          nbrlst[nbrinfo[i].nbrstart+j] = i;
        }
      }
    }
    dl_free(nbrlst);

    /* parse the vertex list */
    for (i=0;i<mynvtxs;++i) {
      tiw = tew = 0.0;
      nint = 0;
      nnbrs = 0;
      wgt_set(htwgt,0,nclusters);
      vtx_set(htnum,0,nclusters);
      me = where[myid][i];
     
      /* find external and internal vertices */
      for (j=xadj[myid][i];j<xadj[myid][i+1];++j) {
        k = adjncy[myid][j];
        if (k < mynvtxs) {
          o = myid;
        } else {
          o = gvtx_to_tid(k,graph->dist);
          k = gvtx_to_lvtx(k,graph->dist);
        }
        if (adjwgt) {
          ewgt = adjwgt[myid][j];
        } else {
          ewgt = 1.0;
        }
        other = where[o][k];
        if (me != other) {
          tew += ewgt; 
          if ((m = htwgt[other]) == 0) {
            ++nnbrs;
          }
          htwgt[other] += ewgt;
          ++htnum[other];
        } else {
          tiw += ewgt;
          ++nint;
        }
      }
      if (!dl_near_equal(nbrinfo[i].iw,tiw) || \
          !dl_near_equal(nbrinfo[i].ew,tew)) {
        eprintf("Mismatch of cluster degree of vertex "PF_VTX_T" in " \
            PF_CID_T", nbrinfo.iw = "PF_WGT_T"/"PF_WGT_T"/"PF_WGT_T", " \
            "actual = "PF_WGT_T"/"PF_WGT_T"/"PF_WGT_T"\n",i,where[myid][i], \
            nbrinfo[i].iw,nbrinfo[i].ew,graph->eadjwgt[myid][i],tiw,tew, \
            tiw+tew);
        return 0;
      }
      if (nbrinfo[i].in != nint) {
        eprintf("Mismatch of cluster number of edges of vertex "PF_VTX_T \
            " in "PF_CID_T", nbrinfo.in = "PF_VTX_T"/"PF_VTX_T \
            ", actual = "PF_VTX_T"/"PF_VTX_T"\n",i,where[myid][i], \
            nbrinfo[i].in,(vtx_t)(xadj[myid][i+1]-xadj[myid][i]),nint,nint);
        return 0;
      }
      if (nbrinfo[i].nnbrs != nnbrs) {
        eprintf("Mismatch of number of neighbors for vertex "PF_VTX_T" in "
            PF_CID_T", nbrinfo.nnbrs = "PF_CID_T", actual = "PF_CID_T"\n",i,
            where[myid][i],nbrinfo[i].nnbrs,nnbrs);
        printf("["); 
        for (other=0;other<nclusters;++other) {
          printf(PF_CID_T":%E,",other,htwgt[other]);
        }
        printf("]\n{");
        for (j=0;j<nbrinfo[i].nnbrs;++j) {
          other = nadjncy[nbrinfo[i].nbrstart+j];
          ewgt = nadjwgt[nbrinfo[i].nbrstart+j];
          printf(PF_CID_T":%E,",other,ewgt);
        }
        printf("}\n");
        return 0;
      }
      if (nnbrs > 0 && nbrinfo[i].nbrstart == NO_NBRS) {
        eprintf("Vertex "PF_VTX_T" has an unset nbrstart\n",i);
        return 0;
      }
      for (j=0;j<nnbrs;++j) {
        other = nadjncy[nbrinfo[i].nbrstart+j];
        /* check edge numbers */
        if (htnum[other] == 0) {
          if (nadjnum[nbrinfo[i].nbrstart+j] == 0) {
            eprintf("Vertex "PF_VTX_T" (in "PF_CID_T") thinks its connect " \
                "to "PF_CID_T" with zero edges\n",i,me,other);
          } else {
            eprintf("Vertex "PF_VTX_T" (in "PF_CID_T") thinks its connect " \
                "to "PF_CID_T" with "PF_VTX_T" edges (not zero)\n",i, \
                me,other,nadjnum[nbrinfo[i].nbrstart+j]);
          }
          return 0;
        }
        if (htnum[other] != nadjnum[nbrinfo[i].nbrstart+j]) {
          eprintf("Mismatch of neighbor edges for vertex "PF_VTX_T" to "
              "cluster "PF_CID_T" of "PF_VTX_T", actual = "PF_VTX_T"\n",i,
              other,nadjnum[nbrinfo[i].nbrstart+j],htnum[other]);
          return 0;
        }
        /* check edge weights */
        if (htwgt[other] == 0) {
          if (nadjwgt[nbrinfo[i].nbrstart+j] == 0) {
            eprintf("Vertex "PF_VTX_T" (in "PF_CID_T") thinks its connect to "
                PF_CID_T" with a weight of ZERO\n",i,me,other);
          } else {
            eprintf("Vertex "PF_VTX_T" (in "PF_CID_T") thinks its connect to "
                PF_CID_T" with a weight of "PF_WGT_T" (not zero)\n",i,me,other,
                nadjwgt[nbrinfo[i].nbrstart+j]);
          }
          return 0;
        }
        if (!dl_near_equal(htwgt[other],nadjwgt[nbrinfo[i].nbrstart+j])) {
          eprintf("Mismatch of neighbor weight for vertex "PF_VTX_T" to "
              "cluster "PF_CID_T" of "PF_WGT_T", actual = "PF_WGT_T"\n",i,
              other,nadjwgt[nbrinfo[i].nbrstart+j],htwgt[other]);
          return 0;
        }
      }
    }
  }

  dl_free(htwgt);
  dl_free(htnum);

  return 1;
}




#endif
