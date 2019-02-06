/**
 * @file nerstrand_bin.c
 * @brief Main driver function
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014-2015, Regents of the University of Minnesota
 * Copyright 2019, Dominique LaSalle
 * @version 1
 * @date 2014-09-16
 */




#ifndef NERSTRAND_BIN_C
#define NERSTRAND_BIN_C




#include <nerstrand.h>
#include <wildriver.h>
#include <domlib.h>
#include "base.h"
#include "graph.h"
#include "cluster.h"
#include "objective.h"
#include "io.h"



/******************************************************************************
* MACROS **********************************************************************
******************************************************************************/


#define __ARRAY_SIZE(a) \
  (sizeof(a) > 0 ? (sizeof(a) / sizeof((a)[0])) : 0)




/******************************************************************************
* OPTIONS *********************************************************************
******************************************************************************/


static const cmd_opt_pair_t AGGREGATE_CHOICES[] = {
  {NERSTRAND_STR_AGGREGATE_RM,"Random Matching",NERSTRAND_AGGREGATE_RM},
  {NERSTRAND_STR_AGGREGATE_RC,"Random Clustering",NERSTRAND_AGGREGATE_RC},
  {NERSTRAND_STR_AGGREGATE_AGM,"Agglomerative Modularity Matching",
      NERSTRAND_AGGREGATE_AGM},
  {NERSTRAND_STR_AGGREGATE_AGH,"Agglomerative Modularity Matching 2-Hop",
      NERSTRAND_AGGREGATE_AGH},
  {NERSTRAND_STR_AGGREGATE_AGC,"Agglomerative Modularity Clustering",
      NERSTRAND_AGGREGATE_AGC}
};


static const cmd_opt_pair_t REFINEMENT_CHOICES[] = {
  {NERSTRAND_STR_REFINEMENT_GREEDY,"Greedy Modularity Refinement.",
      NERSTRAND_REFINEMENT_GREEDY},
  {NERSTRAND_STR_REFINEMENT_RANDOM,"Random Modularity Refinement.",
      NERSTRAND_REFINEMENT_RANDOM}
};


static const cmd_opt_pair_t INITIAL_CLUSTERING_CHOICES[] = {
  {NERSTRAND_STR_INITIAL_CLUSTERING_RANDOM,"Random Clustering.",
      NERSTRAND_INITIAL_CLUSTERING_RANDOM},
  {NERSTRAND_STR_INITIAL_CLUSTERING_VTX,"Initialize each vertex as a cluster.",
      NERSTRAND_INITIAL_CLUSTERING_VTX},
  {NERSTRAND_STR_INITIAL_CLUSTERING_RVTX,"Initialize each vertex as a cluster "
      "and perform random refinement.",NERSTRAND_INITIAL_CLUSTERING_RVTX},
  {NERSTRAND_STR_INITIAL_CLUSTERING_BFS,"BFS Clustering.",
      NERSTRAND_INITIAL_CLUSTERING_BFS},
  {NERSTRAND_STR_INITIAL_CLUSTERING_SEED,"Seed Based Clustering.",
      NERSTRAND_INITIAL_CLUSTERING_SEED},
  {NERSTRAND_STR_INITIAL_CLUSTERING_GROW,"Greedy Cluster Growing.",
      NERSTRAND_INITIAL_CLUSTERING_GROW},
  {NERSTRAND_STR_INITIAL_CLUSTERING_GROWKL,"Greedy Cluster Growing with hill "
     "climbing.",NERSTRAND_INITIAL_CLUSTERING_GROWKL},
  {NERSTRAND_STR_INITIAL_CLUSTERING_LP,"Generate a clustering using Label "
      "Propagation.",NERSTRAND_INITIAL_CLUSTERING_LP},
  {NERSTRAND_STR_INITIAL_CLUSTERING_NEWMAN,"Generate a dendrogram using "
      "Newman's agglomerative method.",NERSTRAND_INITIAL_CLUSTERING_NEWMAN}
};


static const cmd_opt_pair_t SPARSIFY_CHOICES[] = {
  {NERSTRAND_STR_SPARSIFY_NONE,"Do not sparsen the graph.",
      NERSTRAND_SPARSIFY_NONE},
  {NERSTRAND_STR_SPARSIFY_RANDOM,"Sparsen the graph randomly.",
      NERSTRAND_SPARSIFY_RANDOM},
  {NERSTRAND_STR_SPARSIFY_LIGHT,"Sparsen the graph by removing light edges.",
      NERSTRAND_SPARSIFY_LIGHT},
  {NERSTRAND_STR_SPARSIFY_HEAVY,"Sparsen the graph by removing heavy edges.",
      NERSTRAND_SPARSIFY_HEAVY},
  {NERSTRAND_STR_SPARSIFY_DEGREE,"Sparsen the graph by removing high degree "
      "edges.",NERSTRAND_SPARSIFY_DEGREE}
};


static const cmd_opt_pair_t EDGEREMOVAL_CHOICES[] = {
  {NERSTRAND_STR_EDGEREMOVAL_DROP,"Drop the edge weight.",
      NERSTRAND_EDGEREMOVAL_DROP},
  {NERSTRAND_STR_EDGEREMOVAL_LOOP,"Add the edge weight as internal.",
      NERSTRAND_EDGEREMOVAL_LOOP},
  {NERSTRAND_STR_EDGEREMOVAL_DISTRIBUTE,"Add the edge weight evenly to "
      "remaining edges.",NERSTRAND_EDGEREMOVAL_DISTRIBUTE},
  {NERSTRAND_STR_EDGEREMOVAL_PHANTOM,"Remove the edge but keep track of the "
      "external edge weight.",NERSTRAND_EDGEREMOVAL_PHANTOM}
};

static const cmd_opt_pair_t VERBOSITY_CHOICES[] = {
  {NERSTRAND_STR_VERBOSITY_MINIMUM,"Print minimal information.",
      NERSTRAND_VERBOSITY_MINIMUM},
  {NERSTRAND_STR_VERBOSITY_LOW,"Print only aggregate cluster information.",
      NERSTRAND_VERBOSITY_LOW},
  {NERSTRAND_STR_VERBOSITY_MEDIUM,"Print individual cluster information.",
      NERSTRAND_VERBOSITY_MEDIUM},
  {NERSTRAND_STR_VERBOSITY_HIGH,"Print cluster and progress information.",
      NERSTRAND_VERBOSITY_HIGH},
  {NERSTRAND_STR_VERBOSITY_MAXIMUM,"Print everything.",
      NERSTRAND_VERBOSITY_MAXIMUM}
};


static const cmd_opt_pair_t STOPCONDITION_CHOICES[] = {
  {NERSTRAND_STR_STOPCONDITION_VERTICES,"Stop based on the number of vertices "
      "in the coarse graph.",NERSTRAND_STOPCONDITION_VERTICES},
  {NERSTRAND_STR_STOPCONDITION_EDGES,"Stop based on the number of edges in "
      "the coarse graph.",NERSTRAND_STOPCONDITION_EDGES},
  {NERSTRAND_STR_STOPCONDITION_SIZE,"Stop based on the storage size of the "
      "coarse graph.",NERSTRAND_STOPCONDITION_SIZE}
};


static const cmd_opt_pair_t DISTRIBUTION_CHOICES[] = {
  {NERSTRAND_STR_DISTRIBUTION_BLOCK,"Distribute the vertices in continous " \
      "from the initial ordering.",NERSTRAND_DISTRIBUTION_BLOCK},
  {NERSTRAND_STR_DISTRIBUTION_CYCLIC,"Distribute the vertices in a cyclic " \
      "fashion.",NERSTRAND_DISTRIBUTION_CYCLIC},
  {NERSTRAND_STR_DISTRIBUTION_BLOCKCYCLIC,"Distribute the vertices in a " \
    "blockcyclic fashion.",NERSTRAND_DISTRIBUTION_BLOCKCYCLIC}
};


static const cmd_opt_t OPTS[] = {
  {NERSTRAND_OPTION_HELP,'h',"help","Display this help page.",CMD_OPT_FLAG,
      NULL,0},
  {NERSTRAND_OPTION_NCLUSTERS,'c',"clusters","The number of clusters to "
      "generate. Use '0' for any-way clustering.",CMD_OPT_INT,NULL,0},
  {NERSTRAND_OPTION_SEED,'s',"seed","The random seed to use.",CMD_OPT_INT,NULL,
      0},
  {NERSTRAND_OPTION_NRUNS,'n',"runs","The number of clustering solutions to "
      "generate using succesive random seeds.",CMD_OPT_INT,NULL,0},
  {NERSTRAND_OPTION_NINITSOLUTIONS,'N',"initialsolutions","The number of "
      "initial solutions to generate at the coarsest level.",CMD_OPT_INT,NULL,
      0},
  {NERSTRAND_OPTION_AGGTYPE,'a',"aggregation","The aggregation scheme to use.",
      CMD_OPT_CHOICE, AGGREGATE_CHOICES,__ARRAY_SIZE(AGGREGATE_CHOICES)},
  {NERSTRAND_OPTION_REFTYPE,'r',"refinement","The refinement scheme to use.",
      CMD_OPT_CHOICE, REFINEMENT_CHOICES,__ARRAY_SIZE(REFINEMENT_CHOICES)},
  {NERSTRAND_OPTION_INITYPE,'i',"initialclustering","The initial clustering "
      "scheme to use.",CMD_OPT_CHOICE, INITIAL_CLUSTERING_CHOICES,
      __ARRAY_SIZE(INITIAL_CLUSTERING_CHOICES)},
  {NERSTRAND_OPTION_DEGREE_WEIGHT,'w',"degreeweight","The weight given "
      "to external edges during aggregation.",CMD_OPT_FLOAT,NULL,0},
  {NERSTRAND_OPTION_NREFPASS,'P',"refinementpasses","The maximum number of "
      "refinement passes to use at each step.",CMD_OPT_INT,NULL,0},
  {NERSTRAND_OPTION_MAXREFMOVES,'m',"maxrefinementmoves","The maximum number "
      "of moves a thread can make before synchronizing.",CMD_OPT_INT,NULL,0},
  {NERSTRAND_OPTION_VERBOSITY,'v',"verbosity","The verbosity level of "
      "execution",CMD_OPT_CHOICE,VERBOSITY_CHOICES,
      __ARRAY_SIZE(VERBOSITY_CHOICES)},
  {NERSTRAND_OPTION_MODSTATS,'M',"modularitystatistics","Print modularity "
      "information from multiple clusterings (should be used with '-n').",
      CMD_OPT_FLAG,NULL,0},
  {NERSTRAND_OPTION_ICSTATS,'I',"initialclusteringstatistics","Print initial "
      "clustering information from multiple clusterings.",CMD_OPT_FLAG,NULL,0},
  {NERSTRAND_OPTION_LBSTATS,'L',"loadbalancestatistics","Print load balancing "
      "information.",CMD_OPT_FLAG,NULL,0},
  {NERSTRAND_OPTION_AGGSTATS,'A',"aggregationstatistics","Print aggregation "
      "information.",CMD_OPT_FLAG,NULL,0},
  {NERSTRAND_OPTION_REFSTATS,'R',"refinementstatistics","Print refinement "
      "information.",CMD_OPT_FLAG,NULL,0},
  {NERSTRAND_OPTION_TIME,'t',"times","Print timing information",CMD_OPT_FLAG,
      NULL,0},
  {NERSTRAND_OPTION_AGG_RATE,'g',"aggregationrate","Rate at which to limit "
      "the aggregation.",CMD_OPT_FLOAT,NULL,0},
  {NERSTRAND_OPTION_SPATYPE,'S',"sparsify","Method of sparsification to use "
      "if any.",CMD_OPT_CHOICE,SPARSIFY_CHOICES,
      __ARRAY_SIZE(SPARSIFY_CHOICES)},
  {NERSTRAND_OPTION_DISTYPE,'E',"edgeremoval","Method of redistributing the "
      "weight of removed edges.",CMD_OPT_CHOICE,EDGEREMOVAL_CHOICES,
      __ARRAY_SIZE(EDGEREMOVAL_CHOICES)},
  {NERSTRAND_OPTION_SUPERNODE_RATIO,'u',"supernoderatio","An alpha value used "
      "for classifying supernodes (larger means more nodes classified as "
      "supernodes, 0 means none will be classified).",
      CMD_OPT_FLOAT,NULL,0},
  {NERSTRAND_OPTION_STOPRATIO,'o',"stopratio","The vertex ratio between "
      "successive graphs that will cause coarsening to terminate\n",
      CMD_OPT_FLOAT,NULL,0},
  {NERSTRAND_OPTION_STOPCONDITION,'O',"stopcondition","The stopping condition "
      "for coarsening",CMD_OPT_CHOICE,STOPCONDITION_CHOICES,
      __ARRAY_SIZE(STOPCONDITION_CHOICES)},
  {NERSTRAND_OPTION_CNVTXS_PER_CLUSTER,'V',"cnvtxspercluster","The number of "
      "coarse vertices per cluster when making initial clusters",CMD_OPT_INT,
      NULL,0},
  {NERSTRAND_OPTION_BLOCKSIZE,'B',"blocksize","The number of vertices "
      "in block for block-cyclic distribution (no effect when used with "
       "other distributions.",CMD_OPT_INT,NULL,0},
  {NERSTRAND_OPTION_DISTRIBUTION,'D',"distribution","The distribution to use "
      "for assigning vertices to threads.",CMD_OPT_CHOICE,DISTRIBUTION_CHOICES,
      __ARRAY_SIZE(DISTRIBUTION_CHOICES)},
  {NERSTRAND_OPTION_NTHREADS,'T',"threads","The number of threads to use for "
      "clustering",CMD_OPT_INT,NULL,0},
  {NERSTRAND_OPTION_IGNORE_WEIGHT,'W',"ignoreweights","Ignore edge weights "
      "in input file (default=false).",CMD_OPT_FLAG,NULL,0},
  {NERSTRAND_OPTION_RESTEP,'j',"restep","The fraction of edges that will "
      "trigger a re-step.",CMD_OPT_FLOAT,NULL,0}
};

static const size_t NOPTS = __ARRAY_SIZE(OPTS);

#undef __ARRAY_SIZE




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static int __usage(
    const char * const name)
{
  fprintf(stderr,"USAGE:\n");
  fprintf(stderr,"%s [options] <graphfile> [ <clusterfile> | - ]\n",name);
  fprintf(stderr,"\n");
  fprintf(stderr,"Options:\n");
  fprint_cmd_opts(stderr,OPTS,NOPTS);
  return 1;
}


static double * __parse_args(
    cmd_arg_t * args, 
    size_t nargs,
    const char ** r_input, 
    const char ** r_output)
{
  size_t i, xarg;
  double * options = NULL;
  const char * input_file = NULL, * output_file = NULL;

  options = nerstrand_init_options();

  /* set default verbosity to low */
  options[NERSTRAND_OPTION_VERBOSITY] = NERSTRAND_VERBOSITY_LOW;

  for (i=0;i<nargs;++i) {
    switch (args[i].type) {
      case CMD_OPT_CHOICE:
        options[args[i].id] = (double)args[i].val.o;
        break;
      case CMD_OPT_BOOL:
        options[args[i].id] = (double)args[i].val.b;
        break;
      case CMD_OPT_INT:
        options[args[i].id] = (double)args[i].val.i;
        break;
      case CMD_OPT_FLOAT:
        options[args[i].id] = (double)args[i].val.f;
        break;
      case CMD_OPT_FLAG:
        options[args[i].id] = 1.0;
        break;
      default:
        break;
    }
  }

  xarg = 0;
  for (i=0;i<nargs;++i) {
    /* check for help */
    if (args[i].id == NERSTRAND_OPTION_HELP) {
      goto CLEANUP;
    }
    if (args[i].type == CMD_OPT_XARG) {
      switch (xarg) {
        case 0 :
          input_file = args[i].val.s;
          break;
        case 1 :
          output_file = args[i].val.s;
          if (strcmp(output_file,"-") == 0) {
            output_file = NULL;
          }
          break;
        default :
          eprintf("Unknown extra argument '%s'\n",args[i].val.s);
          goto CLEANUP;
      }
      ++xarg;
    }
  }

  if (input_file == NULL) {
    eprintf("Must supply at least an input graph to cluster\n");
    goto CLEANUP;
  }

  *r_output = output_file;
  *r_input = input_file;

  return options;

  CLEANUP:
  dl_free(options);
  *r_output = NULL;
  *r_input = NULL;

  return NULL;
}




/******************************************************************************
* MAIN ************************************************************************
******************************************************************************/


clustering_t * nerstrand_internal_cluster(
    objective_t * objective, 
    const graph_t * graph);


int main(
    const int argc, 
    char ** const argv)
{
  int rv, err;
  size_t nargs;
  adj_t e;
  vtx_t nvtxs, i;
  double * options = NULL;
  cid_t * uwhere = NULL;
  adj_t * xadj = NULL;
  vtx_t * adjncy = NULL;
  wgt_t * adjwgt = NULL;
  wgt_t * iadjwgt = NULL;
  const char * input_file = NULL;
  const char * output_file = NULL;
  cmd_arg_t * args = NULL;
  objective_t * objective = NULL;
  graph_t * graph = NULL;
  clustering_t * clustering = NULL;

  vtx_t * mynvtxs, ** myadjncy;
  adj_t ** myxadj;
  wgt_t ** myadjwgt;

  rv = cmd_parse_args(argc-1,argv+1,OPTS,NOPTS,&args,&nargs);
  if (rv != DL_CMDLINE_SUCCESS) {
    eprintf("Failed with %d\n",rv);
    __usage(argv[0]);
    rv = 1;
    goto CLEANUP;
  }

  options = __parse_args(args,nargs,&input_file,&output_file);

  if (options == NULL) {
    __usage(argv[0]);
    rv = 1;
    goto CLEANUP;
  }

  if (parse_objective(options,&objective) != NERSTRAND_SUCCESS) {
    __usage(argv[0]);
    rv = 1;
    goto CLEANUP;
  }
  dl_free(options);
  options = NULL;

  if (objective->time) {
    dl_start_timer(&(objective->timers.total));
    dl_start_timer(&(objective->timers.io));
  }


  vprintf(objective->verbosity,NERSTRAND_VERBOSITY_LOW,"Reading '%s'\n", \
      input_file);
  err = wildriver_read_graph(input_file, &nvtxs, NULL, NULL, NULL, &xadj, \
      &adjncy, NULL, &adjwgt);

  if (err != 1) {
    eprintf("Error reading graph\n");
    rv = 1;
    goto CLEANUP;
  }

  if (objective->ignore_weight && adjwgt) {
    // replace weight with 1's
    for (e=0;e<xadj[nvtxs];++e) {
      adjwgt[e] = 1;
    }
  }

  if (objective->nthreads > 1) {
    graph = distribute_graph(nvtxs,xadj,adjncy,adjwgt,iadjwgt,NULL,
        objective->nthreads,objective->distribution,objective->blocksize);

    /* free the single arrays */
    dl_free(xadj);
    dl_free(adjncy);
    if (adjwgt) {
      dl_free(adjwgt);
    }
    if (iadjwgt) {
      dl_free(iadjwgt);
    }
  } else {
    /* just use single arrays */
    mynvtxs = vtx_alloc(1);
    mynvtxs[0] = nvtxs;

    myxadj = r_adj_alloc(1);
    myxadj[0] = xadj;

    myadjncy = r_vtx_alloc(1);
    myadjncy[0] = adjncy;

    if (adjwgt) {
      myadjwgt = r_wgt_alloc(1);
      myadjwgt[0] = adjwgt;
    } else {
      myadjwgt = NULL;
    }

    graph = setup_graph(mynvtxs,myxadj,myadjncy,myadjwgt,NULL,NULL,NULL,1);
    graph->free_xadj = 1;
    graph->free_adjncy = 1;
    if (myadjwgt) {
      graph->free_adjwgt = 1;
    }
  }

  vprintf(objective->verbosity,NERSTRAND_VERBOSITY_LOW, \
      "Read '%s' with "PF_VTX_T" vertices and "PF_ADJ_T" edges\n", \
      input_file,graph->nvtxs,graph->nedges);

  if (objective->time) {
    dl_stop_timer(&(objective->timers.io));
    dl_start_timer(&(objective->timers.clustering));
  }

  clustering = nerstrand_internal_cluster(objective,graph);

  if (objective->time) {
    dl_stop_timer(&(objective->timers.clustering));
    dl_start_timer(&(objective->timers.io));
  }

  print_clustering(clustering,graph,0,objective->verbosity);

  if (output_file) {
    uwhere = cid_alloc(nvtxs);
    unify_where((const cid_t * const *)clustering->where,graph,uwhere);

    // write labels
    file_t * file;
    if (dl_open_file(output_file, "w", &file) != DL_FILE_SUCCESS) {
      // output error and continue cleaning up
      eprintf("Unable to open file: '%s'\n", output_file);
    } else {
      for (i=0;i<nvtxs;++i) {
        dl_fprintf(file, PF_CID_T"\n", uwhere[i]);  
      }
      dl_close_file(file);
    }

    dl_free(uwhere);
  }

  if (objective->time) {
    dl_stop_timer(&(objective->timers.io));
    dl_stop_timer(&(objective->timers.total));

    print_timers(objective);
  }

  CLEANUP:

  if (args) {
    dl_free(args);
  }
  if (objective) {
    free_objective(objective);
  }
  if (clustering) {
    free_clustering(clustering);
  }
  if (graph) {
    free_graph(graph);
  }
  if (options) {
    dl_free(options);
  }

  return rv;
}




#endif
