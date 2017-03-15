/* mcmain.c
 *    Parse command line options and dispatch required functions.
 *
 * Copyright (C) 2001-2006 The Regents of the University of California
 * by Glenn Tesler
 *
 * Contains code from main.c in GRAPPA 1.02
 * Copyright (C) 2000-2001  The University of New Mexico and
 *                          The University of Texas at Austin
 * by David A. Bader, Bernard M.E. Moret, Tandy Warnow, Stacia K Wyman, Mi Yan
 *
 * See file COPYRIGHT for details.
 *****************************************************************************
 * This file is part of GRIMM.
 *
 * GRIMM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License, Version 2,
 * dated June 1991, as published by the Free Software Foundation.
 *
 * GRIMM is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

/* Last modified on Wed Aug 2, 2006, by Glenn Tesler
 */

/* Contains excerpts from GRAPPA 1.02: main.c */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <stddef.h>   /* for offsetof */
#include <unistd.h>
#include <string.h>

#include "mcstructs.h"
#include "uniinvdist.h"
#include "mcrdist.h"
#include "mcread_input.h"
#include "write_data.h"
#include "scenario.h"
#include "testrev.h"
#include "opt_scenario.h"
#include "unsigned.h"
#include "e_malloc.h"
#include "circ_align.h"
#include "texgraph.h"
#include "ext_function.h"
#include "countperms.h"
#include "unsignedhc.h"
#include "matrixmisc.h"

/* put quotes around version number */
#define qstr1(s) #s
#define qstr(s) qstr1(s)
#define QVERS qstr(VERS)

FILE *outfile;


void do_fn_distmatrix(int NUM_GENES, int NUM_CHROMOSOMES, int NUM_GENOMES,
		      int circular,
		      int unsigned_dist,
		      int verbose,
                      struct genome_struct *genome_list);
void do_fn_distmatrix_v(int NUM_GENES, int NUM_CHROMOSOMES, int NUM_GENOMES,
			int circular,
			int unsigned_dist,
			struct genome_struct *genome_list);



void print_usage(char *progname) {
  fprintf(stderr,"\nGRIMM %s\n",QVERS);
  fprintf(stderr,"Copyright (C) 2001-2006 The Regents of the University of California\n");
  fprintf(stderr,"Contains code from GRAPPA 1.02 (C) 2000-2001 The University of New Mexico\n");
  fprintf(stderr,"                                         and The University of Texas at Austin\n");
fprintf(stderr,"See file COPYRIGHT for full copyright and authorship details.\n");

  fprintf(stderr,
	  "\nUsage: %s -f datafile [-o outfile] [other options]",
	  progname);

  fprintf(stderr,
	  "\n\n");

  fprintf(stderr,"   -f datafile: name of file containing all the genomes\n");
  fprintf(stderr,"   -o outfile: output file name.  Output goes to STDOUT unless -o used.\n");
  fprintf(stderr,"   -v: Verbose output.\n");
  fprintf(stderr,"       In pairwise mode, outputs all the graph parameters.\n");
  fprintf(stderr,"       In matrix mode, outputs matrices for all graph parameters.\n");
  fprintf(stderr,"       In special functions, outputs extra information.\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"Genome type:\n");
  fprintf(stderr,"   -L: unichromosomal (directed) linear\n");
  fprintf(stderr,"   -C: unichromosomal circular\n");
  fprintf(stderr,"   Multichromosomal (undirected linear chromosomes) used unless -L/-C given\n");
  /* -L: "0" chromosomes -- linear -- the old unichromosomal case */
  /* -C: "0" chromosomes -- circular */
  /* -D2: pad it w/caps, but then do unichromosomal alg only on it;
     This has no mathematical usefulness whatsoever, it is just
     useful for debugging. */
  /* -M #: use user specified caps -- also just useful for debugging */
  fprintf(stderr,"\n");
  fprintf(stderr,"Genome selection:\n");
  fprintf(stderr,"   -g i,j: compare genomes i and j (counting from 1)\n");
  fprintf(stderr,"      Rearrangement scenario will go from genome i to genome j.\n");
  fprintf(stderr,"      In Hannenhalli-Pevzner papers, genome i is gamma and genome j is pi.\n");
  fprintf(stderr,"      For files with 3 or more genomes, unless -g is specified, a distance\n");
  fprintf(stderr,"      matrix is printed and options -d, -c, -z, -s are not available.\n");
  fprintf(stderr,"   -m: force matrix output even for 2 genomes\n");

  fprintf(stderr,"\n");

  fprintf(stderr,"Function: (default -d -c -s for multichromosomal, -d -s for unichromosomal)\n");
  fprintf(stderr,"   -d: distance\n");
  fprintf(stderr,"   -c: show capping\n");

  /* TODO: choose something better than -z */
  fprintf(stderr,"   -z: show capping and chromosome delimeters too\n");


  fprintf(stderr,"   -s: display an optimal scenario (see -S below)\n");
  fprintf(stderr,"   -u: distance of unsigned genomes; find signs giving optimal distance\n");
  fprintf(stderr,"   -U #: unsigned genomes approximation algorithm with # iterations\n");
  fprintf(stderr,"   -W file: optional file with weight matrix for -U\n");
  fprintf(stderr,"   -S #: produce/display optimal scenario in a certain format (#=1,...,7)\n");
  fprintf(stderr,"      #=1: operations may occur in any order\n");
  fprintf(stderr,"      #=2,...,7: operations are prioritized\n");
  fprintf(stderr,"        rev. short->long; flip chromos; fusion; transloc short->long; fission\n");
  fprintf(stderr,"      #=1,2,3,4: scenario respects initial capping\n");
  fprintf(stderr,"      #=5,6,7: scenario may not respect initial capping, so can avoid flips\n");
  fprintf(stderr,"      #=2: 1-line permutations with caps\n");
  fprintf(stderr,"      #=3: multiline permutations, show caps, / at end of chromo\n");
  fprintf(stderr,"      #=4: multiline permutations, hide caps, $ at end of chromo\n");
  fprintf(stderr,"      #=5,6,7: different greedy approaches to reversals first\n");
  fprintf(stderr,"        They allow recapping, so hide caps and show $ at end of each chromo\n");

  fprintf(stderr,"\nSpecial functions:\n");
  fprintf(stderr,"   -t#: test all 1-step reversals, list the ones reducing distance by 1,\n");
  fprintf(stderr,"       and give statistics\n");
  fprintf(stderr,"       #=1,2,3,4 specifies test method (debugging)\n");
  fprintf(stderr,"   -F n: tabulate statistics for given # genes n\n");
  fprintf(stderr,"   -T xxx: output graph in TeX form.  xxx is a combination of\n");
  fprintf(stderr,"      l: left-to-right format  c: component/cycle format  (only use one)\n");
  fprintf(stderr,"      b: bracket chromosomes\n");
  fprintf(stderr,"      s: signed labels\n");
  fprintf(stderr,"      u: unsigned labels 2x-1,2x  a: unsigned labels a,b  (only use one)\n");
  fprintf(stderr,"      r: list of raw components, cycles, paths, not TeX format\n");
  fprintf(stderr,"      1-6: stage in multichromosomal capping algorithm\n");
  fprintf(stderr,"      Some argument must be given, so use 1 if nothing else\n");
  fprintf(stderr,"\nDebugging (not for general use):\n");
  fprintf(stderr,"   -D1: delete caps from multichromosomal genomes, treat as linear\n");
  fprintf(stderr,"   -D2: put caps on multichromosomal genomes, then treat them as unichromosomal\n");
  fprintf(stderr,"   -M n: the data appears in unichromosomal format, but has been manually\n");
  fprintf(stderr,"       capped for n chromosomes\n");
  fprintf(stderr,"   -X arg: demo of extensions\n");

  fprintf(stderr,"\n");
  fprintf(stderr,"Format of datafile: each entry has the form\n");
  fprintf(stderr,"   >human\n");
  fprintf(stderr,"   # Chromosome 1\n");
  fprintf(stderr,"   5 -2 3 -1 $\n");
  fprintf(stderr,"   # Chromosome 2\n");
  fprintf(stderr,"   7 -4 6 $\n");
  fprintf(stderr,"   # Chromosome 3\n");
  fprintf(stderr,"   8\n");
  fprintf(stderr,"   ...\n");
  fprintf(stderr,"   >mouse\n");
  fprintf(stderr,"   ...\n");
  fprintf(stderr,"Each genome starts with a line \">genomename\".\n");
  fprintf(stderr,"The genes are separated by any whitespace, and may be on multiple lines.\n");
  fprintf(stderr,"The chromosomes are separated by \";\" or \"$\".  New lines between\n");
  fprintf(stderr,"chromosomes are recommended for neatness but are not required.\n");
  fprintf(stderr,"Comments begin with \"#\" and go to the end of the line.\n");
  fprintf(stderr,"\nWith the -U option only, strips of unknown orientation may be bracketed:\n");
  fprintf(stderr,"   5 [ -2 3 ] -1 $\n");
  fprintf(stderr,"means it's either -2 3 or -3 2 but the unbracketed genes are fixed.\n");
  fprintf(stderr,"If no [] with the -U option then unknown signs on all genes in genomes 2,3,...\n");
  fprintf(stderr,"Brackets {} () are reserved and currently ignored.\n");
  fprintf(stderr,"\n");
}



/* TODO: This code is monstrous.  Split it up. */
int main(int argc, char *argv[]) {
  /* parsing */
  extern char *optarg;
  extern int optopt;
  int errflg;
  int c;


  /* file i/o */
  FILE *input        = (FILE *) NULL;
  char *inputfname   = (char *) NULL;
  char *outputfname  = (char *) NULL;
  struct stat statbuf;

  /* statistics about the data */
  int NUM_GENES;   /* number of genes in input */
  int NUM_GENOMES; /* number of gene orders in input */
  int NUM_CHROMOSOMES; /* number of chromosomes per genome */

  /* options */

  int debugcaps          = 0;      /* 0: normal;
				      1: delete caps, treat as unichrom
				      2: keep caps, treat as unichrom
				      */
  int unichromosome      = FALSE;  /* TRUE: classical rev dist, lin/circ
				      FALSE: multichromosomal algorithm */

  int circular           = FALSE;  /* unichrom circular genome */

  int unsigned_dist      = FALSE;  /* genomes signed or unsigned? */
  int verbose            = FALSE;
  int precapped          = -1;     /* data already has caps for
				      this # of chromosomes
				      (-1 means not precapped) */

  int fn_distance        = FALSE;  /* function: compute distance */
  int fn_capping         = 0;      /* function: compute capping
				    * 0: don't
				    * 1: display capping on one line
				    * 2: display capping on multiple lines,
				    *    one per chromosome, along with
				    *    chromosome delimeters.
				    */
  int fn_scenario        = FALSE;  /* function: compute optimal scenario */
  int scenario_type      = 4;      /*     method by which to compute it */

  int fn_revstats        = 0;      /* function: test all one step reversals */
  int fn_any             = FALSE;  /* no functions specified */

  int fn_matrix          = FALSE;  /* function: compute distance matrix */
                                   /* Activated by having >2 genomes but
				      not using -g to specify which */
  int fn_texgraph        = FALSE;  /* function: draw graph in TeX */
  char *texgraph_args    = NULL;

  /* code development only */
  int fn_external        = FALSE;  /* function: some external function */
  char *ext_args         = NULL;


  int fn_genfn           = FALSE;  /* function: compute generating function */
  char *fn_genfn_args    = NULL;


  int fn_unsignedhc      = FALSE;  /* approx unsigned */
  int fn_unsignedhc_numit = 0;     /* # iterations to do */
  int nsegs              = 0;      /* # identified segments */
  int *seg_list          = (int *) 0; /* table of segments */
  FILE *weight_file      = (FILE *) NULL;
  char *weight_fname     = (char *) NULL;
  int **weight_matrix    = (int **) NULL;
  int weight_nr          = 0;
  int weight_nc          = 0;



  int gindex1            = -1;     /* which genomes to compare? */
  int gindex2            = -1;     /* which genomes to compare? */
 





  struct genome_struct *genome_list;
  graphstats_t graphstats;         /* parameter details, for verbose option */

  int i;
  int dist;                        /* computed distance */
#ifdef DEBUGc
  int dist2;                       /* unichrom dist of capped perms */
#endif
  int dist_error;                  /* boolean: was distance approximated? */


  outfile            = stdout;



  errflg = 0;
  while ((c = getopt(argc,argv,"f:o:CLdczsmt:uvM:D:g:T:S:X:F:U:W:")) != -1) {
    switch (c) {
    case 'f':
      inputfname = optarg;
      break;
    case 'o':
      outputfname = optarg;
      break;
    case 'L':
      unichromosome = TRUE;
      circular = FALSE;
      break;
    case 'C':
      unichromosome = TRUE;
      circular = TRUE;
      break;

    case 'D':  /* debug caps */
      if (sscanf(optarg,"%d", &debugcaps) != 1 ||
	  debugcaps < 1 || debugcaps > 2) {
	fprintf(stderr, "ERROR: -D #   must be 1 or 2\n");
	errflg++;
      }
      unichromosome = TRUE;
      break;
    case 'M':
      if (sscanf(optarg,"%d",&precapped) != 1 || precapped < 0) {
	fprintf(stderr,"# chromosomes must be at least 0\n");
	errflg++;
      }
      break;

    case 'd':
      fn_distance = TRUE;
      fn_any = TRUE;
      break;
    case 'c': /* display capping on 1 line */
      fn_capping = 1;
      fn_any = TRUE;
      break;
    case 'z': /* display capping on multiple lines w/chrom delims */
      fn_capping = 2;
      fn_any = TRUE;
      break;
    case 's':
      fn_scenario = TRUE;
      fn_any = TRUE;
      break;
    case 'S':
      fn_scenario = TRUE;
      fn_any = TRUE;

      if (sscanf(optarg,"%d",&scenario_type) != 1
	  || scenario_type <= 0
	  || scenario_type > 7
	  ) {
	fprintf(stderr,"invalid option for -S\n");
	errflg++;
      }


      break;
    case 'm':
      fn_matrix = TRUE;
      break;
    case 'X':
      fn_external = TRUE;
      ext_args = e_malloc((strlen(optarg)+1) * sizeof(char),"ext_args");
      strcpy(ext_args,optarg);
      break;
    case 'F':
      fn_genfn = TRUE;
      fn_genfn_args = optarg;
      break;
#if 0
      if (sscanf(optarg,"%d",&fn_genfn_n) != 1 ||
	  fn_genfn_n < 1) {
	fprintf(stderr,"ERROR: -F #   must be at least 1\n");
	errflg++;
      }
      fn_any = TRUE;
      break;
#endif

    case 't':
      if (sscanf(optarg,"%d",&fn_revstats) != 1 ||
	  fn_revstats < 1 || fn_revstats > 4) {
	fprintf(stderr,"ERROR: -t #   must be between 1 and 4\n");
	errflg++;
      }
      fn_any = TRUE;
      break;

    case 'u':
      unsigned_dist = TRUE;
      break;

    case 'v':
      verbose = TRUE;
      break;

    case 'g':
      if (sscanf(optarg,"%d,%d",&gindex1,&gindex2) != 2
	  || gindex1<=0 || gindex2<=0) {
	fprintf(stderr,"-g needs two genome numbers\n");
	errflg++;
      } else {
	/* internally, genomes counted starting at 0 */
	gindex1--; gindex2--;
      }
      break;


      /* TeX graph */
    case 'T':
      fn_texgraph = fn_any = TRUE;
      texgraph_args = optarg;
      break;

    case 'U':
      fn_unsignedhc = fn_any = TRUE;
      sscanf(optarg,"%d",&fn_unsignedhc_numit);
      if (fn_unsignedhc_numit <= 0) {
	fprintf(stderr, "Option -U #  requires #>0\n");
	errflg++;
      }
      break;

    case 'W':
      weight_fname = optarg;
      break;

    case ':':
      fprintf(stderr, "Option -%c requires an operand\n", optopt);
      errflg++;
      break;
    case '?':
      fprintf(stderr, "Unrecognized option: -%c\n", optopt);
      errflg++;
      break;
    default:
      fprintf(stderr, "Option -%c not yet implemented\n", optopt);
      break;
    }
  }


  /* check we have usable parameters */
  if (errflg
      ) {
    print_usage(argv[0]);
    exit(-1);
  }


  /* This function doesn't read any files */
  if (fn_genfn) {
    /*    countperms(fn_genfn_n, circular, unsigned_dist); */
    countperms(fn_genfn_args, circular, unsigned_dist);
    return(0);
  }


  /* continue error checking */

  if (inputfname == (char *)NULL) {
    fprintf(stderr,"ERROR: input filename required\n");
    print_usage(argv[0]);
    exit(-1);
  }
  else {
    if (stat(inputfname, &statbuf)) {
      fprintf(stderr,"ERROR: Input file (%s): ",inputfname);
      perror("");
      print_usage(argv[0]);
      exit(-1);
    }
    if (statbuf.st_mode & S_IFDIR) {
      fprintf(stderr,"ERROR: Input file (%s) is a directory.\n",inputfname);
      print_usage(argv[0]);
      exit(-1);
    }
    input = fopen(inputfname,"r");
    if (input == (FILE*)NULL) {
      fprintf(stderr,"ERROR: Could not open input file (%s): ",inputfname);
      perror("");
      print_usage(argv[0]);
      exit(-1);
    }
  }


  if (outputfname != (char *)NULL) {
    if (!stat(outputfname, &statbuf)) {
      fprintf(stderr,"ERROR: Output file (%s) already exists.\n",outputfname);
      print_usage(argv[0]);
      exit(-1);
    }
    outfile = fopen(outputfname,"a+");
    if (outfile == (FILE *)NULL) {
      fprintf(stderr,"ERROR: Could not open output file (%s): ",outputfname);
      perror("");
      print_usage(argv[0]);
      exit(-1);
    }
  }


  /***********************************************************************/
  /*                       read in the genome data                       */
  /***********************************************************************/
  mcread_data(input,&genome_list,&NUM_GENES,&NUM_GENOMES,&NUM_CHROMOSOMES,
	      &nsegs,&seg_list);

  if (NUM_GENES == 0 || NUM_GENOMES == 0) {
    fprintf(stderr,"No data found in file.\n");
    print_usage(argv[0]);
    exit(-1);
  }


  /***********************************************************************/
  /*                   determine which genomes to compare                */
  /***********************************************************************/

  if (fn_matrix && gindex1 != -1) {
    fprintf(stderr,"The -g and -m options may not be used together.\n");
    print_usage(argv[0]);
    exit(-1);
  }

#if 0
  if (fn_unsignedhc && debugcaps>0) {
    fprintf(stderr,"The -U and -D options may not be used together.\n");
    print_usage(argv[0]);
    exit(-1);
  }
#endif

  if (fn_unsignedhc) {
  } else
  if ((NUM_GENOMES>2 && gindex1 == -1) || fn_matrix) {
    if (fn_any) {
      fprintf(stderr,"A distance matrix is displayed when there are multiple genomes or\n");
      fprintf(stderr,"the -m option is given.  The -d -c -s -t -T options are not available;\n");
      fprintf(stderr,"you must specify two genomes with the -g #,# option.\n");
    }
    fn_scenario = fn_distance = fn_capping = FALSE;
    fn_texgraph = FALSE;
    fn_revstats = 0;
    fn_any = fn_matrix = TRUE;
  }


#if 0
  if (NUM_GENOMES != 2 && gindex1 == -1 && gindex2 == -1) {
    fprintf(stderr,"ERROR: %d genomes, but require exactly 2\n",NUM_GENOMES);
    exit(-1);
  }
#endif

  if (gindex1 == -1 && !fn_unsignedhc) {
    gindex1 = 0; gindex2 = 1;
  }
  if (gindex1 >= NUM_GENOMES || gindex2 >= NUM_GENOMES) {
    fprintf(stderr,
	    "ERROR: specified %d genomes, but requested -g %d,%d\n",
	    NUM_GENOMES,gindex1+1,gindex2+1);
    print_usage(argv[0]);
    exit(-1);
  }
  /* 0 or 1 genomes will be trapped by the above */



  /***********************************************************************/
  /*                   adjust data if pre-capped                         */
  /***********************************************************************/
 

  /* data has caps on it already */
  if (precapped>=0) {
    if (NUM_CHROMOSOMES != 1) {
      fprintf(stderr,"ERROR: -M indicates data is capped with manual cap numbers,\n");
      fprintf(stderr,"ERROR: so it shouldn't have chromosome delimeters.\n");
      print_usage(argv[0]);
      exit(-1);
    }

    for (i=0 ; i<NUM_GENOMES ; i++) {
      delete_caps(&genome_list[i], NUM_GENES, NUM_CHROMOSOMES);
      delete_caps_segs(nsegs, seg_list);
    }

    NUM_GENES -= 2*NUM_CHROMOSOMES;

    NUM_CHROMOSOMES = precapped;

    if (NUM_CHROMOSOMES == 0) unichromosome = TRUE;
    else
      for (i=0 ; i<NUM_GENOMES ; i++) {
	if (!check_caps(&genome_list[i], NUM_GENES, NUM_CHROMOSOMES)) {
	  fprintf(stderr,"ERROR: genome %d is not capped properly.\n", i+1);
	  print_usage(argv[0]);
	  exit(-1);
	}
      }
  }


  /***********************************************************************/
  /*                   adjust data if unichromosomal                     */
  /***********************************************************************/

  /* -D2 option: keep caps in place, treat as unichromosomal */
  if (debugcaps == 2) {
    unichromosome = TRUE;
    NUM_CHROMOSOMES = 0;
  }

  /* if requested circular or linear but have >1 chromosomes, error */
  if (unichromosome && NUM_CHROMOSOMES>1 && debugcaps == 0) {
    fprintf(stderr,"ERROR: -C and -L are only for unichromosomal genomes, but have %d chromosomes\n",NUM_CHROMOSOMES);
    print_usage(argv[0]);
    exit(-1);
  }


  /* -C and -L options require deleting caps */
  if (unichromosome) {
    for (i=0 ; i<NUM_GENOMES ; i++) {
      delete_caps(&genome_list[i], NUM_GENES, NUM_CHROMOSOMES);
      delete_caps_segs(nsegs, seg_list);
    }

    NUM_GENES -= 2*NUM_CHROMOSOMES;
    NUM_CHROMOSOMES = 0;
    unichromosome = TRUE;
  }


  /***********************************************************************/
  /*                   adjust data if circular                           */
  /***********************************************************************/

  if (circular) {
    circular_align(genome_list, NUM_GENOMES, NUM_GENES,
		   0,TRUE);           /* align to 1st gene of 1st genome */

    /* now the linear inversion distance equals the circular inversion dist */
  }




  /***********************************************************************/
  /*                   display report                                    */
  /***********************************************************************/

 if (!fn_texgraph) {
   fprintf(outfile, "Running:                      \t%s\n", argv[0]);
   fprintf(outfile, "Input:                        \t%s\n",inputfname);
   fprintf(outfile, "Genome:                       \t%s\n",
	   !unichromosome ? "Multichromosomal" :
	   circular ? "Circular" : "Linear");
   fprintf(outfile, "Signs:                        \t%s\n",
	   fn_unsignedhc ? (nsegs==0
			    ? "Unsigned (use approx algorithm)"
			    : "Partially signed (use approx algorithm)")
	   :
	   (unsigned_dist ? "Unsigned" : "Signed"));

  fprintf(outfile,  "Number of Genomes:            \t%d\n",NUM_GENOMES);

  if (!unichromosome)
    fprintf(outfile,"Number of Genes:              \t%d + %d caps\n",
	    NUM_GENES-2*NUM_CHROMOSOMES, 2*NUM_CHROMOSOMES);
  else
    fprintf(outfile,"Number of Genes:              \t%d\n",NUM_GENES);


  fprintf(outfile,  "Number of Chromosomes:        \t%d",NUM_CHROMOSOMES);
  if (unichromosome)
    fprintf(outfile, " (unichromosomal)\n");
  else if (NUM_CHROMOSOMES == 1)
    fprintf(outfile, " (multichromosomal)\n");
  else
    fprintf(outfile, "\n");

  fflush(outfile);
 }


#if 0
#ifdef DEBUG
  fprintf(outfile,"options: distance=%d capping=%d scenario=%d\n",
	  fn_distance,fn_capping,fn_scenario);
#endif
#endif




  /***********************************************************************/
  /*                  do the requested functions                         */
  /***********************************************************************/

  if (fn_external)
    do_fn_external(NUM_GENES, NUM_CHROMOSOMES, NUM_GENOMES,
		   circular,
		   unsigned_dist,
		   verbose,
		   genome_list,
		   ext_args);
  else
    if (fn_unsignedhc) {
    fprintf(outfile,"\n======================================================================\n\n");


    /* load weight matrix */

    if (weight_fname != (char *) NULL) {
      if (stat(weight_fname, &statbuf)) {
	fprintf(stderr,"ERROR: weights file (%s): ",weight_fname);
	perror("");
	print_usage(argv[0]);
	exit(-1);
      }
      if (statbuf.st_mode & S_IFDIR) {
	fprintf(stderr,"ERROR: weights file (%s) is a directory.\n",weight_fname);
	print_usage(argv[0]);
	exit(-1);
      }
      weight_file = fopen(weight_fname,"r");
      if (weight_file == (FILE*)NULL) {
	fprintf(stderr,"ERROR: Could not open weights file (%s): ",weight_fname);
	perror("");
	print_usage(argv[0]);
	exit(-1);
      }

      load_matrix(weight_file,
		  &weight_nr,
		  &weight_nc,
		  &weight_matrix);
      if (weight_nr != NUM_GENOMES || weight_nc != NUM_GENOMES) {
	fprintf(stderr, "ERROR: weight matrix is %d by %d, but should be %d by %d\n",
		weight_nr, weight_nc, NUM_GENOMES, NUM_GENOMES);
	errflg++;
      }
      if (!check_sym_mat_int2d(weight_matrix, weight_nr, weight_nc)) {
	fprintf(stderr, "ERROR: weight matrix is not symmetric");
	errflg++;
      }
      if (errflg) {
	print_usage(argv[0]);
	exit(-1);
      }
    }


    /* TODO: decide which parts of -U should be output with vs. w/o -U */
    verbose = TRUE;
    if (gindex1 >= 0) {
      /* -g i1,i2 version */
      unsignedhc_mcdist(NUM_GENES, NUM_CHROMOSOMES,
			circular, verbose,
			genome_list,
			gindex1, gindex2,
			(nsegs==0) ? -1 : nsegs,
			seg_list,
			weight_matrix,
			fn_unsignedhc_numit);
    } else {
      /* multiple genome version */
      unsignedhc_mcdist_mat(NUM_GENES, NUM_CHROMOSOMES, NUM_GENOMES,
			    circular, verbose,
			    genome_list,
			    (nsegs==0) ? -1 : nsegs,
			    seg_list,
			    weight_matrix,
			    fn_unsignedhc_numit);
    }
  }
  else if (fn_matrix)
     do_fn_distmatrix(NUM_GENES, NUM_CHROMOSOMES, NUM_GENOMES,
		      circular,
		      unsigned_dist,
		      verbose,
                      genome_list);
  else if (fn_texgraph)
    draw_texgraph(&genome_list[gindex1], &genome_list[gindex2],
		  NUM_GENES, NUM_CHROMOSOMES,
		  texgraph_args);
 else {





  /* no functions supplied. 
     Default: distance, capping (if applicable), and scenario  */
  if (!fn_any) {
    fn_distance = TRUE;
    if (NUM_CHROMOSOMES>0) fn_capping = TRUE;
    fn_scenario = TRUE;
  }

  /* compute distance and capping, as appopriate */

  dist_error = FALSE;           /* was there an approximation error?
				   only applies to unsigned distance */
  if (unsigned_dist) {
    dist =
      unsigned_mcdist_nomem(&genome_list[gindex1],
			    &genome_list[gindex2],
			    NUM_GENES, NUM_CHROMOSOMES,
			    circular,
			    verbose,
			    &dist_error);
  } else  
  if (unichromosome) {
    dist =  invdist_noncircular_nomem_v(&genome_list[gindex1],
					&genome_list[gindex2],
					0, NUM_GENES,
					&graphstats);
  } else {
    if (fn_capping || fn_scenario || (fn_revstats == 1)) {
      dist =  mcdist_capgraph_nomem(&genome_list[gindex1],
				    &genome_list[gindex2],
				    NUM_GENES,NUM_CHROMOSOMES,
				    &graphstats);
#if DEBUGc
      dist2 = invdist_noncircular_nomem(&genome_list[gindex1],
				    &genome_list[gindex2],
				    0, NUM_GENES);
      if (dist != dist2) {
	fprintf(stderr,"WARNING: capping gave mcdist=%d, uniinvdist=%d\n",
		dist,dist2);
      }
#endif
    } else {
      dist =  mcdist_noncircular_nomem(&genome_list[gindex1],
				       &genome_list[gindex2],
				       NUM_GENES,NUM_CHROMOSOMES,
				       &graphstats);
    }
  }



  /***********************************************************************/
  /*                         Display report sections                     */
  /***********************************************************************/

  if (verbose) {
    if (unsigned_dist) {
      /* verbose output is displayed in the unsigned distance routines */
    } else if (unichromosome) {
      fprintf(outfile,"Number of Breakpoints:         \t%d\n", graphstats.br);
      fprintf(outfile,"Number of Cycles:              \t%d\n", graphstats.c4);
      fprintf(outfile,"Number of Hurdles:             \t%d\n", graphstats.h);
      fprintf(outfile,"Fortress:                      \t%d\n", graphstats.f);
    } else {
      fprintf(outfile,"Number of black edges in G:    \t%d\n", graphstats.bl);
      fprintf(outfile,"Number of cycles and paths:    \t%d\n", graphstats.cp);
      fprintf(outfile,"Number of gamma-gamma paths:   \t%d\n", graphstats.pgg);
#if 0
      fprintf(outfile,"Number of semi-knots:          \t%d\n", graphstats.s);
      fprintf(outfile,"Number of real knots in G-bar: \t%d\n", graphstats.rr);
#endif
      fprintf(outfile,"Number of semi-real-knots:     \t%d\n", graphstats.s);
      fprintf(outfile,"Number of real-knots:          \t%d\n", graphstats.r);
      fprintf(outfile,"Parameter gr:                  \t%d\n", graphstats.gr);
      fprintf(outfile,"Parameter fr:                  \t%d\n", graphstats.fr);
      if (graphstats.badbonds) {
	fprintf(outfile,"Number of bad bonds            \t%d\n",
		graphstats.badbonds);
      }
      fprintf(outfile,"Number of internal breakpoints:\t%d\n", graphstats.bp_int);
      fprintf(outfile,"Number of external breakpoints:\t%d\n", graphstats.bp_ext);
    }
  }



  if (dist_error)
    fprintf(outfile,"WARNING: These unsigned permutations are too complex; answers are approximated.\n");

  /* Display distance */
  if (fn_distance) {
    if (unichromosome)
      fprintf(outfile,"Reversal Distance:            \t%d\n",dist);
    else
      fprintf(outfile,"Multichromosomal Distance:    \t%d\n",dist);
  }




  /* Display capping */
  if (fn_capping) {
    fprintf(outfile,"\n======================================================================\n\n");

    if (unsigned_dist) {
      if (unichromosome) {
	if (circular)
	  fprintf(outfile,"Optimal signs and circular shifts:\n");
	else
	  fprintf(outfile,"Optimal signs:\n");
      } else {
	fprintf(outfile,"Optimal signs and capping:\n");
      }
    } else
    if (unichromosome) {
      if (circular)
	fprintf(outfile,"Optimal circular shifts:\n");
      else
	fprintf(outfile,"No capping is done for unichromosomal genomes:\n");
    } else {
      fprintf(outfile,"Optimal capping:\n");
    }
    if (fn_capping == 1) {
      print_genome_unichrom(&genome_list[gindex1], NUM_GENES);
      print_genome_unichrom(&genome_list[gindex2], NUM_GENES);
    } else {
      print_genome_multichrom_wcaps(&genome_list[gindex1],
				    NUM_GENES, NUM_CHROMOSOMES);
      print_genome_multichrom_wcaps(&genome_list[gindex2],
				    NUM_GENES, NUM_CHROMOSOMES);
    }
  }

  



  /* Display reversal scenario */
  if (fn_scenario) {
    fprintf(outfile,"\n======================================================================\n\n");


    if (unichromosome)
      fprintf(outfile,"An optimal sequence of reversals:\n");
    else
      fprintf(outfile,"An optimal sequence of rearrangements:\n");

    if (scenario_type == 1) {
      print_scenario(&genome_list[gindex1],
		     &genome_list[gindex2],
		     NUM_GENES, NUM_CHROMOSOMES);
    } else {
      print_scenario_2(&genome_list[gindex1],
		       &genome_list[gindex2],
		       NUM_GENES, NUM_CHROMOSOMES,
		       scenario_type);
    }
  }



  /* Display one-step reversal statistics */
  if (fn_revstats) {
    fprintf(outfile,"\n======================================================================\n\n");


    if (unichromosome)
      fprintf(outfile,"Test all single step reversals:\n\n");
    else
      fprintf(outfile,"Test all single step rearrangements:\n\n");

    if (fn_revstats == 1) 
      print_rev_stats_t1(&genome_list[gindex1],
			 &genome_list[gindex2],
			 NUM_GENES, NUM_CHROMOSOMES);
    else if (fn_revstats == 2 || fn_revstats == 4)
      print_rev_stats_t2(&genome_list[gindex1],
			 &genome_list[gindex2],
			 NUM_GENES, NUM_CHROMOSOMES,
			 fn_revstats == 4);
    else if (fn_revstats == 3)
      print_rev_stats_t3(&genome_list[gindex1],
			 &genome_list[gindex2],
			 NUM_GENES, NUM_CHROMOSOMES);
  }
 } /* END  if (!fn_matrix && !fn_texgraph) */



  /***********************************************************************/
  /*                   free allocated memory                             */
  /***********************************************************************/

  for (i=0 ; i<NUM_GENOMES ; i++) {
    free(genome_list[i].genes);
    free(genome_list[i].gnamePtr);
  }
  free(genome_list);

  return(0);
}

/****************************************************************************/
/*            Distance matrix for multiple genomes                          */
/****************************************************************************/


/* Stacia added this so that the ouput can be directly fed into other
 * programs such as the tds suite. It expects # of taxa on the first line
 * by itself, and then a full symmetric matrix with the taxa name as the
 * first thing on a line.
 */
#define GNAME_WIDTH 20
void print_full_distmatrix(int **distmatrix,int num_genomes, 
			   struct genome_struct *genome_list) {
  int i, j;     /* loop over matrix entries */
  int k;        /* loop to print genome name */
  int pad;

  int max_entry = 0;
  int entry_width;
  for (i=0; i<num_genomes; i++) {
    for (j=0; j<i; j++) {
      if (max_entry > distmatrix[j][i]) {
	max_entry = distmatrix[j][i];
      }
    }
  }
  entry_width = num_digits(max_entry);

  fprintf(outfile,"\n   %d\n",num_genomes);

  for (i=0 ; i<num_genomes ; i++) {

    /* print genome name, truncating or padding if necessary */
    for(k=0, pad=FALSE; k<GNAME_WIDTH; k++) {
      if (!pad) {
	if (genome_list[i].gnamePtr[k] != '\0')
	  fputc(genome_list[i].gnamePtr[k], outfile);
	else
	  pad = TRUE;
      }
      if (pad) fputc(' ', outfile);
    }
    fputc('\t', outfile);
    
    /* fprintf(outfile,"%s\t",genome_list[i].gnamePtr); */


    for (j=0; j<i ; j++)
      fprintf(outfile, "%*d ", entry_width, distmatrix[j][i]);
    fprintf(outfile, "%*d ", entry_width, 0);
    for (j=i+1 ; j<num_genomes ; j++) {
      fprintf(outfile,"%*d ", entry_width, distmatrix[i][j]);
    }
    fprintf(outfile,"\n");
  }
  fprintf(outfile,"\n");
  
  return;
}


/* print a matrix of the values of a particular field in the distance
 * structures */
void print_matrix_offset(char *matrix_name,
			 int field_offset,
			 graphstats_t **graphstats_matrix,
			 int num_genomes,
			 struct genome_struct *genome_list) {
  int i, j;     /* loop over matrix entries */
  int k;        /* loop to print genome name */
  int pad;

  int max_entry = 0;
  int e;
  int entry_width;
  for (i=0; i<num_genomes; i++) {
    for (j=0; j<i; j++) {
      e = *(int *)((char *)(&graphstats_matrix[i][j]) + field_offset);
      if (max_entry > e) {
	max_entry = e;
      }
    }
  }
  entry_width = num_digits(max_entry);


  fprintf(outfile,"\n");
  fprintf(outfile,"%s Matrix:\n", matrix_name);
  fprintf(outfile,"   %d\n",num_genomes);

  for (i=0 ; i<num_genomes ; i++) {

    /* print genome name, truncating or padding if necessary */
    for(k=0, pad=FALSE; k<GNAME_WIDTH; k++) {
      if (!pad) {
	if (genome_list[i].gnamePtr[k] != '\0')
	  fputc(genome_list[i].gnamePtr[k], outfile);
	else
	  pad = TRUE;
      }
      if (pad) fputc(' ', outfile);
    }
    fputc('\t', outfile);

    for (j=0 ; j<num_genomes ; j++) {
      fprintf(outfile, "%*d ", entry_width,
	      *(int *)((char *)(&graphstats_matrix[i][j]) + field_offset));
    }
    fprintf(outfile,"\n");
  }
  fprintf(outfile,"\n");
  
  return;
}




/*****************************************************************************/

/* circular option: the genomes have already been aligned so that
   linear signed inv dist = circular signed inv dist
   However, unsigned still needs to do more work */

void do_fn_distmatrix(int NUM_GENES, int NUM_CHROMOSOMES, int NUM_GENOMES,
		      int circular,
		      int unsigned_dist,
		      int verbose,
                      struct genome_struct *genome_list)
{
  int i;
  int **distmatrix;
  distmem_t distmem;
  int dist_error;

  if (verbose && !unsigned_dist) {
    do_fn_distmatrix_v(NUM_GENES, NUM_CHROMOSOMES, NUM_GENOMES,
		       circular, unsigned_dist,
		       genome_list);
    return;
  }


  fprintf(outfile,"\nDistance Matrix:\n");

    /* allocate space for matrix */
    distmatrix = (int **) e_malloc((NUM_GENOMES)*sizeof(int *), "distmatrix");

    for (i=0; i < NUM_GENOMES; i++) {
      distmatrix[i] = (int *) e_malloc((NUM_GENOMES)*sizeof(int), "distmatrix");
    }

    /* allocate memory for distance computations */
    mcdist_allocmem(NUM_GENES,NUM_CHROMOSOMES,&distmem);


    /* compute the appropriate matrix */
    if (unsigned_dist) {
      set_unsigneddist_matrix(distmatrix, genome_list,
			      NUM_GENES, NUM_CHROMOSOMES, NUM_GENOMES,
			      circular,
			      &distmem,
			      &dist_error);
      if (dist_error) {
	fprintf(outfile,"WARNING: These unsigned permutations are too complex;\n");
	fprintf(outfile,"some answers are approximated.\n");
      }

    } else if (NUM_CHROMOSOMES == 0) {
      /* unichromosomal data */ 
      setinvmatrix(distmatrix,genome_list,
		   NUM_GENES,NUM_GENOMES,
		   &distmem,

		   /* CIRCULAR already converted to equiv linear problem */
                   FALSE);
    } else {
      /* multichromosomal data */
      setmcdistmatrix(distmatrix,genome_list,
		      NUM_GENES,NUM_CHROMOSOMES,NUM_GENOMES,
		      &distmem);
    }

    /* display the matrix */
    print_full_distmatrix(distmatrix,NUM_GENOMES,genome_list);

    /* free the memory for distance computations */
    mcdist_freemem(&distmem);

    /* free memory for the distance matrix */
    for (i=NUM_GENOMES-1; i >=0 ; i--)
      free(distmatrix[i]);
    free(distmatrix);
}


/* verbose matrix option: display matrices of all the various parameters */
/* TODO: unsigned_dist not supported at this time */

#define print_matrix_field(display,var) \
    print_matrix_offset(display, \
		    offsetof(graphstats_t, var), \
		    statsmatrix, \
		    NUM_GENOMES, \
		    genome_list)



void do_fn_distmatrix_v(int NUM_GENES, int NUM_CHROMOSOMES, int NUM_GENOMES,
			int circular,
			int unsigned_dist,
			struct genome_struct *genome_list)
{
  int i,j;
  graphstats_t **statsmatrix;
  distmem_t distmem;


  /* allocate space for matrix */
  statsmatrix =
    (graphstats_t **) e_malloc((NUM_GENOMES)*sizeof(graphstats_t *),
			       "statsmatrix");

  for (i=0; i < NUM_GENOMES; i++) {
    statsmatrix[i] = (graphstats_t *) e_malloc((NUM_GENOMES)*sizeof(graphstats_t),
					     "statsmatrix");
  }

  /* allocate memory for distance computations */
  mcdist_allocmem(NUM_GENES,NUM_CHROMOSOMES,&distmem);


  /* compute all pairwise graphs.
   * Some parameters aren't symmetric, so do all pairs (i,j) and (j,i).
   * TODO: compute diagonal ones simply
   */

  for (i = 0 ;  i < NUM_GENOMES ;  i++) {
    for (j = 0 ;  j < NUM_GENOMES ;  j++) {
      mcdist_noncircular(&genome_list[i], &genome_list[j],
			 NUM_GENES, NUM_CHROMOSOMES,
			 &distmem,
			 &statsmatrix[i][j]);
    }
  }


  /* display the matrices */
  if (NUM_CHROMOSOMES == 0) {
    print_matrix_field("Distance", d);
    print_matrix_field("Number of Breakpoints", br);
    print_matrix_field("Number of Cycles", c4);
    print_matrix_field("Number of Hurdles", h);
    print_matrix_field("Number of Fortresses", f);
  } else {
    print_matrix_field("Distance", d);
    print_matrix_field("Number of Black Edges", bl);
    print_matrix_field("Number of Cycles and Paths", cp);
    print_matrix_field("Number of Gamma-Gamma Paths", pgg);
    print_matrix_field("Number of Semi-knots", s);
#if 0
    print_matrix_field("Parameter rr", rr);
#endif
    print_matrix_field("Parameter r", r);
    print_matrix_field("Parameter fr", fr);
    print_matrix_field("Parameter gr", gr);

#if 0
    /* can't do this w/o forming capping/concat, which we didn't do */
    print_matrix_field("Number of bad bonds", badbonds);
#endif

    print_matrix_field("Number of internal breakpoints", bp_int);
    print_matrix_field("Number of external breakpoints", bp_ext);
  }



  /* free the memory for distance computations */
  mcdist_freemem(&distmem);

  /* free memory for the distance matrix */
  for (i=NUM_GENOMES-1; i >=0 ; i--)
    free(statsmatrix[i]);
  free(statsmatrix);
}
