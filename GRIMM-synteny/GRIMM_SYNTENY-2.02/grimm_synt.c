/* grimm_synt.c
 *    GRIMM-Synteny: Create synteny blocks for GRIMM.
 *
 * Copyright (C) 2001-2006 The Regents of the University of California
 * by Glenn Tesler
 *
 * See file COPYRIGHT for details.
 *****************************************************************************
 * This file is part of GRIMM-Synteny.
 *
 * GRIMM-Synteny is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License, Version 2,
 * dated June 1991, as published by the Free Software Foundation.
 *
 * GRIMM-Synteny is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

/* Last modified on Wed Sep 6, 2006, by Glenn Tesler
 */

/*
 * GRIMM-Synteny 2.0
 * Glenn Tesler
 * April 14, 2004
 *
 * File format:
 * score chr1 start1 len1 sign1 chr2 start2 len2 sign2 ...
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
/*#include <values.h>*/
#include <ctype.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include "time_stamp.h"
#include "hash.h"
#include "ckalloc.h"
#include "parsestr.h"
#include "anctable.h"
#include "gsy.h"
#include "gscomp.h"
#include "gsfile.h"
#include "grimm_anch.h"

/* put quotes around version number */
#define qstr1(s) #s
#define qstr(s) qstr1(s)
#define QVERS qstr(VERS)


int main(int argc, char *argv[]);
void print_usage(char *progname);

void calc_num_anc_species(FILE *input,
			  int merge_chrnames,

			  int *nanc, int *nspecies,
			  char ***snames_ret,
			  HTABLE ***chrnames);


void anchors_alloc(ANCTABLE *anctable);
void read_anc(FILE *input, ANCTABLE *anctable, int same_species);
void anchors_destroy(ANCTABLE *anctable);
void chrnames_destroy(int nspecies, HTABLE **chrnames);
void snames_destroy(int nspecies, char **snames);


void compute_orders(FILE *output, ANCTABLE *anctable);
int compute_order_cmp(const void *a, const void *b);

void grimm_synt(FILE *report_file, ANCTABLE *anctable, GSparams *gs_params);
int load_anchors(FILE *output, char *inputfname, ANCTABLE *anctable,
		 int same_species);

int chrstr_cmp(const void *a0, const void *b0);
void sort_chrnames(int nspecies, HTABLE **chrnames);
void parse_chrnames(int nspecies, char *s, HTABLE **chrnames);


int scanlist(const char *s, int **list_ret, int *n_ret);
int resizelist(int **list, int n_in,
	       int n_out, int dflt);
void printlist(FILE *output, int *list, int n);
int sumlist(int *list, int n);
int minlist(int *list, int n);
int maxlist(int *list, int n);

void compute_perm_metric(FILE *report_file, ANCTABLE *anctable,
			 GSparams *gs_params);

void redo_list_anc_by_block(ANCTABLE *anctable, ANCTREE *anctree,
			    int ***blockanc_ret, int **node2blockpos_ret
			    );



/*
 * Needed for a global in GRIMM
 * TODO: Fix GRIMM!
 */
FILE *outfile;

static const char *signstrs0[3] = { "-", "0", "+" };
static const char *signstrs1[3] = { "-1", "0", "1" };

int main(int argc, char *argv[]) {
  /* file i/o */
  char *inputfname   = (char *) NULL;
  FILE *report_file  = (FILE *) NULL;
  struct stat statbuf;

  /* option parsing */
  int errflg = 0;
  int c;

  /* input data */ 
  ANCTABLE anctable_mem, *anctable=&anctable_mem;

  int nspecies;
  int i;

  /* options */
  GSparams gs_params_mem, *gs_params = &gs_params_mem;

  int gapthresh_s_n = 0;
  int min_block_size_s_n = 0;
  int min_block_size2_s_n = 0;
  int min_block_sup_s_n = 0;

  outfile = stdout;

  gs_params->gapthresh = MAXINT;
  gs_params->condense_strips = 0;
  gs_params->outputdir = (char *) NULL;
  gs_params->algorithm = 0;  /* 0 = GRIMM-Synteny, 1=GRIMM-Anchors */
  gs_params->min_block_size_s = (int *) NULL;
  gs_params->min_block_size2_s = (int *) NULL;
  gs_params->min_block_sup_s = (int *) NULL;
  gs_params->min_block_nanc = 1;
  gs_params->gapthresh_s = (int *) NULL;
  gs_params->min_corr = (double) 0.0;
  gs_params->sign_err_del = 0;
#if 0
  gs_params->tableau_d = MAXINT;
#endif
  gs_params->RS_complex_n = 0;
  gs_params->RS_complex_cells = (int *) NULL;
  gs_params->remove_overlap = 1;
  gs_params->format_style = 0;

  gs_params->perm_metric_use = 0;
  gs_params->perm_metric_spacing = 0;
  gs_params->perm_metric_length = 0;

  gs_params->sign_threshold = 2.0/3.0;
  gs_params->same_species = 0;


  while ((c = getopt(argc, argv, "Acf:d:g:G:m:M:n:r:R:S:OpP:QT:z")) != -1) {
    switch (c) {
    case 'f':
      inputfname = optarg;
      break;
#if 0
    /* Multiple output files, so use -d for output directory instead */
    case 'o':
      outputfname = optarg;
      break;
#endif
    case 'G':
      sscanf(optarg,"%d",&gs_params->gapthresh);
      break;
    case 'g':
      scanlist(optarg, &gs_params->gapthresh_s,&gapthresh_s_n);
      break;
    case 'm':
      /* minimum block size per-species */
      if (gs_params->perm_metric_use && *optarg == 'n') {
	scanlist(optarg+1,
		 &gs_params->min_block_size2_s,&min_block_size2_s_n);
      } else {
	scanlist(optarg,
		 &gs_params->min_block_size_s,&min_block_size_s_n);
      }
      break;
    case 'M':
      /* minimum block support per-species */
      scanlist(optarg,
	       &gs_params->min_block_sup_s,&min_block_sup_s_n);
      break;
    case 'n':
      /* minimum # anchors/block */
      sscanf(optarg,"%d",&gs_params->min_block_nanc);
      break;
    case 'c':
      gs_params->condense_strips = 1;
      break;
    case 'd':
      gs_params->outputdir = optarg;
      break;
    case 'r':
      sscanf(optarg, "%lf", &gs_params->min_corr);
      gs_params->sign_err_del = 1;
      break;
    case 'R':
      sscanf(optarg, "%lf", &gs_params->min_corr);
      gs_params->sign_err_del = 0;
      break;
#if 0
    case 'D':
      sscanf(optarg, "%d", &gs_params->tableau_d);
      break;
#endif
    case 'S':
      scanlist(optarg,
	       &gs_params->RS_complex_cells,
	       &gs_params->RS_complex_n);
      if (gs_params->RS_complex_n & 1) {
	errflg=1;
	fprintf(stderr,"-S should have an even number of values listed\n");
	break;
      }
      break;

    case 'O':
      gs_params->remove_overlap = 0;
      break;

    case 'Q':
      gs_params->format_style = 1;
      break;

    case 'p':
      gs_params->perm_metric_use = 1;
      gs_params->perm_metric_spacing = 2;
      gs_params->perm_metric_length = 2;
      break;

    case 'P':
      gs_params->perm_metric_use = 1;
      if (sscanf(optarg,"%d,%d",
	     &gs_params->perm_metric_spacing,
		 &gs_params->perm_metric_length) != 2
	  || gs_params->perm_metric_spacing <= 0
	  || gs_params->perm_metric_length <= 0
	  || gs_params->perm_metric_length > gs_params->perm_metric_spacing) {
	errflg++;
	fprintf(stderr,"-P: parameters invalid\n");
      }
      break;

      /* GRIMM-Anchors */
    case 'A':
      gs_params->algorithm = 1;
      break;

    case 'T': /* Sign threshold for GRIMM-Anchors */
      sscanf(optarg, "%lf", &gs_params->sign_threshold);
      if (gs_params->sign_threshold < 0.5)
	gs_params->sign_threshold = 0.5;
      if (gs_params->sign_threshold > 1)
	gs_params->sign_threshold = 1;
      break;

    case 'z':
      gs_params->same_species = 1;
      break;

      /* more options:
       * signed/unsigned
       * unichrom linear, unichrom circular, multichrom; more?
       * filtering options: min block sizes, block size ratios, etc.
       * input formats:
       *   GRIMM gene order
       *   Coordinatized (chr,start,end,sign)
       *   Coordinatized (chr,start,len,sign)
       *   different column arrangements, etc.?
       *
       * chromsome name and length file?
       * input:
       *    raw alignments, possibly w/conflicts               
       *    anchors w/o conflicts
       *    chromosome name/length files
       * output: files for each part of output
       *    report                              f.txt          d/report
       *    MGR microrearrangement input file   f_micro.mgr    d/mgr_micro
       *    MGR macrorearrangement input files  f_macro.mgr    d/mgr_macro
       *       (all combined into one file)
       *    large block coordinate file         f_blocks       d/blocks
       * GRIMM-anchors output
       *    unique anchors                                     d/anchors
       *    repeats                                            d/repeats
       *    special versions if self-self comparison?
       * repeats
       */

    }
  }

  if (gs_params->same_species && gs_params->algorithm != 1) {
    fprintf(stderr,"ERROR: -z option only available with -A (GRIMM-Anchors)\n");
    errflg = 1;

  }

  if (errflg) {
    print_usage(argv[0]);
    exit(-1);
  }

  /* open input */
  if (inputfname == (char *)NULL) {
    fprintf(stderr,"ERROR: input filename required\n");
    print_usage(argv[0]);
    exit(-1);
  }
  if (stat(inputfname, &statbuf)) {
    fprintf(stderr,"ERROR: Input file (%s): ",inputfname);
    perror(argv[0]);
    print_usage(argv[0]);
    exit(-1);
  }

  /* check output directory */
  if (gs_params->outputdir == (char *) NULL) {
    fprintf(stderr,"ERROR: output directory name required\n");
    print_usage(argv[0]);
    exit(-1);
  }

  if (stat(gs_params->outputdir, &statbuf)) {
    fprintf(stderr,"ERROR: Output directory (%s): ",gs_params->outputdir);
    perror(argv[0]);
    print_usage(argv[0]);
    exit(-1);
  }
  if (!(statbuf.st_mode & S_IFDIR)) {
    fprintf(stderr,"ERROR: Output directory name (%s) is not a directory.\n",
	    gs_params->outputdir);
    print_usage(argv[0]);
    exit(-1);
  }

  gs_params->signstrs =
    (char **) (gs_params->format_style ? signstrs1 : signstrs0);



  /* open report file */
  /* report_file = stdout; */

  report_file = fopen_dir(gs_params->outputdir,
			  gs_params->algorithm
			  ? "report_ga.txt" : "report.txt",
			  "w");
  setvbuf(report_file, (char *) NULL, _IONBF, (size_t) 0); /* unbuffer */


  time_stamp(report_file,"Start job");

  /* read input file */
  if (load_anchors(report_file,inputfname,
		   anctable,
		   gs_params->same_species)) {
    fprintf(stderr,"ERROR: Problem reading input file (%s): ",inputfname);
    perror(argv[0]);
    print_usage(argv[0]);
    exit(-1);
  }
  nspecies = anctable->nspecies;

  if (gs_params->algorithm) {
    /* GRIMM-Anchors */
    fprintf(report_file,"Parameters:\n");
    fprintf(report_file,"Algorithm:             GRIMM-Anchors\n");
    fprintf(report_file,"Input file:            %s\n", inputfname);

    if (gs_params->same_species) {
      fprintf(report_file,"  # species:           %d (species vs. self)\n", nspecies);
      fprintf(report_file,"  # alignments:        %d\n",
	      anctable->nanc/nspecies);
    } else {
      fprintf(report_file,"  # species:           %d\n", nspecies);
      fprintf(report_file,"  # alignments:        %d\n",
	      anctable->nanc);
    }

    fprintf(report_file,"Output directory:      %s\n", gs_params->outputdir);
    fprintf(report_file,"Tandem repeats         %s\n",
	    gs_params->condense_strips ? "Condense" : "Do not condense");
    fprintf(report_file,"Ambiguous sign thresh: %f\n",
	    gs_params->sign_threshold);
    fprintf(report_file,"\n");

    /* Compute order in each species
     * TODO: only need order in 1st species
     */
    time_stamp(report_file,"Sort anchors");
    compute_orders(report_file, anctable);
    grimm_anch(report_file, anctable, gs_params);

  } else {
    /* GRIMM-Synteny */

    /* further validate and optimize parameters, now that # species known */
    if (!resizelist(&gs_params->gapthresh_s, gapthresh_s_n, nspecies, MAXINT)) {
      fprintf(stderr,"ERROR: %d species but %d gap threshold(s)\n",
	      nspecies, gapthresh_s_n);
      print_usage(argv[0]);
      exit(-1);
    }

    if (!resizelist(&gs_params->min_block_size_s,
		   min_block_size_s_n,
		   nspecies,
		   0)) {
      fprintf(stderr,"ERROR: %d species but %d min block size limit(s)\n",
	      nspecies, min_block_size_s_n);
      print_usage(argv[0]);
      exit(-1);
    }

    if (gs_params->perm_metric_use
	&& !resizelist(&gs_params->min_block_size2_s,
		   min_block_size2_s_n,
		   nspecies,
		   0)) {
      fprintf(stderr,"ERROR: %d species but %d min block size (in nuc) limit(s)\n",
	      nspecies, min_block_size2_s_n);
      print_usage(argv[0]);
      exit(-1);
    }

    if (!resizelist(&gs_params->min_block_sup_s,
		   min_block_sup_s_n,
		   nspecies,
		   0)) {
      fprintf(stderr,"ERROR: %d species but %d min block support limit(s)\n",
	      nspecies, min_block_sup_s_n);
      print_usage(argv[0]);
      exit(-1);
    }


    /* optimize gap thresholds */
    gs_params->gapthresh_s_min = minlist(gs_params->gapthresh_s, nspecies);
    gs_params->gapthresh = min2(sumlist(gs_params->gapthresh_s, nspecies),
				gs_params->gapthresh);
    gs_params->min_block_size_s_max = maxlist(gs_params->min_block_size_s,
					      nspecies);
    if (gs_params->perm_metric_use) {
      gs_params->min_block_size2_s_max = maxlist(gs_params->min_block_size2_s,
						 nspecies);
    }

    fprintf(report_file,"Parameters:\n");
    fprintf(report_file,"Algorithm:             GRIMM-Synteny\n");
    fprintf(report_file,"Input file:            %s\n", inputfname);
    fprintf(report_file,"  # species:           %d\n", nspecies);
    fprintf(report_file,"  # anchors:           %d\n", anctable->nanc);

    for (i=0; i<nspecies; i++) {
      fprintf(report_file,
	      "  genome%-d:             %s\n",
	      i+1,
	      anctable->snames[i]);
    }

    fprintf(report_file,"Output directory:      %s\n", gs_params->outputdir);

    if (gs_params->perm_metric_use) {
      fprintf(report_file,"Metric:                permutation (spacing %d, length %d)\n",
	      gs_params->perm_metric_spacing,
	      gs_params->perm_metric_length);
    } else {
      fprintf(report_file,"Metric:                nucleotides\n");
    }

    fprintf(report_file,"Gap threshold:         %d\n", gs_params->gapthresh);

    fprintf(report_file,"  per species:         ");
    printlist(report_file, gs_params->gapthresh_s, nspecies);
    fprintf(report_file,"\n");
    
    fprintf(report_file,"Min # anchors/block:   %d\n", gs_params->min_block_nanc);

    if (gs_params->perm_metric_use) {
      fprintf(report_file,"Min block size/species (p): ");
      printlist(report_file, gs_params->min_block_size_s, nspecies);
      fprintf(report_file,"\n");

      fprintf(report_file,"Min block size/species (n): ");
      printlist(report_file, gs_params->min_block_size2_s, nspecies);
      fprintf(report_file,"\n");
    } else {
      fprintf(report_file,"Min block size/species: ");
      printlist(report_file, gs_params->min_block_size_s, nspecies);
      fprintf(report_file,"\n");
    }

    fprintf(report_file,"Min block support/species: ");
    printlist(report_file, gs_params->min_block_sup_s, nspecies);
    fprintf(report_file,"\n");

    fprintf(report_file,"Min correlation:       %f\n", gs_params->min_corr);
    fprintf(report_file,"  on low cor/bad sign: %s\n",
	    gs_params->sign_err_del
	    ? "Delete block"
	    : "Issue warning");
    fprintf(report_file,"Strips:                %s\n",
	    gs_params->condense_strips ? "Condense" : "Do not condense");
    fprintf(report_file,"Overlaps/containments: %s\n",
	    gs_params->remove_overlap ? "Repair" : "Do not repair");
    fprintf(report_file,"Formatting:            %s\n",
	    gs_params->format_style ? "Matlab friendly" : "Original");
    fprintf(report_file,"Permutation complexity\n");
#if 0
    fprintf(report_file,"  Reject RS diagonal >=%d\n", gs_params->tableau_d);
#endif
    fprintf(report_file,"  Reject RS cells      ");
    printlist(report_file,
	      gs_params->RS_complex_cells, gs_params->RS_complex_n);
    fprintf(report_file,"\n");


    grimm_synt(report_file, anctable, gs_params);
  }

  /* free species names */
  snames_destroy(nspecies, anctable->snames);

  /* free chromosome name hashes */
  chrnames_destroy(gs_params->same_species ? 1 : nspecies, 
		   anctable->chrnames);

  /* free anchors */
  anchors_destroy(anctable);

  /* free parameters */
  free((void *) gs_params->min_block_size_s);
  free((void *) gs_params->min_block_size2_s);
  free((void *) gs_params->gapthresh_s);
  free((void *) gs_params->RS_complex_cells);

  if (report_file != stdout && report_file != stderr)
    fclose(report_file);

  return(0);
}

void print_usage(char *progname)
{

  /* TODO: p|n prefix on -g, -M, -G */

  fprintf(stderr,"\nGRIMM-Synteny %s\n",QVERS);
  fprintf(stderr,"Copyright (C) 2001-2006 The Regents of the University of California\n");
  fprintf(stderr,"Contains code from GRAPPA 1.02 (C) 2000-2001 The University of New Mexico\n");
  fprintf(stderr,"                                         and The University of Texas at Austin\n");
  fprintf(stderr,"See file COPYRIGHT for full copyright and authorship details.\n");

  fprintf(stderr,
	  "\nUsage 1 (GRIMM-Synteny): %s -f inputfile -d outdir [other options]\n",
	  progname);
  fprintf(stderr, "    -f inputfile: file with anchor coordinates\n");
  fprintf(stderr, "    -d outdir: directory in which to write all output files\n");
  fprintf(stderr, "\n  Gap threshold:\n");
  fprintf(stderr, "    -g [#|#,#,...]: gap thresholds per-species\n");
  fprintf(stderr, "    -G #: Total gap threshold, using sum of intraspecies gaps\n");
  fprintf(stderr, "\n  Minimmum block size:\n");
  fprintf(stderr, "    -m [n][#|#,#,...]: minimum block size per-species in current metric,\n");
  fprintf(stderr, "       unless prefixed by letter n (nucleotide metric) to override -p.\n");
  fprintf(stderr, "       When using -p, can use two -m's, one with 'n' and one without.\n");
  fprintf(stderr, "    -M [#|#,#,...]: minimum support per block per-species\n");
  fprintf(stderr, "    -n #: minimum # anchors per block\n");
  fprintf(stderr, "\n  Other settings:\n");
  fprintf(stderr, "    -c: condense strips\n");
#if 0
  fprintf(stderr, " -D #: reject perms whose Robinson-Schensted diagonal length >= #\n");
#endif
  fprintf(stderr, "    -O: skip block overlap/containment analysis & repair\n");
  fprintf(stderr, "    -p: use permutation metric instead of nucleotides (same as -P 2,2)\n");
  fprintf(stderr, "    -P a,b: use permutation metric instead of nucleotides:\n");
  fprintf(stderr, "            anchor #i on a chromosome is at [a*i,a*i+b)\n");
  fprintf(stderr, "    -Q: alternate format of output files: signs 1,-1 instead of +,-;\n");
  fprintf(stderr, "        mgr_micro.txt in altered form to be readable by matlab as a matrix\n");
  fprintf(stderr, "\n  Incomplete/experimental:\n");
  fprintf(stderr, "    -S i1,j1,i2,j2,...: reject perms whose Robinson-Schensted shape has\n");
  fprintf(stderr, "       any cells (i1,j1), (i2,j2),... where these are 0-based matrix style\n");
  fprintf(stderr, "    -r #: minimum correlation coefficient for anchors; delete if smaller\n");
  fprintf(stderr, "    -R #: minimum correlation coefficient for anchors; warn if smaller\n");

  /* RESERVED:
  fprintf(stderr, "    ? i1,i2,...: restrict to input species i1,i2,... (0-based)\n");
  fprintf(stderr, "    -u: RESERVED (unsigned genes)\n");
  fprintf(stderr, "    -C: RESERVED (unichromosomal circular genomes)\n");
  fprintf(stderr, "    -L: RESERVED (unichromosomal linear genomes)\n");
  */

  fprintf(stderr, "\n");
  fprintf(stderr, "Usage 2 (GRIMM-Anchors): %s -A -f inputfile -d outdir [other options]\n",
	  progname);
  fprintf(stderr, "    -A: perform the GRIMM-Anchors algorithm.  Required for this usage.\n");
  fprintf(stderr, "    -f inputfile: file with anchor coordinates\n");
  fprintf(stderr, "    -d outdir: directory in which to write all output files\n");
  fprintf(stderr, "    -c: condense strips of tandem repeats into singletons (EXPERIMENTAL)\n");
  fprintf(stderr, "    -T #: value between .5 and 1.  When repeat signs are ambiguous, choose a\n");
  fprintf(stderr, "          sign if a fraction > # of the inputs are consistent with that sign.\n");
  fprintf(stderr, "    -z: Input file is species-on-self alignment coordinates.\n");

}


/*
 * parse comma-separated list of integers
 *
 * return value:
 *    1: good
 *    0: error
 * return parameters:
 *    *list_ret = {x1,x2,...,xn}
 *    *n_ret = n = size of list
 */
int scanlist(const char *s, int **list_ret, int *n_ret)
{
  int n=1;
  int i;
  int *list;
  const char *s1;

  /* count # fields */
  s1 = s;
  while (*s1 != '\0') {
    /* 1st char must be digit */
    if (!isdigit((int) *s)) {
      *list_ret = (int *) NULL;
      *n_ret = 0;
      return 0;
    }

    /* scan till comma */
    while (*s1 != '\0' && *s1 != ',') s1++;

    if (*s1 == ',') {
      n++; s1++;
    }
  }

  /* parse them */
  *list_ret = list = (int *) ckalloc(n, sizeof(int));
  *n_ret = n;

  s1 = s;
  i=0;
  while (*s1 != '\0') {
    if (!sscanf(s1, "%d", &list[i])) {
      free((void *) list);
      *list_ret = (int *) NULL;
      *n_ret = 0;
      return 0;
    }

    /* scan till comma */
    while (*s1 != '\0' && *s1 != ',') s1++;

    if (*s1 == ',') {
      i++; s1++;
    }
  }

  /* success */
  return 1;
}

/*
 * input: *list is a list of integers of size n_in
 *
 * output: *list either stays as-is or is replaced by one of size n_out
 *
 * how:
 *   if n_in=n_out, return
 *   if n_in=0, make it (dflt,dflt,...)
 *   if n_in=1, make it (list[0],list[0],...)
 *   if n_in>1 and n_in != n_out, return error
 *
 * return value:
 *   1: good
 *   0: error
 */

int resizelist(int **list, int n_in,
	       int n_out, int dflt)
{
  int i;

  if (n_in == n_out)
    return 1;  /* success, nothing to do */

  if (n_in > 1)
    return 0;  /* error, list was given but size is incorrect */

  if (n_in == 1)
    dflt = (*list)[0];

  free((void *) *list);
  *list = (int *) ckalloc(n_out, sizeof(int));
  for (i=0; i<n_out; i++)
    (*list)[i] = dflt;

  return 1;
}

void printlist(FILE *output, int *list, int n)
{
  int i;
  for (i=0; i<n; i++) {
    if (i>0)
      fprintf(output, " ");
    fprintf(output, "%d", list[i]);
  }
}

/* sum of a list */
int sumlist(int *list, int n)
{
  int acc=0;
  int acc2;
  int i;
  for (i=0; i<n; i++) {

    /* deal with sums with overflows */
    acc2 = acc+list[i];

    /* detect overflow
     * TODO: negative overflow
     */
    if (acc>=0 && list[i]>=0 && acc2<0) {
      acc = MAXINT;
    } else {
      acc = acc2;
    }
  }

  return acc;
}

/* minimum or maximum of a list */
int minlist(int *list, int n)
{
  int acc = MAXINT;
  int i;
  for (i=0; i<n; i++)
    if (list[i] < acc) acc=list[i];

  return acc;
}
int maxlist(int *list, int n)
{
  int acc = ~MAXINT;
  int i;
  for (i=0; i<n; i++)
    if (list[i] > acc) acc=list[i];

  return acc;
}



/* read anchors from inputfname into *anctable */
int load_anchors(FILE *report_file, char *inputfname, ANCTABLE *anctable,
		 int same_species)
{
  FILE *input;
  HTABLE **chrnames;
  char **snames;

  input = fopen(inputfname, "r");
  if (input == (FILE *) NULL) {
#if 0
    fprintf(stderr,"ERROR: Could not open input file (%s): ",inputfname);
#endif
    return -1;
  }

  time_stamp(report_file,"Scan anchors");

  /* calculate number of anchors and number of species */
  calc_num_anc_species(input,
		       same_species,
		       &anctable->nanc, &anctable->nspecies,
		       &snames,
		       &chrnames);

  if (same_species) {
#if 0
    if (anctable->nspecies != 2) {
      fprintf(stderr,"ERROR: -z option requires file formatted with exactly 2 species\n");
      fclose(input);
      return 1;
    }
#endif
    /* each anchor will be duplicated/shifted so that each species ends
     * up in 1st position once.
     */
    anctable->nanc *= anctable->nspecies;
  }

#if 0
  fprintf(report_file,"%d anchors in %d species\n",
	  anctable->nanc, anctable->nspecies);
#endif
  anctable->fields_per_species = FIELDR_n_per_species;
  anctable->fields_per_anchor =
    anctable->fields_per_species * anctable->nspecies + FIELD_base;
  anctable->chrnames = chrnames;
  anctable->snames = snames;

  /* allocate space for anchors */
  time_stamp(report_file,"Allocate space");
  anchors_alloc(anctable);

  /* read them in */
  time_stamp(report_file,"Read anchors");
  read_anc(input, anctable, same_species);

  fclose(input);

  return 0;
}

/* grimm synteny */
/* TODO: move to gscomp.c? */
void grimm_synt(FILE *report_file, ANCTABLE *anctable, GSparams *gs_params)
{
  ANCTABLE *blocktable;
  ANCTREE *anctree;
  int **blockanc, *node2blockpos;
  int max_nanc; /* max # anchors/block */
  int redo_olap, redo_signs, recalc, redo_RS;
  int did_RS = 0;

#if 0
  FILE *other_out = report_file;  /* fix later */
#endif


  if (gs_params->perm_metric_use) {
    time_stamp(report_file,"Computing permutation metric");
    compute_perm_metric(report_file, anctable, gs_params);
  }


  /* form components */
  time_stamp(report_file,"Form components");
  gstree(anctable,gs_params,
	 &anctree);


  /* TODO:
   * consolidate all the redo_* flags into one flag plus another giving reason
   */
  redo_signs = redo_olap = 1;
  recalc = 0;
  while (redo_signs || redo_olap || redo_RS || recalc) {
    time_stamp(report_file,"Determine component coordinates");
    gs_comp_coords(anctable, gs_params, anctree,
		   &blocktable, &max_nanc);

    time_stamp(report_file,"Filter out small components");
    gs_filter_small(anctable, gs_params, anctree, blocktable);

    time_stamp(report_file,"Sort anchors by component then species");
    list_anc_by_block(anctable, anctree,
		      &blockanc, &node2blockpos);

    time_stamp(report_file,"Determine signs of blocks");
    redo_signs =
      gs_signs(report_file,
	       gs_params,
	       anctable, anctree, blocktable, max_nanc,
	       blockanc, node2blockpos);
    if (redo_signs) {
      time_stamp(report_file,"Sort anchors by component then species");
      redo_list_anc_by_block(anctable, anctree,
			     &blockanc, &node2blockpos);
      redo_signs = 0;  /* no need to redo rest of analysis */
    }

    redo_RS = 0;
    if (!did_RS) {
#if 0
      if (gs_params->tableau_d < MAXINT) {
#endif
      if (gs_params->RS_complex_n > 0) {
	time_stamp(report_file,"Filter complex permutations");
	redo_RS =
	  gs_filter_perms(report_file,
			  gs_params,
			  anctable, anctree, blocktable, max_nanc,
			  blockanc, node2blockpos);
      }
      did_RS = 1;
    }


    if (gs_params->remove_overlap) {
      time_stamp(report_file,"Block overlap analysis");
      redo_olap = gs_overlap(report_file,
			     gs_params,
			     anctable, anctree,
			     blocktable, blockanc, node2blockpos);
    } else {
      redo_olap = 0;
    }

    if (!redo_signs && !redo_olap && !redo_RS && !recalc && gs_params->condense_strips) {
      time_stamp(report_file,"Condensing strips of consecutive blocks");
      gs_combine_consec_blocks(anctable, anctree, blocktable,
			       blockanc, node2blockpos);
      /* need to do another pass to recompute block coordinates */
      recalc = 1;
    } else {
      recalc = 0;
    }

    if (!redo_signs && !redo_olap && !redo_RS && !recalc) break;

    if (redo_signs || redo_RS) {
      time_stamp(report_file,"Blocks were discarded, must recompute several steps");
      recalc = 0;
    } else if (redo_olap) {
      time_stamp(report_file,"Blocks were split, must recompute several steps");
      recalc = 0;
    } else {
      time_stamp(report_file,"Blocks were condensed, must recompute several steps");
    }

    free((void *) node2blockpos);
    free((void *) blocktable);
    array2d_destroy(anctable->nspecies, (void *) blockanc);
  }

  time_stamp(report_file,"Microrearrangement analysis");
  fprintf(report_file,"Block lengths and support are genome-by-genome.\n");
  fprintf(report_file,"Other parameters are a below-the-diagonal triangular genome vs. genome matrix\nprinted in 1 row.\n");
  gs_micro(report_file,
	   anctable, anctree, blocktable, max_nanc, blockanc, node2blockpos,
	   gs_params->RS_complex_n);

  time_stamp(report_file,"Macrorearrangement analysis");
  gs_macro(report_file,
	   gs_params,
	   anctable, anctree, blocktable, blockanc, node2blockpos);

  time_stamp(report_file,"MGR microrearrangement files");
  gs_micro_mgr(gs_params,
	       anctable, anctree, blocktable, max_nanc,
	       blockanc, node2blockpos);

  time_stamp(report_file,"Breakpoint graph complexity analysis TODO");
  time_stamp(report_file,"Nadeau-Taylor model tests TODO");

#if 0
  /* display all anchors with components */
  time_stamp(report_file,"Write out anchors");
  display_components_byanchor(other_out, anctable, anctree);
#endif

  /* clear space */
  time_stamp(report_file,"Free space");
  gstree_destroy(anctree);
  free((void *) node2blockpos);
  free((void *) blocktable);
  array2d_destroy(anctable->nspecies, (void *) blockanc);

  time_stamp(report_file,"Done");
}


void redo_list_anc_by_block(ANCTABLE *anctable, ANCTREE *anctree,
			    int ***blockanc_ret, int **node2blockpos_ret
			    )
{
  free((void *) (*node2blockpos_ret));
  array2d_destroy(anctable->nspecies, (void *)(*blockanc_ret));
  list_anc_by_block(anctable,anctree,blockanc_ret,node2blockpos_ret);
}



#if 0
void calc_num_anc_species_old(FILE *input, int *nanc, int *nspecies)
{
  char buf[BUFSIZ], *curline;
  int c;

  *nanc = 0;
  *nspecies = 0;

  curline = fgets(buf,BUFSIZ,input);
  if (curline == (char *) NULL) {
    fprintf(stderr,"ERROR: null input file\n");
    exit(-1);
  }
  *nspecies = (count_fields(curline) - 1) / 4;
  ++*nanc;

  while ((c = getc(input)) != EOF) {
    if (c == '\n') ++*nanc;
  }
}
#endif

/*
 * scan input
 *   merge_chrnames = 0: each genome has own list of chromosome names
 *                    1: merged list of chromosome names across all genomes
 * return:
 *   nanc = # anchors
 *   nspecies = # species
 *   (*chrnames_ret)[j] = sorted list of chromosome names in species j,
 *                        stored in hash.
 *                        If merged_chrnames==1, all ptrs to hashes are
 *                        the same.
 *   (*snames_ret)[j] = name of species j
 */
void calc_num_anc_species(FILE *input,
			  int merge_chrnames,

			  int *nanc_ret,
			  int *nspecies_ret,
			  char ***snames_ret,
			  HTABLE ***chrnames_ret
			  )
{
  char buf[BUFSIZ], *curline;

  int nanc = 0;
  int nspecies = 0;

  HTABLE **chrnames;    /* chromosome names per species */

  int nnames = 0;
  char **snames;        /* species names */
  char sname_initial;
  int sname_num;

  int i;
  int nf,nf0;
  int namewidth;
  int incomments;

  /* file may start with comment lines, including genome names */
  /* scan initial comment lines, if any */
  incomments = 1;
  while (incomments) {
    curline = fgets(buf,BUFSIZ,input);
    if (curline == (char *) NULL) {
      fprintf(stderr,"ERROR: null input file\n");
      exit(-1);
    }
    curline = skip_blanks(curline);
    if (*curline == '\0') {
      continue;
    }
    if (*curline == '#') {
      /* see if a species name is defined */
      if (sscanf(curline,"# genome%d: %c", &sname_num, &sname_initial) == 2) {
	nnames++;
      }
    } else {
      incomments = 0;
    }
  }

  /* now count # species on first line */
  nf0 = count_fields(curline);
  nspecies = (nf0 - 1) / 4;
  if (nf0 != nspecies*4+1) {
    fprintf(stderr,
	    "ERROR: input file starts with %d fields per line, which is invalid.\n",
	    nf0);
    exit(-1);
  }

  /* check for agreement in # declared names vs. # speces per line */
  if (nspecies <= 0 || (nnames > 0 && nnames != nspecies)) {
    fprintf(stderr,
	    "ERROR: input file declares %d genome names but data has %d species\n",
	    nnames, nspecies);
    exit(-1);
  }

  /* allocate space for array of species names */
  snames = (char **) ckalloc(nspecies, sizeof(char *));


  /* get the genome names */
  if (nnames == 0) {
    namewidth =
      strlen("genome")
      + ceil(log(nspecies)/log(10)) /* length of number */
      + 2;                          /* for rounding error + terminal null */
    for (i=0; i<nspecies; i++) {
      snames[i] = (char *) ckalloc(1, namewidth * sizeof(char));
      snprintf(snames[i],namewidth,
	       "genome%d", i+1);
    }
  } else {
    /* rescan initial comments */
    rewind(input);
    i = 0;

    incomments = 1;
    while (incomments) {
      curline = fgets(buf,BUFSIZ,input);
      if (curline == (char *) NULL) {
	fprintf(stderr,"ERROR: null input file\n");
	exit(-1);
      }
      curline = skip_blanks(curline);
      if (*curline == '\0') {
	continue;
      }
      if (*curline == '#') {
	/* see if a species name is defined */
	if (sscanf(curline,"# genome%d: %c", &sname_num, &sname_initial) == 2) {
	  /* Skip "# genome%d:": */
	  while (*curline != ':') { curline++; }
	  curline++;
	  curline = skip_blanks(curline);

	  /* get genome name */
	  /* space needed for it, incl. terminating null */
	  namewidth = strlen(curline) + 1;

	  /* strip off trailing whitespace, incl. cr/lf */
	  while (namewidth >= 2
		 && isspace((int) (curline[namewidth-2]))) {
	    curline[namewidth-2] = '\0';
	    namewidth--;
	  }


	  /* allocate space and store genome name */
	  snames[i] = (char *) ckalloc(1, namewidth * sizeof(char));
	  strncpy(snames[i], curline, namewidth);

	  i++;
	}
      } else {
	incomments = 0;
      }
    }
  }

  /* Got species names.  Now do first anchor line + rest of file */

  /* allocate space for hash tables for chromosome names */
  chrnames = (HTABLE **) ckalloc(nspecies, sizeof(HTABLE *));
  if (merge_chrnames) {
    chrnames[0] = hash_create(MAXCHR);
    for (i=1; i<nspecies; i++) {
      chrnames[i] = chrnames[0];
    }
  } else {
    for (i=0; i<nspecies; i++) {
      chrnames[i] = hash_create(MAXCHR);
    }
  }

  /* parse chromosome names in 1st line */
  parse_chrnames(nspecies, curline, chrnames);
  nanc++;

  /* parse chromosome names in remaining lines */
  while ((curline = fgets(buf,BUFSIZ,input)) != (char *) NULL) {
    curline = skip_blanks(curline);
    if (*curline == '\0' || *curline == '#') {
      continue;
    }

    nf = count_fields(curline);
    if (nf != nf0) {
      fprintf(stderr,
	      "ERROR: input file starts with %d fields per line but this line has %d fields:\n%s\n",
	      nf0,nf,curline);
      exit(-1);
    }
    parse_chrnames(nspecies, curline, chrnames);
    nanc++;
  }

  /* re-sort chromosome names */
  sort_chrnames(merge_chrnames ? 1 : nspecies,
		chrnames);

  *nanc_ret = nanc;
  *nspecies_ret = nspecies;
  *chrnames_ret = chrnames;
  *snames_ret = snames;
}

/*
 * parse chromosome names from input line and add to chromosome table
 */
void parse_chrnames(int nspecies, char *s, HTABLE **chrnames)
{
  char chr_str[BUFSIZ];

  int j;

  s = skip_blanks(s);
  s = skip_nonblanks(s); /* ignore score */
  s = skip_blanks(s);

  for (j=0; j<nspecies; j++) {
    sscanf(s,"%s",chr_str);
    /* add chromosome name to table */
    hash_str2num(chrnames[j],chr_str);

    s = skip_nonblanks(s); /* skip over chromosome name */
    s = skip_blanks(s);

    s = skip_nonblanks(s); /* skip over start */
    s = skip_blanks(s);

    s = skip_nonblanks(s); /* skip over length */
    s = skip_blanks(s);

    s = skip_nonblanks(s); /* skip over sign */
    s = skip_blanks(s);
  }
}

/* re-assign chromosome name->number maps */
void sort_chrnames(int nspecies, HTABLE **chrnames)
{
  int j;

  /* re-assign numbers to chromosomes for each species */
  for (j=0; j<nspecies; j++)
    hash_sort(chrnames[j], chrstr_cmp);
}

/* compare chromosome names in special order:
 * 1 < 2 < 3 < ... < A < A1 < ... < B < B1 < ...
 */

int chrstr_cmp(const void *a0, const void *b0)
{
  const char *a = *(const char **) a0;
  const char *b = *(const char **) b0;
  int i=0;
  int n_a, n_b;

  /* find first different character */
  while (a[i] != '\0'
	 && b[i] != '\0'
	 && tolower((int) a[i]) == tolower((int) b[i])) {
    i++;
  }



  if (a[i] == '\0') {
    return
      (b[i] == '\0')
      ? 0    /* equal */
      : -1;  /* a smaller */
  }
  if (b[i] == '\0') {
    return 1;  /* b smaller */
  }

  /* if numeric, compare as numbers;
   * if one is numeric and other isn't, number is smaller;
   * else compare as strings
   */
  if (isdigit((int) a[i])) {
    if (isdigit((int) b[i])) {
      sscanf(&a[i], "%d", &n_a);
      sscanf(&b[i], "%d", &n_b);
      return n_a - n_b;
    }
    return -1;  /* a[i] numeric, b[i] not numeric, a is smaller */
  }
  if (isdigit((int) b[i]))
    return 1;   /* b[i] numeric, a[i] not numeric, b is smaller */

  /* compare as strings */
  return tolower((int) a[i]) - tolower((int) b[i]);
}


/* allocate space for anchors */
void anchors_alloc(ANCTABLE *anctable)
{
  int nanc = anctable->nanc;
  int nspecies = anctable->nspecies;
  int fields_per_anchor = anctable->fields_per_species * nspecies + FIELD_base;


  /* allocate space for anchors */
  anctable->anclist =
    int_array2d_create(nanc, fields_per_anchor);
  anctable->anclist2 = (int **) NULL;
  anctable->anclist3 = (int **) NULL;
  anctable->orders = (int **) NULL;

#if 0
  /* allocate space for hash tables */
  anctable->chrnames =
    chrnames = (HTABLE **) ckalloc(nspecies, sizeof(HTABLE *));
  for (i=0; i<nspecies; i++) {
    chrnames[i] = hash_create(MAXCHR);
  }
#endif
}

/* clear space */
void anchors_destroy(ANCTABLE *anctable)
{
  /* clear space for anchors */
  array2d_destroy(anctable->nanc, (void **) anctable->anclist);
  array2d_destroy(anctable->nanc, (void **) anctable->anclist2);
  array2d_destroy(anctable->nanc, (void **) anctable->anclist3);
  array2d_destroy(anctable->nspecies, (void **) anctable->orders);
}

void chrnames_destroy(int nspecies, HTABLE **chrnames)
{
  int i;

  /* clear hash tables */
  for (i=0; i<nspecies; i++) {
    hash_destroy(chrnames[i]);
  }
  free((void *) chrnames);
}

void snames_destroy(int nspecies, char **snames)
{
  array2d_destroy(nspecies, (void *) snames);
}

/* read in the anchors */
void read_anc(FILE *input, ANCTABLE *anctable, int same_species)
{
  int nanc = anctable->nanc;
  int nspecies = anctable->nspecies;
  int fields_per_species = anctable->fields_per_species;
  int **anclist = anctable->anclist;
  HTABLE **chrnames = anctable->chrnames;

  char buf[BUFSIZ], *s;
  char chr_str[BUFSIZ], sign_str[BUFSIZ];

  int chr_num, start, len, sign;

  int ID;

  int base_field;
  int i,j;

  int max_chr;
  int gotline;

  int dup;
  int *anc_dup,*anc_i;

  rewind(input);
  for (i=0; i<nanc; i++) {

    /* skip comment lines */
    gotline = 0;
    while (!gotline) {
      s = fgets(buf,BUFSIZ,input);
      s = skip_blanks(s);
      if (*s != '#' && *s != '\0') { gotline = 1; }
    }

    /* since we have previously checked that all lines are
     * either comments or else have the same number of fields,
     * we do not need to check for terminating '#' on
     * lines that have coordinates + comments; we will already
     * stop parsing just short of it.
     */

    s = skip_blanks(s);
    sscanf(s,"%d", &ID);
    s = skip_nonblanks(s);
    s = skip_blanks(s);
    anclist[i][FIELD_id] = ID;



    for (j=0; j<nspecies; j++) {
      base_field = j*fields_per_species + FIELD_base;

      /* fields: chromosome, start, length, sign */

      sscanf(s,"%s",chr_str);
      chr_num = hash_str2num(chrnames[j],chr_str);

      s = skip_nonblanks(s);
      s = skip_blanks(s);
      sscanf(s,"%d",&start);

      s = skip_nonblanks(s);
      s = skip_blanks(s);
      sscanf(s,"%d", &len);

      s = skip_nonblanks(s);
      s = skip_blanks(s);
      sscanf(s,"%s", sign_str);
      sign = (sign_str[0] == '-') ? -1 : 1;

      s = skip_nonblanks(s);
      s = skip_blanks(s);

      anclist[i][base_field + FIELDR_chr] = chr_num;
      anclist[i][base_field + FIELDR_start] = start;
      anclist[i][base_field + FIELDR_len] = len;
      anclist[i][base_field + FIELDR_sign] = sign;
    }

    /* normalize so sign in first species is + */
    if (anclist[i][FIELD_base + FIELDR_sign + 0*fields_per_species] < 0) {
      for (j=0; j<nspecies; j++) {
	base_field = j*fields_per_species + FIELD_base + FIELDR_sign;
	anclist[i][base_field] = -anclist[i][base_field];
      }
    }

    /* For species-on-self, make new lines with same coords but shifted
     * input line: ID c0 c1 c2    (c0 = chrom0 start0 len0 sign0, etc.)
     * becomes multiple lines:
     *     ID c0 c1 c2
     *     ID c2 c0 c1
     *     ID c1 c2 c0
     */
    if (same_species) {
      anc_i = anclist[i];
      for (dup=1; dup<nspecies; dup++) {
	anc_dup = anclist[i+dup];
	anc_dup[FIELD_id] = anc_i[FIELD_id];

	/* duplicate other fields but with species shifted */
	memcpy(&anc_dup[FIELD_base + 0*fields_per_species],
	       &anc_i  [FIELD_base + (nspecies-dup)*fields_per_species],
	       dup * fields_per_species * sizeof(int));
	memcpy(&anc_dup[FIELD_base + dup*fields_per_species],
	       &anc_i  [FIELD_base + 0*fields_per_species],
	       (nspecies-dup) * fields_per_species * sizeof(int));
      }

      /* normalize so sign in first species is + */
      if (anc_dup[FIELD_base + FIELDR_sign + 0*fields_per_species] < 0) {
	for (j=0; j<nspecies; j++) {
	  base_field = j*fields_per_species + FIELD_base + FIELDR_sign;
	  anc_dup[base_field] = -anc_dup[base_field];
	}
      }

      i += nspecies-1;     /* the for loop has the final i++ */
    }

  }

  /* determine maximum number of chromosomes over all species */
  max_chr = 0;
  for (j=0; j<nspecies; j++) {
    if (chrnames[j]->len > max_chr)
      max_chr = chrnames[j]->len;
  }
  anctable->max_chr = max_chr;
}



static int compute_order_specno;
static int compute_order_base;
static ANCTABLE *compute_order_anctable;

void compute_orders(FILE *output, ANCTABLE *anctable)
{
  int nanc = anctable->nanc;
  int nspecies = anctable->nspecies;
  int fields_per_species = anctable->fields_per_species;
  int **anclist = anctable->anclist;
  int i,k;
  int *order_i;

  char msg[BUFSIZ];

  anctable->orders = (int **) ckalloc(nspecies, sizeof(int *));

  for (i=0; i<nspecies; i++) {
    sprintf(msg, "Sorting species %d", i+1);
    time_stamp(output,msg);

    /* allocate table */
    order_i = anctable->orders[i] = (int *) ckalloc(nanc, sizeof(int));

    for (k=0; k<nanc; k++)
      order_i[k] = k;

    /* sort anchors into order for species i */
    compute_order_specno = i;
    compute_order_base = i*fields_per_species + FIELD_base;
    compute_order_anctable = anctable;

    qsort((void *) order_i, nanc,  sizeof(int),
	  compute_order_cmp);


    /* store orders in anchor */
    for (k=0; k<nanc; k++) {
      anclist[order_i[k]][compute_order_base + FIELDR_order] = k;
    }
  }
}

int compute_order_cmp(const void *a, const void *b)
{
  int *anc1 = compute_order_anctable->anclist[*(int *)a];
  int *anc2 = compute_order_anctable->anclist[*(int *)b];

  int d =
    anc1[compute_order_base + FIELDR_chr]
    - anc2[compute_order_base + FIELDR_chr];

  if (d == 0)
    d =
      anc1[compute_order_base + FIELDR_start]
      - anc2[compute_order_base + FIELDR_start];

  if (d == 0)
    d =
      anc1[compute_order_base + FIELDR_len]
      - anc2[compute_order_base + FIELDR_len];

  if (d == 0)
    d = *(int *)a - *(int *)b;

  return d;
}

void compute_perm_metric(FILE *report_file, ANCTABLE *anctable,
			 GSparams *gs_params)
{
  int nanc = anctable->nanc;
  int nspecies = anctable->nspecies;
  int fields_per_species = anctable->fields_per_species;
  int fields_per_anchor = anctable->fields_per_anchor;
  int **anclist_n = anctable->anclist;
  int **anclist_p = int_array2d_create(nanc, fields_per_anchor);
  int i,k;
  int *order_i;

  int perm_metric_spacing = gs_params->perm_metric_spacing;
  int perm_metric_length = gs_params->perm_metric_length;

  int lastchr, curchr, start;
  int col_chr, col_start, col_len;
  int r;

  char msg[BUFSIZ];


  /* 1. duplicate the list of anchors
   * 2. change metric in the duplicate
   */
  int_array2d_cpy(anclist_p, (const int **) anclist_n,
		  nanc, fields_per_anchor);

  order_i = (int *) ckalloc(nanc, sizeof(int));

  for (i=0; i<nspecies; i++) {
    sprintf(msg, "Sorting species %d", i+1);
    time_stamp(report_file,msg);

    for (k=0; k<nanc; k++)
      order_i[k] = k;

    /* sort anchors into order for species i */
    compute_order_specno = i;
    compute_order_base = i*fields_per_species + FIELD_base;
    compute_order_anctable = anctable;

    qsort((void *) order_i, nanc,  sizeof(int),
	  compute_order_cmp);

    /* recompute metric in nucleotides */

    /* store orders in anchor */
    lastchr = -1;
    start = 0;

    col_chr = compute_order_base + FIELDR_chr;
    col_start = compute_order_base + FIELDR_start;
    col_len = compute_order_base + FIELDR_len;

    for (k=0; k<nanc; k++) {
      r = order_i[k];
      curchr = anclist_p[r][col_chr];
      if (curchr != lastchr) {
	/* restart coordinates at 0 in each chromosome */
	lastchr = curchr;
	start = 0;
      } else {
	start += perm_metric_spacing;
      }
      anclist_p[r][col_start] = start;
      anclist_p[r][col_len] = perm_metric_length;
    }
  }

  free((void *) order_i);

  anctable->anclist = anclist_p;
  anctable->anclist2 = anclist_n;
}
