/* unsignedhc.c
 *    Fast approximation algorithm to find signs in partially signed genomes.
 *
 * Copyright (C) 2001-2006 The Regents of the University of California
 * by Glenn Tesler
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

/* Last modified on Tue Aug 1, 2006, by Glenn Tesler
 */

/* Glenn Tesler
 * Sep 5, 2004
 * Fast approximation to unsigned problem with multiple genomes,
 * partial sign info, and multiple segments:
 *
 *    >genome1
 *    1 2 3 4 5 6 7 $
 *    8 9 10 11 12 13 $
 *    >genome2
 *    10 [ 4 5 6 ] 11 12 $
 *    [1] [2 3] [7] -9 8 $
 *    13
 *    >genome3
 *    ...
 *
 * The regions identified by [...] have unknown overall sign,
 * could be [4 5 6] or [-6 -5 -4].
 * Use hill-climbing algorithm to find minimal score:
 *     D = sum_{i<j} d(genome[i],genome[j])
 * using the appropriate distance function d
 * (multichrom, unichrom circ, unichrom lin)
 * and also find sign patterns in all solutions tied for best score.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcstructs.h"
#include "uniinvdist.h"
#include "unsignedhc.h"
#include "scenario.h"
#include "mcrdist.h"
#include "e_malloc.h"
#include "write_data.h"
#include "matrixmisc.h"

struct uhc_mem {
  int nsegs;
  int run_no;       /* count runs */

  /* memory for pairwise distance calculations */
  distmem_t *distmem;
  struct genome_struct *g1, *g2;

  /* memory to work on copies of the genomes */
  struct genome_struct *genome_list;

  /* distance matrix, num_genomes x num_genomes */
  int *dmat;

  /* row provisionally being modified */
  int *dmat_r;

  /* list of which segment #'s have been tested at current best distance
   * options x=0,1,2,...,nsegs-1: flip segment x
   * options x=nsegs,nsegs+1,...,2*nsegs-1:
   *         flip segment x and segment x+1  (only do for certain x)
   */
  int *options;
  int num_options;


  int *signs;       /* signs[i] indicates sign of segment i:
		     *           0: -;  1: +
		     */
  int *bestsigns;   /* signs on one best-scoring instance */

  int *dist_counts; /* dist_counts[i] = # runs at distance i */
  int *dist_corr;   /* 2d array, nsegs x nsegs;
		     * dist_corr[nsegs*i + j]
		     *   = # runs at best d with sign of seg i = sign of seg j
		     *     (provided i < j);
		     *   = # runs at best d with sign of seg i = +1
		     *     (provided i == j);
		     *   = ?
		     *     (provided i > j)
		     */
  int best_dist;    /* best d so far */
  int max_dist_possible;


  /* formatting line-by-line report */
  int width_run;          /* width of run # */
  int width_run2;         /* width of run # in header column*/
  int width_score;        /* width of score */
  int width_score2;       /* width of score in header column*/
  int width_mat_entry;    /* width of matrix entry */
  int width_segnum;       /* width of segment # */
};


void default_segs(
		  int num_genes,
		  int num_chromosomes,
		  int num_genomes,
		  struct genome_struct *genome_list,
		  int *nsegs_ret,
		  int **seg_list_ret
		  );

int **default_weight_mat(int num_genomes);

int unsignedhc_mcdist_mat2(
			  int num_genes,
			  int num_chromosomes,
			  int num_genomes,
			  int circular,
			  int verbose,
			  struct genome_struct *genome_list_in,
			  int nsegs,
			  int *seg_list,
			  int **weight_matrix,
			  int num_iterations
			  );

void uhc_init_random_seed();

long uhc_random_range(long n);

void uhc_run(
	     int num_genes,
	     int num_chromosomes,
	     int num_genomes,
	     int circular,
	     int verbose,
	     struct genome_struct *genome_list_in,
	     int nsegs,
	     int *seg_list,
	     int **weight_matrix,

	     /* output */
	     struct uhc_mem *uhcmem);

int uhc_dist(
	     struct uhc_mem *uhcmem,
	     int i1,
	     int i2,
	     int num_genes,
	     int num_chromosomes,
	     int circular);



void uhc_summarize_results(
	     int num_genes,
	     int num_chromosomes,
	     int num_genomes,
	     int circular,
	     struct genome_struct *genome_list_in,
	     int nsegs,
	     int *seg_list,

	     /* output */
	     struct uhc_mem *uhcmem);


void calc_sign_groups(int *sign_groups,
		      struct uhc_mem *uhcmem);


void uhc_setsigns(
	     int num_genes,
	     int num_chromosomes,
	     int num_genomes,
	     struct genome_struct *genome_list_in,
	     int nsegs,
	     int *seg_list,
	     int *signs,
	     struct uhc_mem *uhcmem);


inline int max2(int a, int b)
{
  return (a<b) ? b : a;
}



int unsignedhc_mcdist(
		      int num_genes,
		      int num_chromosomes,
		      int circular,
		      int verbose,
		      struct genome_struct *genome_list,
		      int gindex1, int gindex2,
		      int nsegs,
		      int *seg_list,
		      int **weight_matrix,
		      int num_iterations
		      )
{
  int d;
  int num_genomes2 = 2;
  struct genome_struct *genome_list2
    = (struct genome_struct *) e_malloc(num_genomes2
					*sizeof(struct genome_struct),
					"unsignedhc_mcdist: genome_list");
  int **weight_matrix2 = (int **) NULL;

  /* TODO: nsegs is upper bound could compute exactly */
  int *seg_list2;

  int nsegs2 = 0;
  int i;
  int new_num;
  int *curseg, *curseg2;

  /* make subsetted list of segments */
  if (nsegs >= 0) {
    seg_list2 =
      (int *) e_malloc(3 * nsegs * sizeof(int),
		       "unsignedhc_mcdist: seg_list2");
  } else {
    seg_list2 = (int *) 0;
  }


  /* make list of just 2 genomes */
  memcpy(genome_list2, genome_list+gindex1, sizeof(struct genome_struct));
  memcpy(genome_list2+1, genome_list+gindex2, sizeof(struct genome_struct));

  /* copy segments in these genomes */
  if (nsegs < 0) {
    nsegs2 = nsegs;
  } else {
    nsegs2 = 0;
    curseg2 = seg_list2;
    for (i=0; i<nsegs; i++) {
      curseg = USEG(seg_list, i);
      if (curseg[0] == gindex1) {
	new_num = 0;
      } else if (curseg[0] == gindex2) {
	new_num = 1;
      } else {
	continue;
      }
      
      *curseg2++ = new_num;     /* change genome number */
      *curseg2++ = curseg[1];   /* copy offset */
      *curseg2++ = curseg[2];   /* copy length */
      nsegs2++;
    }
  }


  /* subsetted weight matrix */
  weight_matrix2 = default_weight_mat(2);
  if (weight_matrix != (int **) NULL) {
    weight_matrix2[0][1] = weight_matrix2[1][0] =
      weight_matrix[gindex1][gindex2];
  }


  d = unsignedhc_mcdist_mat(num_genes, num_chromosomes, num_genomes2,
			    circular, verbose,
			    genome_list2,
			    nsegs2, seg_list2,
			    weight_matrix2,
			    num_iterations);

  destroy_mat_int2d(weight_matrix2,2);
  free((void *) genome_list2);
  free((void *) seg_list2);

  return d;
}

int unsignedhc_mcdist_mat(
			  int num_genes,
			  int num_chromosomes,
			  int num_genomes,
			  int circular,
			  int verbose,
			  struct genome_struct *genome_list,
			  int nsegs,
			  int *seg_list,
			  int **weight_matrix,
			  int num_iterations
			  )
{
  int nsegs2;
  int *seg_list2;
  int **weight_matrix2;
  int d2;

#if 0
  if (nsegs == 0) {
    fprintf(outfile, "No segments!\n");
    return 0;
  }
#endif

  /* fill in default weight matrix if necessary */
  if (weight_matrix == (int **) NULL) {
    weight_matrix2 = default_weight_mat(num_genomes);
  } else {
    weight_matrix2 = weight_matrix;
  } 


  /* fill in default segment list if necessary */
  if (nsegs == -1) {
    default_segs(num_genes, num_chromosomes, num_genomes,
		 genome_list,
		 &nsegs2,
		 &seg_list2);
  } else {
    nsegs2 = nsegs;
    seg_list2 = seg_list;
  }

  d2 = unsignedhc_mcdist_mat2(
			      num_genes, num_chromosomes, num_genomes,
			      circular, verbose,
			      genome_list,
			      nsegs2, seg_list2,
			      weight_matrix2,
			      num_iterations
			      );

  /* if filled in defaults, destroy them */
  if (nsegs == -1) {
    free((void *) seg_list2);
  }

  if (weight_matrix == (int **) NULL) {
    destroy_mat_int2d(weight_matrix2, num_genomes);
  }

  return d2;
}


/* assign default segments:
 * genome 0 has no segments
 * genomes 1,2,...: each gene is its own segment, caps are not in segments
 */
void default_segs(
		  int num_genes,
		  int num_chromosomes,
		  int num_genomes,
		  struct genome_struct *genome_list,
		  int *nsegs_ret,
		  int **seg_list_ret
		  )
{
  int unichrom = (num_chromosomes == 0);
  int lowcap = num_genes - 2*num_chromosomes + 1;
  int nsegs = 0;
  int nsegs_alloc = (lowcap-1) * (num_genomes-1);
  int *seg_list = e_malloc(3 * nsegs_alloc * sizeof(int),
			   "default_segs");

  int genome_num, gnum;
  int *genes;


  *seg_list_ret = seg_list;

  for (genome_num = 1; genome_num < num_genomes; genome_num++) {
    genes = genome_list[genome_num].genes;
    for (gnum = 0; gnum < num_genes; gnum++) {

      /* unichrom case: every gene becomes its own segment
       * multichrom case: genes (not caps) become their own seg;
       *   exception: single gene chromosomes.
       */
      if (unichrom
	  || (genes[gnum] < lowcap
	      && !(genes[gnum-1]>=lowcap && genes[gnum+1]>=lowcap))) {

	*seg_list++ = genome_num;   /* genome number */
	*seg_list++ = gnum;         /* start position */
	*seg_list++ = 1;            /* length */
	nsegs++;                    /* count # segments */
      }
    }
  }
  *nsegs_ret = nsegs;
}

int **default_weight_mat(int num_genomes)
{
  int i,j;
  int **m = alloc_mat_int2d(num_genomes, num_genomes);
  for (i=0; i<num_genomes; i++) {
    for (j=0; j<num_genomes; j++) {
      m[i][j] =  (i != j);    /* 1 off diagonal, 0 on diagonal */
    }
  }
  return m;
}


/* Do the main algorithm */
int unsignedhc_mcdist_mat2(
			  int num_genes,
			  int num_chromosomes,
			  int num_genomes,
			  int circular,
			  int verbose,
			  struct genome_struct *genome_list_in,
			  int nsegs,
			  int *seg_list,
			  int **weight_matrix,
			  int num_iterations
			  )
{

  distmem_t distmem_mem;
  struct uhc_mem uhcmem_mem, *uhcmem=&uhcmem_mem;
  struct genome_struct uhcmem_g1, uhcmem_g2;

  /* extreme upper bound on max distance possible */
  /* TODO: improve for use with arbitrary weight matrices */
  int max_dist_possible = num_genomes*(num_genomes-1)/2 * (num_genes+1);

  int i,k,s1;
  int cur_genome_num;
  struct genome_struct *cur_genome = (struct genome_struct *) 0;
  int *curseg;
  int cur_start, cur_len;
  int num_options;

  /**********************************************************************
   * allocate memory
   **********************************************************************/

  uhcmem->max_dist_possible = max_dist_possible;
  uhcmem->best_dist = max_dist_possible;
  uhcmem->distmem = &distmem_mem;
  uhcmem->run_no = 0;

  uhcmem->genome_list =
    (struct genome_struct *) e_malloc(num_genomes*sizeof(struct genome_struct),

				      "unsignedhc_mcdist_mat2: genome_list");

  uhcmem->g1 = &uhcmem_g1;
  alloc_simple_genome(num_genes, &uhcmem_g1);
  uhcmem->g2 = &uhcmem_g2;
  alloc_simple_genome(num_genes, &uhcmem_g2);

  /* allocate space to work on a copy of all the genomes;
   * keep the originals intact
   */
  for (i=0; i<num_genomes; i++) {
    uhcmem->genome_list[i].gnamePtr = genome_list_in[i].gnamePtr;
    uhcmem->genome_list[i].genome_num = genome_list_in[i].genome_num;
    uhcmem->genome_list[i].encoding = genome_list_in[i].encoding;
    uhcmem->genome_list[i].genes =
      (int *) e_malloc(num_genes * sizeof(int), "unsignedhc_mcdist_mat2: genes");
  }


  /* allocate memory for distance computations */
  mcdist_allocmem(num_genes, num_chromosomes, uhcmem->distmem);


  /* allocate memory for distance matrix */
  uhcmem->dmat = (int *) e_malloc(num_genomes*num_genomes*sizeof(int),
				  "unsignedhc_mcdist_mat2: dmat");
  /* allocate memory for modified row of distance matrix */
  uhcmem->dmat_r = (int *) e_malloc(num_genomes*sizeof(int),
				  "unsignedhc_mcdist_mat2: dmat_r");

  /* keep track of which flips have been tested and which have not */
  uhcmem->options = (int *) e_malloc(2 * nsegs * sizeof(int),
				     "unsignedhc_mcdist_mat2: options");

  /* allocate memory for statistical summary of results */

  /* # runs with each score */
  uhcmem->dist_counts = (int *) e_calloc(max_dist_possible+1, sizeof(int),
				 "unsignedhc_mcdist_mat2: dist_counts");

  /* correlations */
  uhcmem->dist_corr = (int *) e_calloc(nsegs * nsegs, sizeof(int),
			       "unsignedhc_mcdist_mat2: dist_corr");
  uhcmem->nsegs = nsegs;


  uhcmem->signs = (int *) e_malloc(nsegs * sizeof(int),
				   "unsignedhc_mcdist_mat2: signs");
  uhcmem->bestsigns = (int *) e_malloc(nsegs * sizeof(int),
				       "unsignedhc_mcdist_mat2: bestsigns");


  uhcmem->width_run = num_digits(num_iterations);
  uhcmem->width_run2 = max2(uhcmem->width_run, 5);
  uhcmem->width_score = num_digits(max_dist_possible);
  uhcmem->width_score2 = max2(uhcmem->width_score, 5);
  uhcmem->width_mat_entry = num_digits(num_genes+1);
  uhcmem->width_segnum = num_digits(nsegs);


  uhc_init_random_seed();

  /**********************************************************************
   * initialize list of allowed options
   * 1. all segments individually
   * 2. "antistrips" (just to increase the probability of finding them)
   **********************************************************************/

  /* individual segments */
  num_options = 0;
  for (s1=0; s1<nsegs; s1++) {
    uhcmem->options[num_options++] = s1;
  }

  /* antistrips */
  for (s1=0; s1<nsegs-1; s1++) {
    curseg = USEG(seg_list,s1);
    /* are segs s1 & s1+1 adjacent? skip if not. */
    if (curseg[0] != curseg[3]                /* genome numbers */
	|| curseg[1]+curseg[2] != curseg[4]   /* adjacent positions */
	) {
      /* not adjacent */
      continue;
    }
    /* are s1 & s1-1 adjacent? skip if so (bigger than 2-strip) */
    if (s1>0
	&& curseg[0]==curseg[-3]
	&& curseg[1]==curseg[-2]+curseg[-1]) {
      /* adjacent, 3+ strip, skip */
      continue;
    }
    uhcmem->options[num_options++] = s1 + nsegs;
  }
  uhcmem->num_options = num_options;

  /**********************************************************************
   * print list of segments
   **********************************************************************/

  if (verbose) {
    fprintf(outfile,"\nSegments:\n");
    cur_genome_num = -1;
    for (s1 = 0; s1 < nsegs; s1++) {
      curseg = USEG(seg_list, s1);

      /* print genome name if switched genomes */
      if (curseg[0] != cur_genome_num) {
	cur_genome_num = curseg[0];
	cur_genome = &genome_list_in[cur_genome_num];
	if (cur_genome->gnamePtr != (char *) 0) {
	  fprintf(outfile, ">%s\n", cur_genome->gnamePtr);
	} else {
	  fprintf(outfile, ">genome%d\n", cur_genome_num + 1);
	}
      }
      
      cur_start = curseg[1];
      cur_len = curseg[2];

      fprintf(outfile, "S%-*d  [", uhcmem->width_segnum, s1);
      for (k=0; k<cur_len; k++) {
	fprintf(outfile, " %d", cur_genome->genes[cur_start+k]);
      }
      fprintf(outfile," ]\n");
    }
  }

  /**********************************************************************
   * do lots of iterations
   **********************************************************************/

  if (verbose) {
    fprintf(outfile, "\nTrials:\n");
    fprintf(outfile, "%-*s %-*s  matrix %*s signs\n",
	    uhcmem->width_run2, "run",
	    uhcmem->width_score2, "score",
	    ((uhcmem->width_mat_entry+1)*num_genomes + 2)*num_genomes - 5, "");
  }

  for (i=0; i<num_iterations; i++) {
    uhc_run(num_genes,num_chromosomes,num_genomes,
	    circular,verbose,
	    genome_list_in,
	    nsegs, seg_list,
	    weight_matrix,
	    uhcmem);
  }

  /**********************************************************************
   * report on results
   **********************************************************************/

  if (verbose) {
    uhc_summarize_results(num_genes, num_chromosomes, num_genomes,
			  circular,
			  genome_list_in,
			  nsegs, seg_list,
			  uhcmem);
  }

  /**********************************************************************
   * cleanup memory
   **********************************************************************/

  free((void *) uhcmem->bestsigns);
  free((void *) uhcmem->signs);
  free((void *) uhcmem->options);
  free((void *) uhcmem->dmat_r);
  free((void *) uhcmem->dmat);
  free((void *) uhcmem->dist_corr);
  free((void *) uhcmem->dist_counts);

  /* free memory for distance computations */
  mcdist_freemem(uhcmem->distmem);

  free_simple_genome((void *) uhcmem->g2);
  free_simple_genome((void *) uhcmem->g1);

  /* free memory for copy of genomes */
  for (i=0; i<num_genomes; i++) {
    free((void *) uhcmem->genome_list[i].genes);
  }
  free((void *) uhcmem->genome_list);

  return uhcmem->best_dist;
}

void uhc_init_random_seed()
{
  srandom(1);
}

/* return random integer in range 0,1,2,...,n-1
 * There is a slight bias because random() returns one of 2^31 values,
 * and 2^31 is not usually divisible by n
 */
long uhc_random_range(long n)
{
  return random() % n;
}

void uhc_run(
	     int num_genes,
	     int num_chromosomes,
	     int num_genomes,
	     int circular,
	     int verbose,
	     struct genome_struct *genome_list_in,
	     int nsegs,
	     int *seg_list,
	     int **weight_matrix,

	     /* output */
	     struct uhc_mem *uhcmem)
{
  int *signs = uhcmem->signs;
  int *dmat = uhcmem->dmat;
  int *dmat_r = uhcmem->dmat_r;
  int *options = uhcmem->options;
  int noptions;

  int nrandbits0 = 31;
  int nrandbits = 0;
  long randbits = 0;

  int i,j;
  int *curseg;
  int score, dscore1, dscore2;
  int i1, i2, s1, s2;
  int opt_index, opt_no;
  int seg_no;

  int w;

  uhcmem->run_no++;

  /**********************************************************************
   * choose random initial signs
   **********************************************************************/

  for (i=0; i<nsegs; i++) {
    if (nrandbits == 0) {
      randbits = random();
      nrandbits = nrandbits0;
    }
    signs[i] = randbits & 1;
    randbits >>= 1;
    nrandbits--;
  }

  /**********************************************************************
   * copy genomes and do the sign flips
   **********************************************************************/

  uhc_setsigns(num_genes, num_chromosomes, num_genomes,
	       genome_list_in,
	       nsegs, seg_list, signs,
	       uhcmem);

  /**********************************************************************
   * initialize distance matrix
   **********************************************************************/

  score = 0;
  for (i1=0; i1<num_genomes; i1++) {
    dmat[num_genomes*i1 + i1] = 0;
    for (i2=0; i2<i1; i2++) {
      w = weight_matrix[i1][i2];
      if (w == 0) continue;

      dmat[num_genomes*i1 + i2] =
	dmat[num_genomes*i2 + i1] =
	uhc_dist(uhcmem, i1, i2,
		 num_genes, num_chromosomes, circular);

      score += w * dmat[num_genomes*i1 + i2];
    }
  }

  /**********************************************************************
   * hill climb init
   **********************************************************************/

#if 0
  /* initialize list of options */
  noptions = nsegs;
  for (j=0; j<noptions; j++) {
    options[j] = j;
  }
#endif
  noptions = uhcmem->num_options;

  /**********************************************************************
   * hill climb
   **********************************************************************/

  while (noptions>0) {
    opt_index = uhc_random_range(noptions);
    opt_no = options[opt_index];
    seg_no = (opt_no >= nsegs) ? opt_no-nsegs : opt_no;
    
    curseg = USEG(seg_list, seg_no);
    i1 = curseg[0];

    reverse_in_place_genes(uhcmem->genome_list[i1].genes,
			   curseg[1],
			   curseg[1] + curseg[2] - 1);
    if (opt_no >= nsegs) {
      reverse_in_place_genes(uhcmem->genome_list[i1].genes,
			     curseg[4],
			     curseg[4] + curseg[5] - 1);
    }

    dscore1 = 0; /* score before flip */
    dscore2 = 0; /* score after flip */
    for (i2 = 0; i2 < num_genomes; i2++) {
      if (i1==i2) {
	dmat_r[i2] = 0;
	continue;
      }

      w = weight_matrix[i1][i2];
      if (w == 0) continue;

      dmat_r[i2] = uhc_dist(uhcmem, i1, i2,
			    num_genes, num_chromosomes, circular);
      dscore2 += w * dmat_r[i2];
      dscore1 += w * dmat[num_genomes*i1 + i2];
    }

    /* did score improve, tie, or worsen? */
    if (dscore2 >= dscore1) {
      /* tied or worsened, undo flip */
      reverse_in_place_genes(uhcmem->genome_list[i1].genes,
			     curseg[1],
			     curseg[1] + curseg[2] - 1);
      if (opt_no >= nsegs) {
	reverse_in_place_genes(uhcmem->genome_list[i1].genes,
			       curseg[4],
			       curseg[4] + curseg[5] - 1);
      }
    } else {
      /* improved, keep flip */
      signs[seg_no] = !signs[seg_no];   /* record sign change */
      if (opt_no >= nsegs) {
	signs[seg_no+1] = !signs[seg_no+1];
      }

      /* update score */
      score = score - dscore1 + dscore2;

      /* update distance matrix */
      for (i2 = 0; i2 < num_genomes; i2++) {
	dmat[num_genomes*i1 + i2] =
	  dmat[num_genomes*i2 + i1] = dmat_r[i2];
      }

      /* re-open all other hill-climb options */
#if 0
      noptions = nsegs;
#endif
      noptions = uhcmem->num_options;
    }

    /* swap this option with last possible option,
     * so that remaining untried options are all listed first
     */

    noptions--;
    options[opt_index] = options[noptions];
    options[noptions] = opt_index;
  }

  /**********************************************************************
   * output result
   **********************************************************************/

  if (verbose) {
    /* run# score matrix signs */
    fprintf(outfile,"%*d %*d: [",
	    uhcmem->width_run2, uhcmem->run_no,
	    uhcmem->width_score2, score);
    for (i1=0; i1<num_genomes; i1++) {
      if (i1 > 0) {
	fprintf(outfile, " ;");
      }
      for (i2=0; i2<num_genomes; i2++) {
	fprintf(outfile," %*d",
		uhcmem->width_mat_entry,
		dmat[i1*num_genomes + i2]);
      }
    }
    fprintf(outfile, " ]  ");

    for (j=0; j<nsegs; j++) {
      putc(signs[j] ? '+' : '-', outfile);
    }
    fprintf(outfile, "\n");
  }

  /**********************************************************************
   * update statistical summary
   **********************************************************************/

  /* # at each score */
  if (score < uhcmem->max_dist_possible && score >= 0)
    uhcmem->dist_counts[score]++;
  else
    uhcmem->dist_counts[uhcmem->max_dist_possible]++;  /* under/overflow */

  /* new best score */
  if (uhcmem->best_dist > score) {
    /* new best score */
    uhcmem->best_dist = score;

    /* reset all counts */
    memset(uhcmem->dist_corr, 0, nsegs*nsegs*sizeof(int));

    /* store best results */
    memcpy(uhcmem->bestsigns, uhcmem->signs, nsegs*sizeof(int));
  }

  /* best score */
  if (uhcmem->best_dist == score) {
    /* individual sign counts */
    for (s1=0; s1<nsegs; s1++) {
      uhcmem->dist_corr[s1*nsegs+s1] += signs[s1];
    }

    /* correlations */
    for (s1=1; s1<nsegs; s1++) {
      for (s2=0; s2<s1; s2++) {
	if (signs[s1]==signs[s2]) {
	  uhcmem->dist_corr[s2*nsegs+s1]++;
	}
      }
    }
  }
}


int uhc_dist(
	     struct uhc_mem *uhcmem,
	     int i1,
	     int i2,
	     int num_genes,
	     int num_chromosomes,
	     int circular)
{
  int d;
  struct genome_struct *g1, *g2;

  /* some of the distance routines modify the gene data, so we must
   * work on copies
   */

  if (num_chromosomes == 0) {
    /* unichromosomal does not modify the input genomes */
    g1 = &uhcmem->genome_list[i1];
    g2 = &uhcmem->genome_list[i2];
  } else {
    /* multichromosomal does modify the genomes */
    g1 = uhcmem->g1;
    g2 = uhcmem->g2;
    copy_genes(uhcmem->genome_list[i1].genes, g1->genes, num_genes);
    copy_genes(uhcmem->genome_list[i2].genes, g2->genes, num_genes);
  }

  if (num_chromosomes == 0) {
    if (circular) {
      /* unichromosomal circular */
      d = invdist_circular(g1, g2, num_genes, uhcmem->distmem);
    } else {
      /* unichromosomal linear */
      d = invdist_noncircular(g1, g2, 0, num_genes, uhcmem->distmem);
    }
  } else {
    /* multichromosomal */
    d = mcdist_noncircular(g1, g2,
			   num_genes, num_chromosomes,
			   uhcmem->distmem,
			   (graphstats_t *) NULL);
  }

  return d;
}


void uhc_summarize_results(
	     int num_genes,
	     int num_chromosomes,
	     int num_genomes,
	     int circular,
	     struct genome_struct *genome_list_in,
	     int nsegs,
	     int *seg_list,

	     /* output */
	     struct uhc_mem *uhcmem)
{
  int d;

  int c, count;
  int s1, s2;
  int bestd, num_best;
  int k;

  int *curseg;
  int cur_start, cur_len;
  int cur_genome_num;
  struct genome_struct *cur_genome = (struct genome_struct *) 0;
  int genome_num;

  int *sign_groups;

  int i1,i2,pad,max_entry,entry_width;
  int *dmat = uhcmem->dmat;

  /**********************************************************************
   * Calc sign groups
   **********************************************************************/

  sign_groups = (int *) e_malloc(nsegs * sizeof(int),
				 "uhc_summarize_results: sign_groups");
  calc_sign_groups(sign_groups,
		   uhcmem);



  /**********************************************************************
   * # runs at each distance
   **********************************************************************/

  fprintf(outfile, "\n# trials giving each score:\n");
  fprintf(outfile, "%*s  # times\n",
	  uhcmem->width_score2, "score");
  for (d=0; d < uhcmem->max_dist_possible; d++) {
    if (uhcmem->dist_counts[d] > 0) {
      fprintf(outfile, "%*d  %d\n",
	      uhcmem->width_score2, d,
	      uhcmem->dist_counts[d]);
    }
  }
  if (uhcmem->dist_counts[uhcmem->max_dist_possible] != 0) {
    /* TODO: improve handling of this for arbitrary weight matrices */
    fprintf(outfile, "%*s  %d\n",
	    uhcmem->width_score2, "other",
	    uhcmem->dist_counts[uhcmem->max_dist_possible]);
  }


  /**********************************************************************
   * at best dist: sign pattern
   **********************************************************************/

  bestd = uhcmem->best_dist;
  num_best = uhcmem->dist_counts[bestd];

  fprintf(outfile, "\nBest score: %d\n", bestd);

  fprintf(outfile, "\nSimple sign pattern:  ");

  for (s1=0; s1<nsegs; s1++) {
    c = '?';
    count = sign_groups[s1];
    if (count == -1-nsegs) {
      c = '-';
    } else if (count == nsegs) {
      c = '+';
    }
    putc(c, outfile);
  }
  fprintf(outfile, "\n");
    
  /**********************************************************************
   * at best dist: perfect correlation groups
   **********************************************************************/

  fprintf(outfile, "\nGrouped sign pattern:  ");
  for (s1=0; s1<nsegs; s1++) {
    c = '?';
    count = sign_groups[s1];
    if (count == -1-nsegs) {
      fprintf(outfile, "-");
    } else if (count == nsegs) {
      fprintf(outfile, "+");
    } else if (count >= 0) {
      fprintf(outfile, "[S%d]", count);
    } else {
      fprintf(outfile, "[-S%d]", -1-count);
    }
  }
  fprintf(outfile, "\n");

  /**********************************************************************
   * at best dist: sign counts
   **********************************************************************/

  fprintf(outfile,"\nSign counts and perfect correlations:\n");
  cur_genome_num = -1;
  for (s1 = 0; s1 < nsegs; s1++) {
    curseg = USEG(seg_list, s1);

    /* print genome name if switched genomes */
    if (curseg[0] != cur_genome_num) {
      cur_genome_num = curseg[0];
      cur_genome = &genome_list_in[cur_genome_num];
      if (cur_genome->gnamePtr != (char *) 0) {
	fprintf(outfile, ">%s\n", cur_genome->gnamePtr);
      } else {
	fprintf(outfile, ">genome%d\n", cur_genome_num + 1);
      }
    }
      
    cur_start = curseg[1];
    cur_len = curseg[2];

    fprintf(outfile, "S%-*d  [", uhcmem->width_segnum, s1);
    for (k=0; k<cur_len; k++) {
      fprintf(outfile, " %d", cur_genome->genes[cur_start+k]);
    }
    fprintf(outfile," ]   ");

    count = uhcmem->dist_corr[s1*nsegs+s1];
    fprintf(outfile,
	    "%d+  %d-  ",
	    count, num_best-count);

    count = sign_groups[s1];
    if (count == -1-nsegs) {
      fprintf(outfile, "sign -");
    } else if (count == nsegs) {
      fprintf(outfile, "sign +");
    } else if (count == s1) {
    } else if (count >= 0) {
      fprintf(outfile, "sign same as S%d", count);
    } else if (count < 0) {
      fprintf(outfile, "sign same as -S%d", -1-count);
    }

    fprintf(outfile, "\n");
  }

  /**********************************************************************
   * at best dist: detailed correlations
   **********************************************************************/

  fprintf(outfile,"\nCorrelations: Number of times segments have same sign:\n");
  for (s1=0; s1<nsegs; s1++) {
    fprintf(outfile, "S%-*d  ", uhcmem->width_segnum, s1);
    for (s2=0; s2<nsegs; s2++) {
      if (s1 == s2) {
	count = num_best;
      } else if (s2<s1) {
	count = uhcmem->dist_corr[s2*nsegs+s1];
      } else {
	count = uhcmem->dist_corr[s1*nsegs+s2];
      }
      fprintf(outfile, " %*d", uhcmem->width_run, count);
    }
    fprintf(outfile, "\n");
  }

  /**********************************************************************
   * at best dist: one particular signage
   **********************************************************************/

  fprintf(outfile,"\nA best scoring solution:\n");

  uhc_setsigns(num_genes, num_chromosomes, num_genomes,
	       genome_list_in,
	       nsegs, seg_list, uhcmem->bestsigns,
	       uhcmem);
  for (genome_num=0; genome_num<num_genomes; genome_num++) {
    if (num_chromosomes == 0) {
      print_genome_unichrom(uhcmem->genome_list + genome_num,
			    num_genes);
    } else {
      print_genome_multichrom_3(uhcmem->genome_list + genome_num,
				num_genes,
				num_chromosomes,
				0, /* use min width for each gene */
				1, /* 1 chromo per line */
				0, /* omit caps, chromo delim = $ */
				1  /* show genome name */
				);
    }
  } 

  /**********************************************************************
   * and the pairwise matrix for that signage
   **********************************************************************/

  max_entry = 0;
  for (i1=0; i1<num_genomes; i1++) {
    dmat[num_genomes*i1 + i1] = 0;
    for (i2=0; i2<i1; i2++) {
      d =
	dmat[num_genomes*i1 + i2] =
	dmat[num_genomes*i2 + i1] =
	uhc_dist(uhcmem, i1, i2,
		 num_genes, num_chromosomes, circular);
      if (d > max_entry) { max_entry = d; }
    }
  }

  /* TODO: should integrate this with the other matrix printing routines */

#define GNAME_WIDTH 20

  fprintf(outfile,"\nDistance matrix for that specific solution:\n");

  entry_width = num_digits(max_entry);
  for (i1=0 ; i1<num_genomes ; i1++) {

    /* print genome name, truncating or padding if necessary */
    for(i2=0, pad=FALSE; i2<GNAME_WIDTH; i2++) {
      if (!pad) {
	if (genome_list_in[i1].gnamePtr[i2] != '\0')
	  fputc(genome_list_in[i1].gnamePtr[i2], outfile);
	else
	  pad = TRUE;
      }
      if (pad) fputc(' ', outfile);
    }
    fputc('\t', outfile);
    
    /* fprintf(outfile,"%s\t",genome_list_in[i1].gnamePtr); */

    for (i2=0; i2<i1 ; i2++)
      fprintf(outfile, "%*d ", entry_width, dmat[num_genomes*i2 + i1]);
    fprintf(outfile, "%*d ", entry_width, 0);
    for (i2=i1+1 ; i2<num_genomes ; i2++) {
      fprintf(outfile,"%*d ", entry_width, dmat[num_genomes*i2 + i1]);
    }
    fprintf(outfile,"\n");
  }
  fprintf(outfile,"\n");






}

/* sign_groups[i] = sign group of segment i
 *                  group nsegs: constant sign +
 *                  group -1-nsegs: constant sign -
 *                  group j>=0: same as segment Sj
 *                  group j<0: same as segment Sk, k=-1-j
 *                     where j (or k) is minimum possible
 */

void calc_sign_groups(int *sign_groups,
		      struct uhc_mem *uhcmem)
{
  int nsegs = uhcmem->nsegs;
  int bestd = uhcmem->best_dist;
  int num_best = uhcmem->dist_counts[bestd];
  int count1, count2;
  int s1, s2;

  for (s1=0; s1<nsegs; s1++) {
    count1 = uhcmem->dist_corr[s1*nsegs+s1];
    if (count1 == 0) {
      sign_groups[s1] = -1-nsegs;   /* always - */
    } else if (count1 == num_best) {
      sign_groups[s1] = nsegs;      /* always + */
    } else {
      /* find least segment # to whose sign it is perfectly correlated */
      sign_groups[s1] = s1;         /* default: no prior correlated segment */
      for (s2=0; s2<s1; s2++) {
	count2 = uhcmem->dist_corr[s2*nsegs+s1];
	if (count2 == 0) {
	  sign_groups[s1] = -1-s2;  /* segs s1 & s2 always have opp. sign */
	  break;
	} else if (count2 == num_best) {
	  sign_groups[s1] = s2;
	  break;
	}
      }
    }
  }
}

void uhc_setsigns(
	     int num_genes,
	     int num_chromosomes,
	     int num_genomes,
	     struct genome_struct *genome_list_in,
	     int nsegs,
	     int *seg_list,
	     int *signs,
	     struct uhc_mem *uhcmem)
{
  int genome_num, i;
  int *curseg;

  /* copy genomes */
  for (genome_num = 0; genome_num < num_genomes; genome_num++) {
    copy_genes(genome_list_in[genome_num].genes,
	       uhcmem->genome_list[genome_num].genes,
	       num_genes);
  }

  /* initial sign flips */
  for (i=0; i<nsegs; i++) {
    if (!signs[i]) {
      curseg = USEG(seg_list,i);
      reverse_in_place_genes(uhcmem->genome_list[curseg[0]].genes,
			     curseg[1],
			     curseg[1] + curseg[2] - 1);
    }
  }
}
