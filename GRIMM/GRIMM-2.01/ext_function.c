/* ext_function.c
 *    Template for new feature in GRIMM.  For testing only.
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

#include <stdio.h>
#include <stdlib.h>
#include "mcstructs.h"
#include "uniinvdist.h"
#include "mcrdist.h"
#include "unsigned.h"
#include "e_malloc.h"
#include "texgraph.h"
#include "ext_function.h"

/* sorry, I'm not going to change the header files now just for this */
extern
void print_full_distmatrix(int **distmatrix,int num_genomes, 
			   struct genome_struct *genome_list);


/* example of an additional function in MCDIST
 *   mcdist -f genome_file -X "args" ...
 *
 * -X "args"
 * produces the call
 *   do_fn_external(num_genes, num_chromosomes, num_genomes,
 *                  circular, unsigned_dist, verbose,
 *                  genome_list,
 *                  ext_args)
 *
 * where genome_file contains:
 * num_genomes:       number of genomes in genome_file
 * num_chromosomes:   number of chromosomes per genome
 *                    padded with null chromosomes if necessary
 *                    0 = unichromosomal mode
 * num_genes:         number of "genes" per genome
 *                    = # real markers + 2 * number of chromosomes
 * circular:          TRUE if -C option specified
 * unsigned_dist:     TRUE if -u option specified
 * verbose:           TRUE if -v option specified
 * genome_list:       array of genomes read in from genome_file
 * ext_args:          string "args" given with -X option
 */

/* For now, I just adapted code from do_fn_distmatrix, setinvmatrix,
 * setmcdistmatrix
 */

void do_fn_external(int num_genes, int num_chromosomes, int num_genomes,
		    int circular,
		    int unsigned_dist,
		    int verbose,
		    struct genome_struct *genome_list,
		    char *ext_args)
{
  int i, j;
  int **distmatrix;
  /*  int dist_error; */
  distmem_t *distmem;

  int aborted=0;

  fprintf(outfile,"\nExternal function %s\n", ext_args);

  if (verbose && !unsigned_dist) {
    /*
        do_fn_distmatrix_v(num_genes, num_chromosomes, num_genomes,
    		       circular, unsigned_dist,
    		       genome_list);
    */
    fprintf(outfile,"verbose, signed distance.  Code deleted.\n");
    return;
  }


  fprintf(outfile,"\nDistance Matrix:\n");

  /* allocate space for matrix */
  distmatrix = (int **) e_malloc((num_genomes)*sizeof(int *), "distmatrix");

  for (i=0; i < num_genomes; i++) {
    distmatrix[i] = (int *) e_malloc((num_genomes)*sizeof(int), "distmatrix");
  }

  /* allocate memory for distance computations */
  distmem = (distmem_t *) e_malloc(sizeof(distmem_t), "distmem");
  mcdist_allocmem(num_genes,num_chromosomes,distmem);


  /* compute the appropriate matrix */
  if (unsigned_dist) {
    fprintf(outfile,"Unsigned distance.  Code deleted.\n");
    aborted++;

#if 0
          set_unsigneddist_matrix(distmatrix, genome_list,
    			      num_genes, num_chromosomes, num_genomes,
    			      circular,
    			      distmem,
    			      &dist_error);
          if (dist_error) {
    	fprintf(outfile,"WARNING: These unsigned permutations are too complex;\n");
    	fprintf(outfile,"some answers are approximated.\n");
          }
#endif /* 0 */

  } else if (num_chromosomes == 0) {
    /* unichromosomal data */ 
#if 0
        setinvmatrix(distmatrix,genome_list,
    		 num_genes,num_genomes,
    		 distmem,
    
    		 /* circular already converted to equiv linear problem */
    		 FALSE);
#endif /* 0 */
    /* inline the code to demonstrate explicitly how distmem is used */

    for (i=0 ; i<num_genomes ; i++) {
      distmatrix[i][i] = 0;
      for (j=i+1 ; j<num_genomes ; j++) {
	if (circular)
	  distmatrix[i][j] = distmatrix[j][i] =
	    invdist_circular(genome_list+i,genome_list+j,
			     num_genes,distmem);
	else 
	  distmatrix[i][j] = distmatrix[j][i] =
	    invdist_noncircular(genome_list+i,genome_list+j,
				0,num_genes,distmem);
      }
    }



  } else {
#if 0
        /* multichromosomal data */
        setmcdistmatrix(distmatrix,genome_list,
    		    num_genes,num_chromosomes,num_genomes,
    		    distmem);
#endif /* 0 */

    /* inline the code to demonstrate explicitly how distmem is used */

    for (i=0 ; i<num_genomes ; i++) {
      distmatrix[i][i] = 0;
      for (j=i+1 ; j<num_genomes ; j++) {
	distmatrix[i][j] = distmatrix[j][i] =
	  mcdist_noncircular(genome_list+i,genome_list+j,
			     num_genes,num_chromosomes,distmem,
			     (graphstats_t *)NULL); /* not verbose */
      }
    }

  }

  if (!aborted) {
    /* display the matrix */
    print_full_distmatrix(distmatrix,num_genomes,genome_list);
  }

  /* free the memory for distance computations */
  mcdist_freemem(distmem);   /* frees memory pointed to by fields of distmem */
  free(distmem);             /* frees the distmem structure itself */

  /* free memory for the distance matrix */
  for (i=num_genomes-1; i >=0 ; i--)
    free(distmatrix[i]);
  free(distmatrix);
}
