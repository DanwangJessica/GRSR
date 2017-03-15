/* gs_grimm.c
 *    Interface to GRIMM.
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

/* Last modified on Tue Aug 1, 2006, by Glenn Tesler
 */

#include <stddef.h>
#include "ckalloc.h"
#include "mcstructs.h"
#include "mcrdist.h"
#include "gs_grimm.h"


/* Allocate space
 * It's not mandatory to allocate all structures
 */
void gs_grimm_alloc(int num_chromosomes,
		    int num_genes,
		    int num_genomes,

		    distmem_t **distmem_ret,
		    graphstats_t ***statsmatrix_ret,
		    struct genome_struct **genome_list_ret)
{
  distmem_t *distmem;
  graphstats_t **statsmatrix;
  struct genome_struct *genome_list;
  int n = num_genes + 2*num_chromosomes;  /* space for genes & caps */
  int i;

  if (genome_list_ret != (struct genome_struct **) NULL) {
    genome_list =
      (struct genome_struct *)
      ckalloc(num_genomes, sizeof(struct genome_struct));

    /* allocate space for each genome */
    for (i=0; i<num_genomes; i++) {
      genome_list[i].gnamePtr = (char *) NULL;
      genome_list[i].encoding = (char *) NULL;
      genome_list[i].genes = (int *) ckalloc(n, sizeof(int));
    }

    *genome_list_ret = genome_list;
  }

  if (distmem_ret != (distmem_t **) NULL) {
    distmem = (distmem_t *) ckalloc(1, sizeof(distmem_t));
    mcdist_allocmem(n, num_chromosomes, distmem);
    *distmem_ret = distmem;
  }

  if (statsmatrix_ret != (graphstats_t ***) NULL) {
    statsmatrix =
      (graphstats_t **) ckalloc(num_genomes, sizeof(graphstats_t *));
    for (i=0; i<num_genomes; i++) {
      statsmatrix[i] =
	(graphstats_t *) ckalloc(num_genomes, sizeof(graphstats_t));
    }

    *statsmatrix_ret = statsmatrix;
  }
}

void gs_grimm_destroy(int num_chromosomes,
		      int num_genes,
		      int num_genomes,
		      distmem_t *distmem,
		      graphstats_t **statsmatrix,
		      struct genome_struct *genome_list)
{
  int i;

  if (statsmatrix != (graphstats_t **) NULL) {
    for (i=0; i<num_genomes; i++) {
      free((void *) statsmatrix[i]);
    }
    free((void *) statsmatrix);
  }

  if (distmem != (distmem_t *) NULL) {
    mcdist_freemem(distmem);
    free((void *) distmem);
  }

  if (genome_list != (struct genome_struct *) NULL) {
    /* free space for each genome */
    for (i=0; i<num_genomes; i++) {
      if (genome_list[i].gnamePtr != (char *) NULL) {
	free((void *) genome_list[i].gnamePtr);
      }
      if (genome_list[i].encoding != (char *) NULL) {
	free((void *) genome_list[i].encoding);
      }
      if (genome_list[i].genes != (int *) NULL) {
	free((void *) genome_list[i].genes);
      }
    }
    free((void *) genome_list);
  }
}






