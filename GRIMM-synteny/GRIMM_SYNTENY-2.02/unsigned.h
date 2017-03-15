/* unsigned.h
 *    Exact algorithm to find a most parsimonious assignment of signs
 *    between two unsigned genomes.
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

#ifndef UNSIGNED_H
#define UNSIGNED_H

int unsigned_mcdist_nomem(struct genome_struct *g1,
			  struct genome_struct *g2,
			  int num_genes, int num_chromosomes,
			  int circular,
			  int verbose,
			  int *dist_error);

int unsigned_mcdist(struct genome_struct *g1, struct genome_struct *g2,
		    int num_genes,
    		    int num_chromosomes,
		    int circular,
		    distmem_t *distmem,
		    int verbose,
		    int *dist_error);


void set_unsigneddist_matrix(int **distmatrix, struct genome_struct *genomes,
			     int num_genes, int num_chromosomes,
			     int num_genomes,
			     int circular,
			     distmem_t *distmem,
			     int *dist_error);


/*
  used for determining the optimal "spin" (signs)
  of unsigned permutations
*/

typedef struct {
  int num_genes, num_chromosomes;

  int num_sing;         /* # of singletons */
  int num_doub;         /* # of 2-strips */
  int *singletons;      /* array of locations of singletons w/in g2 */

  int *doubletons;      /* array of locations of 2-strips in g2;
			   record positions of left member */

  struct genome_struct *g1; /* reference genome Gamma */
  struct genome_struct *g2; /* destination genome Pi */
  struct genome_struct
   loop_g2,            /* Pi as we iterate over it */
   superspin_g2,       /* modifications to that to get a superspin */
#if 0
   best_g1,            /* best signed & capped Gamma found so far */
#endif
   best_g2;            /* best signed & capped Pi found so far */

  int best_d;           /* lowest distance found so far */

  distmem_t *distmem;      /* memory for distance computations */
  
} spins_t;

#endif /* UNSIGNED_H */
