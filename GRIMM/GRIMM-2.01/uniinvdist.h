/* uniinvdist.h
 *    1-line description
 *
 * Copyright (C) 2001-2006 The Regents of the University of California
 * by Glenn Tesler
 *
 * Contains code from invdist.h in GRAPPA 1.02
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

/* Last modified on Tue Aug 1, 2006, by Glenn Tesler
 */

/* Modified from GRAPPA 1.02: invdist.h */

#ifndef UNIINVDIST_H
#define UNIINVDIST_H

#include <sys/time.h>
/*#include "structs.h"*/
#include "mcstructs.h"


extern double time_linear;

void calc_invmatrix(struct genome_struct *genomes, int num_genes, int num_genomes,
		    distmem_t *distmem, int CIRCULAR);

void setinvmatrix(int **distmatrix,struct genome_struct *genomes,
		  int num_genes,int num_genomes,distmem_t *distmem,int CIRCULAR);

int invdist_noncircular_G(struct genome_struct *g1,
			  struct genome_struct *g2,
			  int offset, int num_genes,
			  distmem_t *distmem,
			  graph_t *G,
			  pathcounts_t *pcounts_G,
			  graphstats_t *graphstats);
int invdist_noncircular(struct genome_struct *g1, struct genome_struct *g2,
			int offset, int num_genes, distmem_t *distmem);
int invdist_noncircular_v(struct genome_struct *g1, struct genome_struct *g2,
			  int offset, int num_genes, distmem_t *distmem,
			  graphstats_t *graphstats);
int invdist_circular(struct genome_struct *g1, struct genome_struct *g2,
		     int num_genes, distmem_t *distmem);

int invdist_noncircular_nomem(struct genome_struct *g1,
			      struct genome_struct *g2,
			      int offset, int num_genes);
int invdist_noncircular_nomem_v(struct genome_struct *g1,
				struct genome_struct *g2,
				int offset, int num_genes,
				graphstats_t *graphstats);
int invdist_circular_nomem(struct genome_struct *g1,
			   struct genome_struct *g2,
			   int num_genes);

int calculate_offset(struct genome_struct *g1, struct genome_struct *g2,
		     int num_genes);

#endif /* UNIINVDIST_H */
