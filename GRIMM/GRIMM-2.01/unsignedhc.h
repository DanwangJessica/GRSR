/* unsignedhc.h
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


#ifndef UNSIGNEDHC_H
#define UNSIGNEDHC_H


#define USEG(seg,i) (&((seg)[3*(i)]))

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
		      );

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
			  );

#endif /* UNSIGNEDHC_H */
