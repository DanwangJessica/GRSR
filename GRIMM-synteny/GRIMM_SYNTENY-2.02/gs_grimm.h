/* gs_grimm.h
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

#ifndef GS_GRIMM_H
#define GS_GRIMM_H

void gs_grimm_alloc(int num_chromosomes,
		    int num_genes,
		    int num_genomes,

		    distmem_t **distmem_ret,
		    graphstats_t ***statsmatrix_ret,
		    struct genome_struct **genome_list_ret);

void gs_grimm_destroy(int num_chromosomes,
		      int num_genes,
		      int num_genomes,
		      distmem_t *distmem,
		      graphstats_t **statsmatrix,
		      struct genome_struct *genome_list);

#endif /* GS_GRIMM_H */
