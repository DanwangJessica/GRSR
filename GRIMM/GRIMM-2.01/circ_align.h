/* circ_align.h
 *    Rotate all circular genomes to align them on a particular gene.
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

#ifndef CIRC_ALIGN_H
#define CIRC_ALIGN_H

#include "mcstructs.h"



void circular_align(struct genome_struct *genomes,
		    int num_genomes, int num_genes,
		    int align_on,
		    int align_front);

void circular_align_2(struct genome_struct *g1, struct genome_struct *g2,
		      int num_genes,
		      int align_on,
		      int align_front);

void circular_align_1(struct genome_struct *g1,
		      int num_genes,
		      int align_on,
		      int align_front);

void reverse_in_place_genes_unsigned(int *source_perm,
				     int start, int end);

void circular_align_position(int *genes, int start, int num_genes);

#endif /* CIRC_ALIGN_H */
