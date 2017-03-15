/* texgraph.h
 *    Output a graph drawing for use with LaTeX package bpgraph.sty.
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

#ifndef TEXGRAPH_H
#define TEXGRAPH_H

void draw_texgraph(struct genome_struct *g1,
		   struct genome_struct *g2,
		   int num_genes, int num_chromosomes,
		   char *options);

void texgraph_nomem(struct genome_struct *g1,
		   struct genome_struct *g2,
		   int num_genes, int num_chromosomes,
		   int lr_or_cc,
		   int showlabels,
		   int stage);

void texgraph(struct genome_struct *g1,
	     struct genome_struct *g2,
	     int num_genes, int num_chromosomes,
	     int lr_or_cc,
	     int showlabels,
	     int stage,
	     distmem_t *distmem);



#define TEXGR_LR 0
#define TEXGR_CC 1
#define TEXGR_RAW 2     /* not TeX output, but easiest to put here:
			 * list the verts on each path/cycle
			 */

#define TEXGR_S 1       /* show signed (orig) entries of pi */
#define TEXGR_B 2       /* show chromosome brackets */
#define TEXGR_U 4       /* show unsigned (doubled) entries of pi */ 
#define TEXGR_AB 8      /* doubled entries of pi: show as #a, #b */

#endif /* TEXGRAPH_H */
