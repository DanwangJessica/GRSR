/* testrev.h
 *    Find operations that move 1 step in most parsimonious scenario.
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

#ifndef TESTREV_H
#define TESTREV_H

void print_rev_stats_t1(struct genome_struct *g1, struct genome_struct *g2,
			int num_genes, int num_chromosomes);

void print_rev_stats_t2(struct genome_struct *g1, struct genome_struct *g2,
			int num_genes, int num_chromosomes,
			int output_coded);

void print_rev_stats_t3(struct genome_struct *g1, struct genome_struct *g2,
			int num_genes, int num_chromosomes);

void init_cbounds_wmem(int num_genes, int num_chromosomes,
		       struct genome_struct *g1,
		       cbounds_t *cb);
void init_cbounds(int num_genes, int num_chromosomes,
		  struct genome_struct *g1,
		  cbounds_t *cb);

void free_cbounds(cbounds_t *cb);



void init_isBP_mem(int num_genes, int num_chromosomes,
		   char **breakpoints);
void clean_isBP_mem(char *breakpoints);
void init_isBP_array_2(int num_genes, int num_chromosomes,
		       struct genome_struct *g1, struct genome_struct *g2,
		       char *breakpoints);
void init_isBP_array(int num_genes, int num_chromosomes,
		     int num_genomes,
		     struct genome_struct *genome_list,
		     int gindex, /* which genome to compute b.p. against */
		     char *breakpoints);
int isBP(int g, char *breakpoints);




/* private declarations */

#define BP_FIS     " /"
#define BP_BP      " :"
#define BP_ABFUS     " | + |"
#define BP_ABNOFUS   " | . |"
#define BP_BAFUS     " | +"
/* #define BP_BANOFUS   " | ." */
#define BP_BANOFUS   ""


#define CREV 0         /* counts for reversals */
#define CTRA 1         /* counts for translocations */
#define CFIS 2         /* fissions */
#define CFUS 3         /* fusions */
#define CTOT0 4        /* total */
#define CNF0 5         /* number of fields */

#define CDOWN1 0       /* events lowering distance by 1 */
#define CSAME  1       /* events keeping distance the same */
#define CUP1   2       /* events raising distance by 1 */
#define CERROR 3       /* another amount, error */
#define CTOT1  4       /* total */
#define CFIRST 5       /* count # tries till first success */
#define CNF1   6       /* number of fields */

#endif /* TESTREV_H */
