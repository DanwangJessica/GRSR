/* circ_align.c
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

#include "mcstructs.h"
#include "scenario.h"
#include "circ_align.h"

/* routines for realigning circular genomes */

/* invdist_circular does not require realigning genomes,
   but most everything else does */

/* circular_align: align list of genomes so they start or end w/same gene
   circular_align_2: same for 2 genomes

   align_on: 0 = use the gene that's at the start or end of 1st genome
            1,2,...,n or -1,...,-n: use this gene
   align_front: TRUE=align on first gene of genomes, FALSE=on last gene


   circular_align_1
       same for 1 genome, except align_on can't be 0
   
*/

void circular_align(struct genome_struct *genomes,
		    int num_genomes, int num_genes,
		    int align_on,
		    int align_front)
{
  int genome_num=0;

  if (align_on == 0) {
    /* get first or last gene of first genome */
    align_on = align_front
      ? genomes[0].genes[0]
      : genomes[0].genes[num_genes-1];

      /* and bump counter so we don't waste time aligning first genome */
      genome_num++;
  }

  /* align all genomes to start or end with gene align_on */
  for ( ; genome_num < num_genomes; genome_num++)
    circular_align_1(&genomes[genome_num], num_genes,
		     align_on, align_front);
}

void circular_align_2(struct genome_struct *g1, struct genome_struct *g2,
		      int num_genes,
		      int align_on,
		      int align_front)
{
  if (align_on == 0) {
    /* get first or last gene of first genome */
    align_on = align_front
      ? g1->genes[0]
      : g1->genes[num_genes-1];
  } else {
    /* align first genome */
    circular_align_1(g1, num_genes, align_on, align_front);
  }

  /* align 2nd genome */
  circular_align_1(g2, num_genes, align_on, align_front);
}

void circular_align_1(struct genome_struct *g1,
		      int num_genes,
		      int align_on,
		      int align_front)
{
  int g_pos;                     /* position of gene align_on in genome */
  int g_sign=1;                  /* +1 or -1 depending if it's +/- align_on */
  int *genes = g1->genes;


  /* locate gene align_on in the genome */
  for (g_pos=0; g_pos<num_genes; g_pos++) {
    if (genes[g_pos] == align_on) {
      g_sign = 1; break;
    }
    if (genes[g_pos] == -align_on) {
      g_sign = -1; break;
    }
  }



  /* how to reorder the data:
     a = gene to align on
     X, Y = other genes before and after a
     ' = reverse sequence but don't change signs
     - = reverse sequence, do change signs
     / = divides into blocks

     align front, sign +: X a Y = X / a Y -> X' / a Y -> X' / Y' a -> a Y / X
     align end,   sign +: X a Y = X a / Y -> a X' / Y -> a X' / Y' -> Y / X a
     (both are circular shifts)

     align front, sign -: X -a Y = X -a / Y -> a -X / Y -> a -X / -Y
     align end,   sign -: X -a Y = X / -a Y -> -X / -a Y -> -X / -Y a

  */


  if (g_sign == +1) {
    /* alignment involves 3 unsigned reversals */

    /* g_pos may be adjusted by 1 so that the two blocks are in positions
       0..(g_pos-1) and g_pos..(num_genes-1)
       */
    if (align_front) {
      if (g_pos == 0) return;            /* already aligned */
    } else {
      if (g_pos == num_genes-1) return;  /* already aligned */
      g_pos++;                           /* adjust start of block */
    }
    circular_align_position(genes, g_pos, num_genes);

  } else {
    /* alignment involves 2 signed reversals */

    /* g_pos may be adjusted by 1 so that the two blocks are in positions
       0..(g_pos-1) and g_pos..(num_genes-1)
       */
    if (align_front) {
      g_pos++;                           /* adjust start of block */
    }
    circular_align_position(genes, g_pos + num_genes, num_genes);
  }
}


void reverse_in_place_genes_unsigned(int *source_perm,
				     int start, int end)
{
  int k,l,temp;
  for (k=start,l=end; k<l; k++,l--) {
    temp = source_perm[k];
    source_perm[k] = source_perm[l];
    source_perm[l] = temp; 
  }
}


/* circular_align_position(genes, start, num_genes)
   if 0<=start<num_genes, circularly left shift genes array by start positions
   if start=num_genes+x (0<=x<num_genes),
       genes[0..x-1] genes[x..num_genes-1] ->
       -genes[x-1] .. -genes[0] -genes[num_genes-1] ... -genes[x]
*/
      
void circular_align_position(int *genes, int start, int num_genes)
{
  if (start < num_genes) {
    reverse_in_place_genes_unsigned(genes, 0, start-1);
    reverse_in_place_genes_unsigned(genes, start, num_genes-1);
    reverse_in_place_genes_unsigned(genes, 0, num_genes-1);
  } else {
    reverse_in_place_genes(genes, 0, start-num_genes-1);
    reverse_in_place_genes(genes, start-num_genes, num_genes-1);
  }
}
