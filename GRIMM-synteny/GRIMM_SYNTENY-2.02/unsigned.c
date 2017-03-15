/* unsigned.c
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


/* Glenn Tesler */
/* Unsigned permutation distance algorithm

   Unichromosomal: proved in "To cut ... or not to cut", Hannenhalli & Pevzner

      Enhancement:
      If g1, g2 both start or end with a k-strip (k >= 1)
      then canonically sign it.  In particular, do not count it as
      a 1-strip or 2-strip with undetermined sign.

      This enhancement applies to unichrom linear/circular,
      but not to the multichromosomal alg below.

   Circular:
      If there is a common 3-strip, rotate/flip both genomes to start
      with it and then do normal unichromosomal algorithm.

      Else, must align both genomes to start (or end) w/some gene x,
      do the unichrom alg on it, then reverse one genome in positions 2..n
      (keeping x at start) and do unichrom alg on that too.

   Multichromosomal:

      Locate all singletons/2-strips, except
         caps break strips and do not count as singletons/2-strips, since
	 they will be reconnected;

	 1-gene chromos in either genome do not count as singletons,
	 because the capping alg gives their sign.

      Form G(Pi,Gamma), do unichrom alg on it but use components in IU
      instead of U, and use paths/cycles with 2-strips instead of cycles.
*/

/* TODO: try ALLCROSS/NOCROSS idea */

#include <stdio.h>
#include <stdlib.h>
#include "uniinvdist.h"
#include "mcstructs.h"
#include "scenario.h"
#include "mcrdist.h"
#include "unsigned.h"
#include "graph_edit.h"
#include "graph_components.h"
#include "e_malloc.h"
#include "circ_align.h"

void loop_unsigned_perm(spins_t *sp,
			int sing_index       /* which singleton we're at */
			);

void sign_2strips(spins_t *sp);



/* canonically sign a pair of genomes,
   and return pointers to singletons and 2-strips

     g1=gamma will stay as-is.  Used as reference for the "identity".
     g2=pi    will be changed to canonical signs on >=2-strips.

     locate/classify/point to:
     0. where in g1 each elt occurs (i.e., inverse perm, unsigned)
     1. the singletons (array of positions in g2)
     2. the 2-strips
        array of positions in g2; must distinguish left vs. right number.
*/

void assign_canon_signs(spins_t *sp, distmem_t *distmem, int verbose) {
  int i, i1, i2;
  int gene;

  int *iperm1 = distmem->perm1;

  int *genes1     = sp -> g1 -> genes;
  int *genes2     = sp -> loop_g2.genes;
  int *singletons = sp -> singletons;
  int *doubletons = sp -> doubletons;
  int num_genes   = sp -> num_genes;
  int num_sing=0, num_doub=0;    /* count # singletons & 2-strips */

  int strip_len;               /* running length of strip */
  int strip_dir, cur_dir;      /* direction
				  +1 = increasing,
				  -1 = decreasing,
				   0 = new strip, direction not yet known
				*/
  int strip_done;             /* boolean flag */
  int cur_gene, lowcap;
  int cur_pos_g1, last_pos_g1;

  int unichrom = sp->num_chromosomes == 0;  /* 0=multichrom, 1=unichrom */



  /* form the inverse perm of g1=gamma,
     ignoring signs on g1.  Want
     genes1[i] = j   =>   iperm1[|j|] = i
   */

#if 0
  /* 1-based index; could save one int, but then have to offset index
     every time */
  iperm1 = (int *) e_malloc((num_genes+1) * sizeof(int),
			    "iperm1 (assign_canon_signs)");
#endif

  for (i=0 ; i < num_genes ; i++) {
    gene = genes1[i];
    if (gene<0)
      gene = -gene;
    iperm1[gene] = i;
  }


  /* scan g2 left-to-right for strips.
     Terminate at caps or at start/end of perm.
     Canonically sign strips of length >= 2.
     Record positions of singletons and 2-strips:

     singletons = array { i1, i2, ..., ik }
            of increasing positions where singletons are located

     doubletons = array { j1, j2, ..., jm }
            of increasing positions where 2-strips are located.
	    If a 2-strip is located in position r,r+1, record r,
	    whether it's increasing or decreasing.

     */


  /* initialize loop_g2 with the genes of g2 */
  copy_genes(sp->g2->genes, genes2, num_genes);

  /* lowest cap #; caps break strips. */
  lowcap = num_genes - 2 * sp->num_chromosomes + 1;



#ifdef DEBUGu
  fprintf(outfile,"assign_canon_signs.  lowcap=%d\n",lowcap);
  fprintf(outfile,"input  g1: ");
  print_genes(genes1,num_genes);
  fprintf(outfile,"input  g2: ");
  print_genes(genes2,num_genes);
#endif


  /* unichrom mode:
   * if genomes start or end with identical strip of any length, incl.
   * 1-strip or 2-strip, then it should be identically signed.
   */
  i1 = 0;
  i2 = num_genes;
  if (unichrom) {
    /* canonically sign strip, if any, at start of perm */
    for (i=i1 ; i < i2; i++) {
      cur_gene = abs(genes2[i]);

      /* locate +/- same gene in g1 */
      cur_pos_g1 = iperm1[cur_gene];

      if (cur_pos_g1 == i) {
	/* canonically sign current element */
	genes2[i] = genes1[i];
      } else {
	break;
      }
    }
    i1 = i;

    /* canonically sign strip, if any, at end of perm */
    for (i=i2-1 ; i > i1; i--) {
      cur_gene = abs(genes2[i]);

      /* locate +/- same gene in g1 */
      cur_pos_g1 = iperm1[cur_gene];

      if (cur_pos_g1 == i) {
	/* canonically sign current element */
	genes2[i] = genes1[i];
      } else {
	break;
      }
    }
    i2 = i+1;
  }


  strip_len = 0;                         /* no strip yet */

  /* No init needed for these, but the compiler complains.
     The loop sets these up before they're used. */
  strip_dir = cur_dir = last_pos_g1 = 0;

  for (i=i1 ; i < i2 ; i++) {
    /* input signs in g2 are irrelevant */
    cur_gene = genes2[i];
    if (cur_gene<0)
      cur_gene = -cur_gene;

    /* locate +/- same gene in g1 */
    cur_pos_g1 = iperm1[cur_gene];
#ifdef DEBUGu
    /*    fprintf(outfile,"[%d]",cur_pos_g1); */
    fprintf(outfile,"[%d@(%d,%d)]",cur_gene,cur_pos_g1,i);
#endif

    /* are we continuing a strip? */
    strip_done = TRUE;                   /* assume it will be terminated */

    /* caps end strips */
    /* noncaps may or may not */
    if (cur_gene < lowcap) {
      cur_dir = cur_pos_g1 - last_pos_g1;

      /* determine if starting >=2-strip, and the direction of the strip */
      if (strip_len == 1) {
	if (cur_dir == 1 || cur_dir == -1) {
	  strip_dir = cur_dir;                    /* set direction */
#ifdef DEBUGu
	  fprintf(outfile,"<%d>",strip_dir);
#endif
	  genes2[i-1] = strip_dir * genes1[last_pos_g1];
	                                         /* canon sign 1st elt */
	}
      }

      /* determine if in a strip */
      if (strip_dir == cur_dir   &&  strip_len >= 1) {
	genes2[i] = strip_dir * genes1[cur_pos_g1];
	                                          /* canon sign cur elt */

	strip_done = FALSE;                       /* still in a strip */
      }
    }


    if (strip_done) {

      if (strip_len == 1) {                       /* record singletons */

	if (unichrom) {
	  singletons[num_sing++] = i-1;
#ifdef DEBUGu
	  fprintf(outfile,"s%d ",i-1);
#endif
	} else {
	  /* filter out caps and 1-gene chromos */

	  if (
	      abs(genes2[i-1] < lowcap)  /* filter out caps */

	      /* filter out 1-gene chromos on g2 */
	      && i >= 2
	      && (abs(genes2[i-2]) < lowcap || abs(genes2[i]) < lowcap)

	      /* filter out 1-gene chromos on g1 */
	      /* following not necessary because we filtered out caps */
	      /*   && last_pos_g1 >= 1 && last_pos_g1 < num_genes-1   */

	      && (abs(genes1[last_pos_g1-1]) < lowcap ||
		  abs(genes1[last_pos_g1+1]) < lowcap)) {

	    singletons[num_sing++] = i-1;
#ifdef DEBUGu
	    fprintf(outfile,"s%d ",i-1);
#endif
	  }
	}


      } else if (strip_len == 2) {                /* record 2-strips */
	doubletons[num_doub++] = i-2;

#ifdef DEBUGu
	fprintf(outfile,"d%d ",i-2);
#endif

      }                                           /* don't record others */
#ifdef DEBUGu
      else {
	fprintf(outfile,"o(%d,%d) ",i-strip_len,strip_len);
      }
#endif

      /* If current gene is a cap, it is not part of a strip.
         Else, current gene counts as 1st elt of a new strip. */
      strip_len = (cur_gene >= lowcap) ? 0 : 1;
    } else {
      /* continue strip */
      strip_len++;
#ifdef DEBUGu
      fprintf(outfile,".");
#endif
    }

    last_pos_g1 = cur_pos_g1;
  }

  /* close & record final strip, if necessary */
  if (strip_len == 1) {                       /* record singletons */
    /* singleton @ end only happens in unichrom case, so don't need to
       filter out 1-gene chromos */
    singletons[num_sing++] = i-1;
#ifdef DEBUGu
    fprintf(outfile,"s%d ",i-1);
#endif

  } else if (strip_len == 2) {                /* record 2-strips */
    doubletons[num_doub++] = i-2;

#ifdef DEBUGu
    fprintf(outfile,"d%d ",i-2);
#endif

  }                                           /* don't record others */
#ifdef DEBUGu
  else {
    fprintf(outfile,"o(%d,%d) ",i-strip_len,strip_len);
  }
  fprintf(outfile,"\n");
#endif

#if 0
  free(iperm1);
#endif


#ifdef DEBUGu
  fprintf(outfile,"output g1: ");
  print_genes(genes1,num_genes);
  fprintf(outfile,"output g2: ");
  print_genes(genes2,num_genes);
  fprintf(outfile,"%d singletons, %d 2-strips\n",
	  num_sing,num_doub);

  fprintf(outfile,"singletons at ");
  for (i=0 ; i<num_sing ; i++)
    fprintf(outfile,"%d ",singletons[i]);
  fprintf(outfile,"\n2-strips at ");
  for (i=0 ; i<num_doub ; i++)
    fprintf(outfile,"%d ",doubletons[i]);
  fprintf(outfile,"\n");
#endif

  if (verbose) {
    fprintf(outfile,"Number of Singletons:         \t%d\n",num_sing);
    fprintf(outfile,"Number of 2-strips:           \t%d\n",num_doub);
    fflush(outfile);
  }




  /* record outputs */
  sp -> num_sing = num_sing;
  sp -> num_doub = num_doub;
}



/* find the length of a cyclic strip
   starting @ positions pos1, pos2 (between 0 & num_genes-1)
   in arrays genes1, genes2 of length num_genes
   It goes in directions dir1, dir2 = +/-1
   The maximum length of the strip is max_len.
   Positions pos1, pos2 are assumed to agree already.
   Returns the length.
   */


int test_strip(int *genes1, int *genes2,
	       int pos1, int pos2,
	       int dir1, int dir2,
	       int num_genes,
	       int max_len) {
  int strip_len;
  int gene1, gene2;

#if 0
  fprintf(outfile,"[p1=%d,p2=%d, d1=%d,d2=%d, max_len=%d,", pos1,pos2,dir1,dir2, max_len);
#endif

  for (strip_len = 1; strip_len < max_len; strip_len++) {
    /* at start of iteration, strip_len positions have agreed */
    pos1 += dir1;
    if (pos1 == -1) pos1 = num_genes-1;
    else if (pos1 == num_genes) pos1 = 0;
    gene1 = genes1[pos1];

    pos2 += dir2;
    if (pos2 == -1) pos2 = num_genes-1;
    else if (pos2 == num_genes) pos2 = 0;
    gene2 = genes2[pos2];

    if (gene1 != gene2 && gene1 != -gene2) break;
  }

#if 0
  fprintf(outfile,"len=%d]\n", strip_len);
#endif


  return strip_len;
}


/* determine if an unsigned strip of length >= 3 exists between g1,g2
   in either direction.
   If not, return g1pos=-1.
   If so, g1pos, g2pos point to it, as follows:
      g2pos between 0 & num_genes-1: it starts in g2 @ pos g2pos & goes fwd
      g2pos = num_genes+x,  x between 0 & num_genes-1:
              it ends at x and goes backwards
      similar for g1pos      
*/
void find_3strip(struct genome_struct *g1,
		 struct genome_struct *g2,
		 int num_genes,
		 distmem_t *distmem,
		 int *g1pos,
		 int *g2pos) {

  int *iperm1 = distmem->perm1;

  int *genes1 = g1->genes;
  int *genes2 = g2->genes;

  int i;
  int gene;
  int cur_gene;
  int strip_len;


  /* default: no 3-strip */
  *g1pos = *g2pos = -1;

  /* no 3-strip if too few genes */
  if (num_genes<4) return;


  /* form the inverse perm of g1=gamma,
     ignoring signs on g1.  Want
     genes1[i] = j   =>   iperm1[|j|] = i
   */

  for (i=0 ; i < num_genes ; i++) {
    gene = genes1[i];
    if (gene<0)
      gene = -gene;
    iperm1[gene] = i;
  }

  for (i=0 ; i < num_genes ; i++) {
    /* i = possible value of g2pos = possible value at which to start strip */
    cur_gene = genes2[i];
    if (cur_gene<0)
      cur_gene = -cur_gene;


    /* test 3-strip fwds in g1 & g2 */
    strip_len = test_strip(genes1,genes2,
			   iperm1[cur_gene],i, /* start positions */
			   1,1,                /* directions */
			   num_genes,
			   3);
    if (strip_len == 3) {
      *g1pos = iperm1[cur_gene];
      *g2pos = i;
      break;
    }

    /* test 3-strip backwd in g1, fwd in g2 */
    strip_len = test_strip(genes1,genes2,
			   iperm1[cur_gene],i, /* start positions */
			   -1,1,               /* directions */
			   num_genes,
			   3);
    if (strip_len == 3) {
      *g1pos = iperm1[cur_gene] + num_genes;
      *g2pos = i;
      break;
    }
  }


  /* if 3-strip starts in g2 at position 0, it might be a longer strip
     starting earlier cyclically.
     */
  if (*g2pos == 0) {

    if (*g1pos < num_genes) {
      /* strip fwd in g1 & g2, see if extends other direction */
      strip_len = test_strip(genes1,genes2,
			     *g1pos,*g2pos,
			     -1,-1,
			     num_genes,
			     num_genes-2);  /* back up w/o rereading strip */

      /* the genomes are identical.  No need to extend the strip. */
      if (strip_len == num_genes-2)
	strip_len = 1;

      /* don't count the first character in the extension length */
      strip_len--;

      /* back up start of strip in g1 */
      *g1pos -= strip_len;
      if (*g1pos < 0) *g1pos += num_genes;

      /* back up start of strip in g2 */
      *g2pos -= strip_len;
      if (*g2pos < 0) *g2pos += num_genes;

    } else {
      /* strip backwd in g1, fwd in g2; see if extends other direction */
      strip_len = test_strip(genes1,genes2,
			     *g1pos,*g2pos,
			     1,-1,
			     num_genes,
			     num_genes-2);  /* back up w/o rereading strip */

      /* the genomes are identical.  No need to extend the strip. */
      if (strip_len == num_genes-2)
	strip_len = 1;

      /* don't count the first character in the extension length */
      strip_len--;

      /* advance start of strip in g1 */
      *g1pos += strip_len;
      if (*g1pos >= 2*num_genes) *g1pos -= num_genes;

      /* back up start of strip in g2 */
      *g2pos -= strip_len;
      if (*g2pos < 0) *g2pos += num_genes;
    }

  }

  /* need to adjust g1pos by 1 if backwards */
  if (*g1pos >= num_genes) {
    if (++*g1pos == 2*num_genes) *g1pos = num_genes;
  }

  return;


}






/* distmem should have been allocated already */
/* this allocates the fields of sp */
void alloc_strip_mem(struct genome_struct *g1,
		     struct genome_struct *g2,
		     int num_genes, int num_chromosomes,
		     distmem_t *distmem,
		     spins_t *sp)
{
  sp->num_genes = num_genes;
  sp->num_chromosomes = num_chromosomes;

  sp->best_d = 2*num_genes + 1;  /* bigger than any possible distance */

  sp->num_sing = sp->num_doub = 0; /* no singletons/2-strips yet */

  sp->singletons = (int *) e_malloc(num_genes * sizeof(int),
				    "sp->singletons");

  sp->doubletons = (int *) e_malloc(num_genes * sizeof(int),
				  "sp->doubletons");

  sp->g1 = g1;
  sp->g2 = g2;


  alloc_simple_genome(num_genes,&sp->loop_g2);
  alloc_simple_genome(num_genes,&sp->superspin_g2);
#if 0
  alloc_simple_genome(num_genes,&sp->best_g1);
#endif
  alloc_simple_genome(num_genes,&sp->best_g2);
  /* gene data for these will be initialized in assign_canon_signs */

  sp->distmem = distmem;

  
}

void free_strip_mem(spins_t *sp)
{
  free_simple_genome(&sp->best_g2);
#if 0
  free_simple_genome(&sp->best_g1);
#endif
  free_simple_genome(&sp->superspin_g2);
  free_simple_genome(&sp->loop_g2);


  free(sp->doubletons);
  free(sp->singletons);
}



/* g1 = gamma, g2 = pi */

/* works for uni and multichromosomal genomes;
   num_chromosomes=0 for uni, >0 for multi */

/* g1, g2 with the best signs & capping,
   are returned in distmem->cappedp1, distmem->cappedp2.
*/
int unsigned_mcdist(struct genome_struct *g1, struct genome_struct *g2,
		    int num_genes,
    		    int num_chromosomes,
		    int circular,
		    distmem_t *distmem,
		    int verbose,
		    int *dist_error) {

  spins_t sp;
  int dist = -1;   /* function always should set dist >= 0 */


  /* copy of genomes */
  struct genome_struct Dup_g1, Dup_g2;
  struct genome_struct *dup_g1,*dup_g2;

  int circular_tries = 1;     /* if circular, may need 2nd pass */
  int g1pos, g2pos;           /* where to circularly align g1, g2 */


#if 0
  fprintf(outfile,"unsigned_mcdist\n");
#endif
  if (circular) {
#if 0
    fprintf(outfile,"unsigned_mcdist circular\n");
#endif

    /* duplicate genomes */
    dup_g1 = &Dup_g1;
    alloc_simple_genome(num_genes, dup_g1);
    copy_genes(g1->genes, dup_g1->genes, num_genes);

    dup_g2 = &Dup_g2;
    alloc_simple_genome(num_genes, dup_g2);
    copy_genes(g2->genes, dup_g2->genes, num_genes);

    find_3strip(dup_g1, dup_g2, num_genes, distmem, &g1pos, &g2pos);
#if 0
    fprintf(outfile,"3-strip: g1 %d, g2 %d\n", g1pos, g2pos);
#endif
    if (g1pos >= 0) {
      circular_align_position(dup_g1->genes, g1pos, num_genes);
      circular_align_position(dup_g2->genes, g2pos, num_genes);
    }
    else {
      /* no 3 strips! */
      circular_align_2(dup_g1, dup_g2, num_genes,
		       0, TRUE);  /* align genomes on 1st gene of g1 */

      /* need to do unsigned dist and then
	 flip the 2nd genome and try again! */
      circular_tries = 2;
    }

  } else {
    /* linear or multichromosomal, don't need to duplicate genomes */
    dup_g1 = g1;
    dup_g2 = g2;
  }



  /* create own memory for the perms & trial perms, as in scenario generator */
  alloc_strip_mem(dup_g1, dup_g2,
		  num_genes, num_chromosomes,
		  distmem,
		  &sp);


  for ( /* circular_tries, already set */
       ;
       circular_tries > 0;
       circular_tries--) {


    /* assign canonical signs, locate singletons & 2-strips,
       and initialize the iterator over all ways of signing singletons */
    assign_canon_signs(&sp, distmem, verbose);



    /* CRUDE APPROXIMATION FOR NOW: only take the initial singletons,
       ignore the rest
       */
    *dist_error = FALSE;
    if (sp.num_sing > MAX_UNSIGNED_SING) {
      *dist_error = TRUE;
      sp.num_sing = MAX_UNSIGNED_SING;
    }



    /* loop over all signs on singletons */
    loop_unsigned_perm(&sp,
		       0); /* loop starting with 0th singleton */

    dist = sp.best_d;  /* get best distance */


    if (circular_tries == 2) {
      /* need to flip 2nd genome (after 1st gene) & try again */
      /*      reverse_in_place_genes(dup_g2->genes, 1, num_genes-1); */

      /* flip 2nd genome around 1st gene, keeping it in place */
      circular_align_position(dup_g2->genes, num_genes+1, num_genes);
    }
  }


  /* store best signs in distmem->cappedp1 & cappedp2 */
  if (num_chromosomes == 0) {
    copy_genes(dup_g1->genes, distmem->cappedp1, num_genes);
    copy_genes(sp.best_g2.genes, distmem->cappedp2, num_genes);
  } else {
    mcdist_capgraph(dup_g1, &sp.best_g2,
		    num_genes, num_chromosomes,
		    distmem,
		    FALSE);
  }


#ifdef DEBUGu
  fprintf(outfile,"Best signed perms:\n");
  fprintf(outfile,"g1: ");
  print_genes(distmem->cappedp1, num_genes);
  fprintf(outfile,"g2: ");
  print_genes(distmem->cappedp2, num_genes);
#endif

  /* free memory */

  if (circular) {
    free_simple_genome(dup_g2);
    free_simple_genome(dup_g1);
  }

  free_strip_mem(&sp);


  return dist;
}

/* return the optimal distance between g1 & g2 w/unsigned perms */
/* genes of g1, g2 overwritten with the
   optimal signs & cappings

   use unsigned_mcdist & manage it differently to avoid the overwriting
   */
int unsigned_mcdist_nomem(struct genome_struct *g1,
			  struct genome_struct *g2,
			  int num_genes, int num_chromosomes,
			  int circular,
			  int verbose,
			  int *dist_error) {
  distmem_t distmem;
  int dist;

  /* allocate memory */
  mcdist_allocmem(num_genes, num_chromosomes,
		  &distmem);


  /* compute distance */
  dist = unsigned_mcdist(g1, g2,
			 num_genes, num_chromosomes,
			 circular,
			 &distmem,
			 verbose,
			 dist_error);

  /* return the properly capped & signed genomes
     unichromosomal linear case: g1 same as input, g2 changes
         (except null bug may have swapped these at the beginning)
     unichromosomal circular, or multichromosomal: both change
   */
  if (num_chromosomes>0 || circular)
    copy_genes(distmem.cappedp1, g1->genes, num_genes);

  copy_genes(distmem.cappedp2, g2->genes, num_genes);



  /* free memory */
  mcdist_freemem(&distmem);

  return(dist);
}





/* loop over all signs on singletons */

void loop_unsigned_perm(spins_t *sp,
			int sing_index       /* which singleton we're at */
			)
{
  /* haven't yet set all signs, so recurse */
  if (sing_index < sp->num_sing) {
    loop_unsigned_perm(sp,
		       sing_index+1);

    /* toggle sign and try again */
    sp->loop_g2.genes[sp->singletons[sing_index]] =
      - sp->loop_g2.genes[sp->singletons[sing_index]];


    loop_unsigned_perm(sp,
		      sing_index+1);

#ifdef DEBUGu
    /* not necessary to toggle back.
       Useful for debugging to keep track of sign changes
       in a sequence */
    sp->loop_g2.genes[sp->singletons[sing_index]] =
      - sp->loop_g2.genes[sp->singletons[sing_index]];
#endif



  } else {
    sign_2strips(sp);
  }
}



/* now we have a specific spin on singletons
 * need to assign signs to the 2-strips & get the distance
 *
 * make the graph
 * loop over unoriented components or intrachrom. unor. comp.
 * change 1st 2-strip in it, if it has one, to anticanonical
 * compute new distance, update best signs & dist if applicable
 */



void sign_2strips(spins_t *sp)
{
#ifdef DEBUGu
  static int trialno = 0;
#endif


  int num_genes = sp -> num_genes;
  int num_chromosomes = sp -> num_chromosomes;
  struct genome_struct *g1 = sp -> g1;

  int num_2strips_toflip;       /* # 2-strips remaining to make anticanon */

  int v;                        /* position of 2-strip */
  int num_doub = sp -> num_doub;
  struct genome_struct *loop_g2 = &sp->loop_g2;
  struct genome_struct *ss_g  = &sp->superspin_g2;

  int cur_d;                    /* distance in current permutations */
#ifdef DEBUGu
  int canon_d;                  /* distance for canonical signing */
  int cur_d2;                   /* verification distance */
#endif

  int i;

  graph_t G;                    /* graph with canonical signs */
  pathcounts_t pcounts_G;       /* path counts in such graph */

  distmem_t *distmem
                  = sp -> distmem;
  component_t *components
                  = distmem -> components;
  int *cc         = distmem -> cc;
  int *cycle      = distmem -> labeled;


  int comp_no;
  int *doubletons = sp -> doubletons;




  /* Duplicate loop iterator g2.
     Changes to 2-strips are in the copy ss_g, not the iterator. */

  copy_genes(loop_g2->genes, ss_g->genes, num_genes);

  /* get distance and graph, prior to dealing w/antistrips */

  /* Build breakpoint graph to proper stage 
   * unichrom: full breakpoint graph
   * multichrom: G(Pi,Gamma) or Gbar(Pi,Gamma)
   */
  if (num_chromosomes == 0) {    /* unichromosomal distance */
    cur_d = invdist_noncircular_G(g1, ss_g,
				  0, /* offset */
				  num_genes,
				  distmem,
				  &G, &pcounts_G,
				  FALSE);

  } else {                           /* multichromosomal distance */
#if 0
    cur_d = mcdist_capgraph_connections(g1, ss_g,
					num_genes, num_chromosomes, distmem,
					&G, &pcounts_G,
					FALSE);
#endif
    mcdist_partial_G(g1, ss_g,
		     num_genes, num_chromosomes,
		     distmem,
		     &G, &pcounts_G);
  }
#ifdef DEBUGu
  /*  canon_d = cur_d;  */                 /* store original distance */
#endif


  /* Scan through all 2-strips.
   * First mark all (intrachrom) unoriented components PENDING,
   * other components DONE
   */
#ifdef DEBUGu
  fprintf(outfile,"(c#,or,ic,ss): ");
#endif

  num_2strips_toflip = 0;
  for (i=0 ; i < G.num_components ; i++) {
    if (!components[i].oriented && components[i].intrachrom) {
      /* intrachrom, unoriented component */
      num_2strips_toflip++;                   /* count # such components */
      components[i].scanstate = SCAN_PENDING; /* need to flip a strip */

    } else {
      components[i].scanstate = SCAN_DONE;
    }

#ifdef DEBUGu
    fprintf(outfile,"(%d,%d,%d,%d) ",
	    i, components[i].oriented, components[i].intrachrom,
	    components[i].scanstate);
#endif
  }

#ifdef DEBUGu
  fprintf(outfile,"\n");
#if 0
  fprintf(outfile,"cc:  ");
  for (i=0 ; i < 2*num_genes+2 ; i++)
    fprintf(outfile,"%4d ",cc[i]);
  fprintf(outfile,"\n");
#endif /* 0 */
#endif

  /* Then scan 2-strips, re-signing the first one encountered in each
     (intrachrom) unoriented component */
  if (num_2strips_toflip > 0)
  for (i=0 ; i < num_doub; i++) {
    v = doubletons[i];             /* position of 2-strip in g2 */
                                   /* in graph, it's at 2v+1,...,2v+4 */
    comp_no = cc[2*v+1];

    /* Caps don't count for strips for the unsigned
       perm algorithm, but counting the caps might have made it a
       3-strip, in which case the vertex @ 2r-1 could be part of
       an adjacency or tail etc. instead of an actual component. */
    if (comp_no>=0 && /* see above */
#if 0
	comp_no == cc[2*v+4] &&    /* test both ends of 2-strip in same comp */
#endif
	/* test both ends of 2-strip in same cycle */
	cycle[2*v+1] == cycle[2*v+4] &&
	components[comp_no].scanstate == SCAN_PENDING) {
      /* found a new unoriented component w/2-strip */

#ifdef DEBUGu
      fprintf(outfile,"Anticanonical strip @ %d,%d, component %d\n",
	      v,v+1,comp_no);
#endif

      /* mark the component */
      components[comp_no].scanstate = SCAN_DONE;


      /* turn it into an antistrip in pi */
      /* unichromosomal: no need to edit graph
	 multichrom: perhaps could speed up computaton of capping
	 by editing graph here & then continuing w/steps 18.9-19 */

      ss_g -> genes[v] = - ss_g -> genes[v] ;
      ss_g -> genes[v+1] = - ss_g -> genes[v+1] ;


      /* can stop early if have handled every nonoriented component */
      /* It's possible a nonoriented component doesn't have any 2-strips,
	 in which case this early abort will not occur */
      if (--num_2strips_toflip == 0) break;
    }

  }


  /* have set all signs.  Compute distance. */
  if (num_chromosomes == 0) {    /* unichromosomal distance */
    cur_d = invdist_noncircular(g1, ss_g,
				0,         /* offset */
				num_genes, distmem);
  } else {                           /* multichromosomal distance */
    cur_d = mcdist_noncircular(g1, ss_g,
			       num_genes, num_chromosomes, distmem,
			       FALSE);
  }

#ifdef DEBUGu
  fprintf(outfile,"Dist %d: ", cur_d);
  print_genes(ss_g->genes,num_genes);
#endif


  
  if (cur_d < sp->best_d) {
    /* have new best distance */

    sp->best_d = cur_d;

    copy_genes(ss_g->genes, sp->best_g2.genes, num_genes);
    /* get the capping later, if desired */

#if 0
    /* must do the capping now while we still have the graph */
    if (num_chromosomes > 0) {
      mcdist_capgraph_finalperms(g1, ss_g,
				 num_genes, num_chromosomes,
				 distmem,
				 &G,
				 (graphstats_t *)NULL);  /* not verbose */
      copy_genes(distmem->cappedp1, sp->best_g1.genes, num_genes);
      copy_genes(distmem->cappedp2, sp->best_g2.genes, num_genes);
    } else {
      copy_genes(ss_g->genes, sp->best_g2.genes, num_genes);
    }
#endif
  }
    
}







/***************************************************************************/

/* like setinvmatrix, but for multiple genomes */
void set_unsigneddist_matrix(int **distmatrix, struct genome_struct *genomes,
			     int num_genes, int num_chromosomes,
			     int num_genomes,
			     int circular,
			     distmem_t *distmem,
			     int *dist_error)
{
  int i, j;
  int new_dist_error;      /* set to TRUE if any distance is approximated */

  *dist_error = FALSE;

  for (i=0 ; i<num_genomes ; i++) {
    distmatrix[i][i] = 0;
    for (j=i+1 ; j<num_genomes ; j++) {
	distmatrix[i][j] = distmatrix[j][i] =
		unsigned_mcdist(genomes+i,genomes+j,
				num_genes,num_chromosomes,
				circular,
				distmem,
				FALSE,   /* not verbose */
				&new_dist_error);
	if (new_dist_error) *dist_error = TRUE;
    }
  }
  return;
}
