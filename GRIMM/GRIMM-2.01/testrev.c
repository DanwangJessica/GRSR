/* testrev.c
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

#include <stdio.h>
#include <stdlib.h>
#include "uniinvdist.h"
#include "mcstructs.h"
#include "scenario.h"
#include "mcrdist.h"
#include "mcread_input.h"
#include "e_malloc.h"
#include "write_data.h"
#include "testrev.h"


typedef struct {
  struct genome_struct *trial_genome;
  struct genome_struct *dest_genome;
  int num_genes;
  int num_chromosomes;
  distmem_t *distmem;

  int counts[CNF0][CNF1];
} tryrevparams_t;



void calc_bp_array(int *genes, int num_genes,
		   int s, int e,
		   int bptype,
		   int *successors,
		   int *bplist,
		   int *bplen);


void print_bp_seg(int *genes,
		  int s,
		  int e);

void print_bp_delim(char *delim);

void init_try_operation(tryrevparams_t *p);

int try_operation(int optype,
		  int d,
		  tryrevparams_t *p);

void print_tryops_counts(tryrevparams_t *p);




/*****************************************************************************/
/*            Statistics on how many reversals decrease distance by 1        */
/*                                Method I                                   */
/* For multichromosomal genomes, this method is inaccurate: it gives the     */
/* ops which decrease distance by 1 FOR THIS CAPPING/CONCATENATION.  But     */
/* other capping/concatenations may allow other operations to decrease dist  */
/* by 1.                                                                     */
/*****************************************************************************/

void print_rev_stats_t1(struct genome_struct *g1, struct genome_struct *g2,
			int num_genes, int num_chromosomes) {
  distmem_t distmem;
  struct genome_struct current_genome, trial_genome;
  struct genome_struct *dest_genome = g2;

  int num_up1,num_down1,num_same,num_error;
  int num_bkpt;
  int num_tillfound;        /* how many reversals tested till 1st success? */
  int foundany;             /* has a reversal been found at all? */
  int foundany_i;           /* has a reversal been found for this start? */
  int i,j;
  int d, d2;
  int lowcap;

  int *successors;  /* array of successors in final genome */


  init_scenario_mem(num_genes,num_chromosomes,&distmem,
		    &current_genome,&trial_genome,
		    g1,
		    0);
  init_scenario_bp(num_genes,num_chromosomes,
		   /*&distmem,*/
		    g2,
		    0,
		    &successors);


  /* count # breakpoint positions */
  num_bkpt = 0;
  for (i=0; i<num_genes; i++) {
      if (i==0) {
	if (num_chromosomes > 0) continue;
      } else {
	if (SUCCESSOR(current_genome.genes[i-1]) == current_genome.genes[i])
	  continue;
      }
      num_bkpt++;
  }
  /* does END count as a breakpoint? Yes for unichrom, no for multichrom */
  if (num_chromosomes == 0) num_bkpt++;

  fprintf(outfile,"Breakpoint positions (%d):",num_bkpt);
  for (i=0; i<num_genes; i++) {
      if (i==0) {
	if (num_chromosomes > 0) continue;
      } else {
	if (SUCCESSOR(current_genome.genes[i-1]) == current_genome.genes[i])
	  continue;
      }
      fprintf(outfile," %d",i);
  }
  /* does END count as a breakpoint? Yes for unichrom, no for multichrom */
  if (num_chromosomes == 0) {
      fprintf(outfile," %d",num_genes);
  }
  fprintf(outfile,"\n\n");



  /* Now try ALL reversals starting/ending at breakpoints.
     Count how many change the distance by -1, 0, +1, other (error).
     Output coded list of all the ones that decrease the distance by 1.
     */

  num_down1 = num_up1 = num_same = num_error = 0;
  num_tillfound = 0;
  foundany = FALSE;

  /* UNIchromosomal reversal distance computation */
  d = invdist_noncircular(&current_genome, dest_genome,
			  0, /* offset */
			  num_genes, &distmem);
  /* the lowest cap number */
  lowcap = num_genes - 2*num_chromosomes + 1; 


  for (i=0 ; i<num_genes ; i++) {
    foundany_i = FALSE;

    /* don't split signed 2-strip at start of reversal */
    if (i==0) {
      if (num_chromosomes > 0) continue;
    } else {
      if (SUCCESSOR(current_genome.genes[i-1]) == current_genome.genes[i])
	continue;
    }

  
    for (j=i ; j<num_genes ; j++) {

      /* don't split signed 2-strip at end of reversal */
      if (j==num_genes-1) {
	if (num_chromosomes > 0) continue;
      } else {
	if (SUCCESSOR(current_genome.genes[j]) == current_genome.genes[j+1])
	  continue;
      }


      copy_and_reverse_genes(current_genome.genes,trial_genome.genes,i,j, num_genes);
	
      d2 = invdist_noncircular(&trial_genome, dest_genome,
			       0, /* offset */
			       num_genes, &distmem);

      if (!foundany) num_tillfound++;

      if (d2 == d-1) {
	num_down1++;

	if (!foundany_i) {
	  fprintf(outfile,"%d: ",i);
	  foundany = foundany_i = TRUE;
	}

	fprintf(outfile," %d",j);

      } else if (d2 == d) {
	num_same++;
      } else if (d2 == d+1) {
	num_up1++;
      } else {
	num_error++;
	fprintf(outfile,"\nERROR: Potential distance caculation error;\n");
	fprintf(outfile," Reversal in positions %d..%d changed distance from %d to %d\n",i,j,d,d2);
      }
    }
    if (foundany_i) fprintf(outfile,"\n");
  }




  fprintf(outfile,"\nNumber of reversals changing distance by\n   -1 is %d;  0 is %d;  +1 is %d;  other is %d;  total: %d\n",
	  num_down1,num_same,num_up1,num_error,
	  num_down1+num_same+num_up1+num_error);
  fprintf(outfile,"First success after %d trials\n",num_tillfound);



  clean_scenario_bp_mem(successors);
  clean_scenario_mem(&current_genome,&trial_genome,&distmem);
}

/*****************************************************************************/
/*            Method II (& variation IV)                                     */
/*            Do all possible multichromosomal ops, and                      */
/*            print a report of how many change dist by +1, 0, -1,           */
/*            as well as report all the ops changing it by -1.               */
/*****************************************************************************/

/*
   (t2) output_coded = FALSE: print out each genome in full
   (t4) output_coded = TRUE: print out codes saying which operations were
                             performed and in which positions
*/


/*****************************/
/* Subroutines for Method II */
/*****************************/

/* init chromosome boundary arrays */
/* init_cbounds: allocates memory for ptrs in cb
 * init_cbounds_wmem: the pointers in cb already point to
 * blocks of memory large enough
 */

void init_cbounds_wmem(int num_genes, int num_chromosomes,
		       struct genome_struct *g1,
		       cbounds_t *cb) {

  int i;
  int inchrom, cnum;
  int lowcap;
  int *genes, g;

  lowcap = num_genes - 2*num_chromosomes + 1;

  cnum = 0;
  inchrom = FALSE;
  genes = g1->genes;

  for (i=0 ; i<num_genes ; i++) {
    g = genes[i];
    if (g<0)   g = -g;

    if (g >= lowcap) {
      if (inchrom) {
	/* end of chromosome */
	inchrom = FALSE;
      } else {
	/* start of chromosome */
	inchrom = TRUE;
	cb->cBound[cnum] = i;
	cnum++;
      }
    }
    cb->cNum[i] = cnum;
  }

  if (num_chromosomes > 0)
    cb->cBound[cnum] = num_genes;  /* end of last chromosome */

#if 0
  if (num_chromosomes == 0) {
    cb->cBound[0] = 0;
    cb->cBound[1] = num_genes;
  } else {
    cb->cBound[cnum] = num_genes;
  }
#endif

#if 0
  fprintf(outfile,"cb->cNum: ");
  for (i=0 ; i < num_genes ; i++) {
    fprintf(outfile,"%d ",cb->cNum[i]);
  }
  fprintf(outfile,"\n");

  fprintf(outfile,"cb->cBound: ");
  for (i=0 ; i <= num_chromosomes ; i++) {
    fprintf(outfile,"%d ",cb->cBound[i]);
  }
  fprintf(outfile,"\n");
#endif
  
}

void init_cbounds(int num_genes, int num_chromosomes,
		  struct genome_struct *g1,
		  cbounds_t *cb) {
  cb->cNum     = (int *) e_malloc(num_genes*sizeof(int), "cb->cNum");
  cb->cBound   = (int *) e_malloc((num_chromosomes+1)*sizeof(int),
				  "cb->cBound");

  init_cbounds_wmem(num_genes, num_chromosomes, g1, cb);
}



void free_cbounds(cbounds_t *cb) {
  free(cb->cBound);
  free(cb->cNum);
}




/*********************************************************/
/* create/maintain isBP "array":                        */
/* If a genome has consecutive genes (a,b) or (-b,-a),  */
/* then the position between the two genes is called:   */
/* "after a", "before b", "after -b", "before -a"       */
/* and isBP(-a)=isBP(b) is boolean valued:              */
/*             =1  if the position is a breakpoint      */
/*             =0  if it is not a breakpoint            */
/*                                                      */
/* Special gene values:                                 */
/*     isBP(num_genes+1) = 1      if end is b.p.        */
/*     isBP(-(num_genes+1)) = 1   if start is b.p.      */
/*                                                      */
/* maintain it as a bit array                           */
/* for genes g=x,   x=1,...,num_genes, use bit 2x-2     */
/* for genes g=-x,  x=1,...,num_genes, use bit 2x-1     */
/********************************************************/

/* allocate memory for breakpoint array.
   Return pointer  (char *) *breakpoints
   */

#define BPsize(num_genes)   (num_genes/4 + 1)

void init_isBP_mem(int num_genes, int num_chromosomes,
		   char **breakpoints)
{
  /* # bits required = 2*(num_genes+1), rounded up to a byte */
  /* # bytes = floor((2*(num_genes+1) + 7)/8) = floor(num_genes/4) + 1 */

  *breakpoints = (char *) e_malloc( BPsize(num_genes) * sizeof(char),
				   "breakpoints");
  
}

void clean_isBP_mem(char *breakpoints)
{
  free(breakpoints);
}



/* set or clear whole b.p. array */
/* value = 0 or 1 */
void set_isBP_array_whole(int num_genes, char *breakpoints, int value)
{
  int bp_bytes = BPsize(num_genes);
  int i;
  char setto = value ? 0xff : 0;

  for (i=0; i<bp_bytes; i++)
    breakpoints[i] = setto;
}


/* set_isBP_array_1: set bit for gene g in b.p. array
   set_isBP_array_0: clear bit for gene g
   */
void set_isBP_array_1(int g, char *breakpoints)
{
  /* for genes g=x,   x=1,...,num_genes, use bit 2x-2     */
  /* for genes g=-x,  x=1,...,num_genes, use bit 2x-1     */

  int bitno = (g>0)  ?   (g-1)<<1  :  (((-g) - 1)<<1)  + 1 ;


  /* break it up into byte number and bit w/in byte */
  int byteno = bitno >> 3;
  bitno = bitno & 7;


  breakpoints[byteno] |=   1<<bitno;
}


void set_isBP_array_0(int g, char *breakpoints)
{
  /* for genes g=x,   x=1,...,num_genes, use bit 2x-2     */
  /* for genes g=-x,  x=1,...,num_genes, use bit 2x-1     */

  int bitno = (g>0)  ?   (g-1)<<1  :  (((-g) - 1)<<1)  + 1 ;


  /* break it up into byte number and bit w/in byte */
  int byteno = bitno >> 3;
  bitno = bitno & 7;


  breakpoints[byteno] &=   0xff ^ (1<<bitno);
}


/* Methods for determining what are breakpoints:
   Method 1. All points are breakpoints.
   Method 2. All caps are 2-sided b.p.
             Non-caps: if adjacency for genomes  g[i0],g[j]  for all j,
                       it's not a b.p., else it is a b.p.
   Method 3. The fancy b.p. graph / hurdle analysis
*/
void init_isBP_array_m1(int num_genes, int num_chromosomes,
			int num_genomes,
			struct genome_struct *genome_list,
			int gindex,
			char *breakpoints)
{
  /* set all genes to b.p. */
  set_isBP_array_whole(num_genes, breakpoints, 1);
}




/* breakpoint determination method 2:
   Caps are 2-sided b.p.
   For adjacencies not involving caps, if (a,b)/(-b,-a) is in all genomes,
   there's no b.p. before b or -a, else there is.

   In unichromo case, genome starts (ends) w/b.p. iff any genomes differ in
   starting (ending) gene
   */
void init_isBP_array_m2(int num_genes, int num_chromosomes,
			int num_genomes,
			struct genome_struct *genome_list,
			int gindex,
			char *breakpoints)
{
  int lowcap = num_genes - 2*num_chromosomes + 1;

  int gindex2;         /* secondary genome */
  int i;               /* loop over gene positions */
  int gene1, gene2;

  int *successors;     /* SUCCESSOR array for primary genome */


  /* determine successor array of primary genome */
  /* The memory allocated for this could be done just once
     if this routine is repeatedly called */
  init_scenario_bp(num_genes, num_chromosomes,
		   &genome_list[gindex],
		   0,             /* no pad, it's already accounted for */
		   &successors);

  /* clear all b.p.'s */
  set_isBP_array_whole(num_genes, breakpoints, 0);



  /* if multichrom, the start and end will always be considered b.p.
     These values are never examined, but to be consistent with
     how the positions between internal caps are treated, call them b.p.
     */

  if (num_chromosomes>0) {
      /* b.p. after start position */
      set_isBP_array_1(-(num_genes+1), breakpoints);

      /* b.p. before end position */
      set_isBP_array_1(num_genes+1, breakpoints);
  }




  /* Form union of b.p. over all pairs (primary genome, other genome) */
  for (gindex2 = 0; gindex2 < num_genomes; gindex2++) {
    /* skip if this is the primary genome */
    if (gindex2 == gindex)
      continue;


    /* add b.p. for this genome to the breakpoint list */
    /* first check internal b.p., and then check start/end b.p. */
    for (i=0; i<num_genes-1; i++) {
      gene1 = genome_list[gindex2].genes[i];
      gene2 = genome_list[gindex2].genes[i+1];

#if 0
      fprintf(outfile,
	      "g1=%d, g2=%d, SUCC(g1)=%d, SUCC(-g2)=%d",
	      gene1, gene2, SUCCESSOR(-gene1), SUCCESSOR(-gene2));
#endif

      /* caps are always b.p.
	 adjacencies not involving caps are ignored */

      if (SUCCESSOR(gene1) == gene2              /* pair is adjacency */
	  &&
	  gene1 < lowcap && gene1 > -lowcap      /* gene1 not a cap */
	  &&
	  gene2 < lowcap && gene2 > -lowcap      /* gene2 not a cap */
	  ) {
#if 0
	fprintf(outfile,"\n");
#endif
	continue;
      }
#if 0
      fprintf(outfile,"   B.P.\n");
#endif


      /* add b.p. to list */
      set_isBP_array_1(gene2, breakpoints);  /* b.p. before gene2 */
      set_isBP_array_1(-gene1, breakpoints); /* b.p. after gene1 */
    }


    /* now check start/end b.p. */

    if (num_chromosomes > 0) {
      /* multichrom: these values should never be examined.
	 To be consistent w/the positions between internal cocaps,
	 call it a breakpoint.
	 */


      /* b.p. before leading cap */
      gene2 = genome_list[gindex2].genes[0];
      set_isBP_array_1(gene2, breakpoints);

      /* b.p. after start position: already set */
      /* set_isBP_array_1(-(num_genes+1), breakpoints); */


      /* b.p. after ending cap */
      gene2 = genome_list[gindex2].genes[num_genes-1];
      set_isBP_array_1(-gene2, breakpoints);

      /* b.p. before end position: already set */
      /* set_isBP_array_1(num_genes+1, breakpoints); */

    } else {
      /* unichrom: b.p. at start (end) if genomes start (end) w/different genes
       */

      /* check start */
      gene1 = genome_list[gindex].genes[0];
      gene2 = genome_list[gindex2].genes[0];

      if (gene1 != gene2) {
	/* b.p. before first gene */
	set_isBP_array_1(gene2, breakpoints);

	/* b.p. after start position */
	set_isBP_array_1(-(num_genes+1), breakpoints);
      }



      /* check end */
      gene1 = genome_list[gindex].genes[num_genes-1];
      gene2 = genome_list[gindex2].genes[num_genes-1];

      if (gene1 != gene2) {
	/* b.p. after last gene */
	set_isBP_array_1(-gene2, breakpoints);

	/* b.p. before end position */
	set_isBP_array_1(num_genes+1, breakpoints);
      }
    }


  } /* for (gindex2...) */



  clean_scenario_bp_mem(successors);
}


/* THIS ROUTINE IS INCOMPLETE.  DO NOT ACTIVATE IT AT THIS TIME. */
/* breakpoint determination method 3:

   All positions are b.p. _unless_ they form an "unbreakable adjacency"
   for some secondary genome w.r.t. the primary genome.


   In unichromo case, genome starts (ends) w/b.p. iff any genomes differ in
   starting (ending) gene
   */
void init_isBP_array_m3(int num_genes, int num_chromosomes,
			int num_genomes,
			struct genome_struct *genome_list,
			int gindex,
			char *breakpoints)
{
  int lowcap = num_genes - 2*num_chromosomes + 1;

  int gindex2;         /* secondary genome */
  int i;               /* loop over gene positions */
  int gene1, gene2;

  int *successors;     /* SUCCESSOR array for primary genome */


  /* determine successor array of primary genome */
  /* The memory allocated for this could be done just once
     if this routine is repeatedly called */
  init_scenario_bp(num_genes, num_chromosomes,
		   &genome_list[gindex],
		   0,             /* no pad, it's already accounted for */
		   &successors);


  /* set all positions to b.p. */
  set_isBP_array_whole(num_genes, breakpoints, 1);


#if 0
  /* if multichrom, the start and end will always be considered b.p.
     These values are never examined, but to be consistent with
     how the positions between internal caps are treated, call them b.p.
     */

  if (num_chromosomes>0) {
      /* b.p. after start position */
      set_isBP_array_1(-(num_genes+1), breakpoints);

      /* b.p. before end position */
      set_isBP_array_1(num_genes+1, breakpoints);
  }
#endif




  /* Clear b.p. if they form an "unbreakable adjacency" */
  for (gindex2 = 0; gindex2 < num_genomes; gindex2++) {
    /* skip if this is the primary genome */
    if (gindex2 == gindex)
      continue;



    /* form the breakpoint graph */
    if (num_chromosomes == 0) {
    } else {
      /* multichromosomal */
      /* TODO */
    }


    /* add b.p. for this genome to the breakpoint list */
    /* first check internal b.p., and then check start/end b.p. */
    for (i=0; i<num_genes-1; i++) {
      gene1 = genome_list[gindex2].genes[i];
      gene2 = genome_list[gindex2].genes[i+1];

      /* caps are always b.p.
	 adjacencies not involving caps are ignored */

      if (SUCCESSOR(gene1) != gene2               /* pair is not adjacency */
	  &&
	  gene1 < lowcap && gene1 > -lowcap       /* gene1 not a cap */
	  &&
	  gene2 < lowcap && gene2 > -lowcap)      /* gene2 not a cap */
	continue;

      /* TODO:
	 check the breakpoint graph to see if this is unbreakable
	 */
      if (1)
	continue;

      /* if we made it this far, it's unbreakable, so remove b.p. from list */

      set_isBP_array_0(gene2, breakpoints);  /* b.p. before gene2 */
      set_isBP_array_0(-gene1, breakpoints); /* b.p. after gene1 */
    }


    /* now check start/end b.p. */

    if (num_chromosomes > 0) {
      /* multichrom: these values should never be examined.
	 To be consistent w/the positions between internal cocaps,
	 call it a breakpoint.
	 */

#if 0
      /* b.p. before leading cap */
      gene2 = genome_list[gindex2].genes[0];
      set_isBP_array_1(gene2, breakpoints);

      /* b.p. after start position: already set */
      /* set_isBP_array_1(-(num_genes+1), breakpoints); */


      /* b.p. after ending cap */
      gene2 = genome_list[gindex2].genes[num_genes-1];
      set_isBP_array_1(-gene2, breakpoints);

      /* b.p. before end position: already set */
      /* set_isBP_array_1(num_genes+1, breakpoints); */
#endif

    } else {
      /* unichrom: b.p. at start (end) if genomes start (end) w/different genes
       */

      /* check start */
      gene1 = genome_list[gindex].genes[0];
      gene2 = genome_list[gindex2].genes[0];
      if (gene1 == gene2) {	  /* no b.p. at start of genome */
	if (0) {                  /* TODO: test if its unbreakable */
	  /* if made it this far, it's unbreakable adjacency, remove b.p. */
	  /* b.p. before first gene */
	  set_isBP_array_0(gene2, breakpoints);

	  /* b.p. after start position */
	  set_isBP_array_0(-(num_genes+1), breakpoints);
	}
      }



      /* check end */
      gene1 = genome_list[gindex].genes[num_genes-1];
      gene2 = genome_list[gindex2].genes[num_genes-1];
      if (gene1 == gene2) {	  /* no b.p. at end of genome */
	if (0) {                  /* TODO: test if its unbreakable */
	  /* b.p. after last gene */
	  set_isBP_array_0(-gene2, breakpoints);

	  /* b.p. before end position */
	  set_isBP_array_0(num_genes+1, breakpoints);
	}
      }
    }


  } /* for (gindex2...) */



  clean_scenario_bp_mem(successors);
}





/* initialize the breakpoint array for a list of genomes */
void init_isBP_array(int num_genes, int num_chromosomes,
		     int num_genomes,
		     struct genome_struct *genome_list,
		     int gindex, /* which genome to compute b.p. against */
		     char *breakpoints)
{
#if 0
  int i;
#endif

  /* pick method _m1, _m2, or _m3 */
  init_isBP_array_m2(num_genes,num_chromosomes,
		     num_genomes, genome_list,
		     gindex,
		     breakpoints);

#if 0
  fprintf(outfile,"Breakpoint array: ");
  for (i=0; i<BPsize(num_genes); i++)
    fprintf(outfile," %02x", (unsigned char) breakpoints[i]);
  fprintf(outfile,"\n");
#endif
}


/* initialize the breakpoint array for 2 genomes */
void init_isBP_array_2(int num_genes, int num_chromosomes,
		       struct genome_struct *g1, struct genome_struct *g2,
		       char *breakpoints)
{
  struct genome_struct    genome_list[2];

  genome_list[0].genes = g1->genes;
  genome_list[1].genes = g2->genes;

  init_isBP_array(num_genes, num_chromosomes,
		  2, /* 2 genomes */
		  genome_list,
		  0, /* compute b.p. relative to genome 0 */
		  breakpoints);
}







/* check if a position is a breakpoint */
int isBP(int g, char *breakpoints)
{
  /* for genes g=x,   x=1,...,num_genes, use bit 2x-2     */
  /* for genes g=-x,  x=1,...,num_genes, use bit 2x-1     */

  int bitno = (g>0)  ?   (g-1)<<1  :  (((-g) - 1)<<1)  + 1 ;


  /* break it up into byte number and bit w/in byte */
  int byteno = bitno >> 3;
  bitno = bitno & 7;

  /* retrieve bit */
  return (breakpoints[byteno] >> bitno) & 1;
}


/******************************/
/* Full routine for Method II */
/******************************/


void print_rev_stats_t2(struct genome_struct *g1, struct genome_struct *g2,
			int num_genes, int num_chromosomes,
			int output_coded) {
  distmem_t distmem;
  struct genome_struct current_genome, trial_genome;
  /*  struct genome_struct *dest_genome = g2; */
  struct genome_struct Dest_genome, *dest_genome=&Dest_genome;


  int optype;

  int num_bkpt;

  int foundany_i;      /* has a reversal been found for this start? */
  int i,j;
  int d;
  int lowcap;

  int c1,c2;           /* chromosomes being examined */
  int sc1, sc2;        /* signed versions, in case need to flip them */
  int f1;              /* is chromo c1 flipped? */

  int s1=0,e1=0,s2=0,e2=0;  /* start & end of chromosomes being examined */
  int s1i, e1i;        /* start/end of interior, dropping caps if multichrom */
  int s2_orig=0, e2_orig=0;  /* original start/end of c2, in case it's moved */
  int move_c2;         /* boolean, indicating if c2 is moved to abut c1 */


  /* FUNCTIONALITY MOVED to init_isBP_array */
  /* test at breakpoints only, or test at all locations? */
  /*  int testbp = TRUE;  */
  /*  int testbp = FALSE; */

  cbounds_t cb;        /* structure with chromosome boundaries */

  tryrevparams_t tryrevparams;


  char *breakpoints;    /* bit array of genes that are breakpoints */

  /* we pad the genomes with a null chromosome to allow fissions */
  int num_genes_p = num_genes + 1;                /* # genes, incl. pad */
  int num_chromosomes_p = num_chromosomes + 1;    /* # chromos, incl. pad */
  int num_chromosomes_nn;                      /* # nonnull chromos in g1 */

  cb.cNum = (int *) NULL;
  cb.cBound = (int *) NULL;


  /* determine # genes in padded genome */
  if (num_chromosomes == 0) {
    num_genes_p = num_genes;
    num_chromosomes_p = num_chromosomes;
  } else {
    num_genes_p = num_genes + 2;
    num_chromosomes_p = num_chromosomes + 1;
  }


  init_scenario_mem(num_genes,num_chromosomes,&distmem,
		    &current_genome,&trial_genome,
		    g1,

		    /* if multichromosomal, add a null chromosome
		       so we can do fissions */
		    num_chromosomes > 0);

#if 0
  /* This was separated out from init_scenario_mem,
     and then was moved down below, AFTER dest_genome created w/pads */
  init_scenario_bp(num_genes,num_chromosomes,
		   /*&distmem,*/
		    g1,g2,

		    /* if multichromosomal, add a null chromosome
		       so we can do fissions */
		    num_chromosomes > 0,

		    &successors);
#endif

  /* copy g2 to destination genome, and pad with null chromosome at end
   to match the null added onto copy of g1 in previous call */
  alloc_simple_genome(num_genes_p, dest_genome);
  copy_genes(g2->genes,
	     dest_genome->genes,
	     num_genes);

  /* Need nulls at end regardless of how capping was done,
   * reformat them now */
  reformat_caps(&current_genome, num_genes, num_chromosomes);
  reformat_caps(dest_genome, num_genes, num_chromosomes);


  if (num_chromosomes > 0) {
    /* add null chromo for doing fissions */
    dest_genome->genes[num_genes_p-2] = num_genes_p - 1;
    dest_genome->genes[num_genes_p-1] = num_genes_p;

    /* get the chromosome boundaries */
    init_cbounds(num_genes_p,
		 num_chromosomes_p,
		 &current_genome,
		 &cb);

    /* get the number of nonnull chromosomes */
    /* find last nonnull chromosome in g1 */
    num_chromosomes_nn = num_chromosomes;
    i = num_genes;
    while (num_chromosomes_nn > 0 && current_genome.genes[i-2] == i-1) {
      /* ends in two caps, which is our standardized null */
      i -= 2;
      num_chromosomes_nn--;
    }
  } else {
    /* unichromosomal case: number of nonnull chromosomes considered 1 */
    num_chromosomes_nn = 1;
  }    



  /* initialize breakpoint array */
  /* functionality separated out and changed from init_scenario_mem above */
  init_isBP_mem(num_genes_p,num_chromosomes_p,
		&breakpoints);
  init_isBP_array_2(num_genes_p,num_chromosomes_p,
		    &current_genome,dest_genome,
		    breakpoints);




  /************************************************************/
  /* count # breakpoint positions                             */
  /* complexities of caps not fully appreciated in this tally */
  /************************************************************/

  if (output_coded) {
    fprintf(outfile,"Breakpoint positions:");

    /* count # breakpoint positions */
    num_bkpt = 0;
    for (i=0; i<num_genes; i++) {
      if (!isBP(current_genome.genes[i], breakpoints))
	continue;
      fprintf(outfile," %d",i);
      num_bkpt++;
    }
    /* does END count as a breakpoint? */
    if (isBP(num_genes+1, breakpoints)) {
      fprintf(outfile," %d",num_genes+1);
      num_bkpt++;
    }

    fprintf(outfile," (%d total)\n", num_bkpt);
  }


  /************************************************************************/
  /*           init remaining params to test all operations               */
  /************************************************************************/


  /* initial distance between the genomes */
  d = mcdist_noncircular(&current_genome, dest_genome,
			 num_genes, num_chromosomes,
			 &distmem,
			 FALSE);

  /* the lowest cap number */
  lowcap = num_genes - 2*num_chromosomes + 1; 

  /* initialize parameters used for every distance trial */
  init_try_operation(&tryrevparams);   /* zero counters in all categories */
  tryrevparams.trial_genome = &trial_genome;
  tryrevparams.dest_genome = dest_genome;
  tryrevparams.num_genes = num_genes_p;
  tryrevparams.num_chromosomes = num_chromosomes_p;
  tryrevparams.distmem = &distmem;


  /* categorize by operation type */

  /************************************************************************/
  /*               Test all 1 chromo operations                           */
  /************************************************************************/

  /* start and end of the first null chromosome;
     used in fission operations */
  if (num_chromosomes>0) {
    s2 = cb.cBound[num_chromosomes_nn+1 -1];
    e2 = cb.cBound[num_chromosomes_nn+1] -1;
  }

  for (c1 = 1; c1 <= num_chromosomes_nn; c1++) {
    if (output_coded) fprintf(outfile,"\nC%d:\n", c1);

    /* start & end of chromosome */
    if (num_chromosomes > 0) {
      s1 = cb.cBound[c1-1];
      e1 = cb.cBound[c1]-1;
    } else {
      s1 = 0;
      e1 = num_genes-1;
    }

    /**************/
    /*  fissions  */
    /**************/
  if (num_chromosomes>0) {


    /* for fissions, need to use padded # genes, chromos */
    tryrevparams.num_genes = num_genes_p;
    tryrevparams.num_chromosomes = num_chromosomes_p;

    foundany_i = FALSE;

    for (i=s1+2; i<e1 ; i++) {
      /* test the fission starting just before position i */
      if (!isBP(current_genome.genes[i], breakpoints))
	continue;

      copy_and_reverse_genes(current_genome.genes,
			     trial_genome.genes,
			     i,s2, num_genes_p);

      if (try_operation(CFIS,d,&tryrevparams)) {
	if (output_coded) {
	  if (!foundany_i) {
	    foundany_i = TRUE;
	    fprintf(outfile,"F: ");
	  }
	  fprintf(outfile,"%d ",i);    /* fissions in chromosome */
	} else
	  print_genome_multichrom(&trial_genome,
				  num_genes_p, num_chromosomes_p);
      }
    }
    if (output_coded && foundany_i) fprintf(outfile,"\n");
  }
    
    /************************/
    /*  internal reversals  */
    /************************/

    /* for reversals, need to use unpadded # genes, chromos */
    tryrevparams.num_genes = num_genes;
    tryrevparams.num_chromosomes = num_chromosomes;

    /* determine start/end of interior of chromosome, by
       removing caps if multichromosomal */
    if (num_chromosomes==0) {
      s1i = s1; e1i = e1;
    } else {
      s1i = s1+1; e1i = e1-1;
    }

    for (i=s1i ; i <= e1i ; i++) {
      foundany_i = FALSE;

      /* only start at breakpoints */
      if (!isBP(current_genome.genes[i], breakpoints))
	continue;


      for (j=i ; j<=e1i ; j++) {

	if (i == s1i && j == e1i && num_chromosomes>0)
	  continue; /* no chromosome flips in multichromosomal mode */

	/* only end at breakpoints */
	if (!isBP(-current_genome.genes[j], breakpoints))
	  continue;

	copy_and_reverse_genes(current_genome.genes,
			       trial_genome.genes,
			       i,j, num_genes);

	if (try_operation(CREV,d,&tryrevparams)) {
	  if (output_coded) {
	    if (!foundany_i) {
	      fprintf(outfile,"R%d: ",i); /* start reversals in chromosome */
	      foundany_i = TRUE;
	    }
	    fprintf(outfile," %d",j);
	  } else
	    print_genome_multichrom(&trial_genome, num_genes, num_chromosomes);
	}
      } /* end for j (ending point within chromosome) */

      if (output_coded && foundany_i) fprintf(outfile,"\n");
    } /* end for i (starting point within chromosome) */
  } /* end for c1 */



  /************************************************************************/
  /*               Test all 2 chromo operations                           */
  /************************************************************************/

  if (num_chromosomes>0) {

    /* for fusions and translocations, need to use unpadded # genes, chromos */
    tryrevparams.num_genes = num_genes;
    tryrevparams.num_chromosomes = num_chromosomes;

    for (c1=1 ; c1 < num_chromosomes_nn ; c1++) {
      for (f1=0 ; f1<=1 ; f1++) { /* do we flip chromo c1 or not? */
	/* start & end of chromosome */
	s1 = cb.cBound[c1-1];
	e1 = cb.cBound[c1]-1;


	if (f1) {
	  /* when we flip this chromosome, we will never be looking
	     at any earlier chromosomes again, so there is no
	     need to restore it later */

	  reverse_in_place_genes(current_genome.genes, s1,e1);
	}


	for (c2=c1+1 ; c2 <= num_chromosomes_nn ; c2++) {

	  /* identify chromosomes involved */
	  /* - indicates chromo is flipped */
	  sc1 = f1 ? -c1 : c1;
	  sc2 = c2;

	  /* start & end of chromosome, pointing at caps */
	  s2 = cb.cBound[c2-1];
	  e2 = cb.cBound[c2]-1;

	  /* speedup: move chromo c2 to just after c1 */
	  move_c2 = (e1+1 < s2);
	  if (move_c2) {
	    s2_orig = s2;
	    e2_orig = e2;

	    /* from:  xxx 1 A 2 yyy 3 B 4 zzz */
	    /* to:    xxx 1 A 2 -4 -B -3 -yyy zzz */
	    reverse_in_place_genes(current_genome.genes,
				   e1+1,e2);

	    /* new location */
	    s2 = e1+1;
	    e2 = s2+(e2_orig-s2_orig);

	    /* indicate chromo2 was flipped */
	    sc2 = -sc2;
	  }
	  


	  if (output_coded) {
	    fprintf(outfile,"\nP%d,%d:\n",
		    sc1, sc2);
	  }



	  /*******************************/
	  /*  fusions & translocations   */
	  /*******************************/

	  for (i=s1+1 ; i <= e1 ; i++) {
	    foundany_i = FALSE;

	    /* only start at breakpoints */
	    if (!isBP(current_genome.genes[i], breakpoints))
	      continue;

	    for (j=s2 ; j<e2 ; j++) {

	      /* must take a nonnull portion of at least one chromosome;
		 otherwise it's just a cap flip */
	      if (i==e1 && j==s2)
		continue;

	      /* also disallow exchanging both chromos */
	      if (i==s1+1 && j==e2-1)
		continue;

	      /* only end at breakpoints */
	      if (!isBP(-current_genome.genes[j], breakpoints))
		continue;


	      copy_and_reverse_genes(current_genome.genes,
				     trial_genome.genes,
				     i,j, 
				     num_genes);

	      optype = CTRA;

	      if ((i == s1+1 && j == s2) ||    /* fusion -chrom1 chrom2 */
		  (i == e1   && j == e2-1)) {  /* fusion chrom1 -chrom2 */
		optype = CFUS;
	      }
		
	      if (try_operation(optype, d, &tryrevparams)) {
		if (output_coded) {
		  if (!foundany_i) {
		    fprintf(outfile,"T%d: ",i); /* start reversals in chromosome */
		    foundany_i = TRUE;
		  }
		  fprintf(outfile," %d",j);
		} else
		  print_genome_multichrom(&trial_genome, num_genes, num_chromosomes);
	      }
	    } /* end for j (ending point within chromosome c2) */

	    if (output_coded && foundany_i) fprintf(outfile,"\n");
	  } /* end for i (starting point within chromosome c1) */


	  /* if moved c2 for speed, move it back */
	  if (move_c2) {
	    /* from:  xxx 1 A 2 -4 -B -3 -yyy zzz */
	    /* to:    xxx 1 A 2 yyy 3 B 4 zzz */
	    reverse_in_place_genes(current_genome.genes,
				   e1+1,e2_orig);

	  }

	} /* end: c2 */
      } /* end: f1 */
    } /* end: c1 */
  } /* end: if (num_chromosomes>0) */


  /********************************************************/

  /* print summary */
  print_tryops_counts(&tryrevparams);


  /* clean up memory */

  clean_isBP_mem(breakpoints);

  if (num_chromosomes>0)
    free_cbounds(&cb);

  free_simple_genome(dest_genome);

  clean_scenario_mem(&current_genome,&trial_genome,&distmem);
}

  


/****************************************************************************/
/*                               Method III                                 */
/* Similar to Method II, but produces operations in an order suitable for   */
/* the graphical web display.                                               */
/****************************************************************************/

void print_rev_stats_t3(struct genome_struct *g1, struct genome_struct *g2,
			int num_genes, int num_chromosomes) {
  distmem_t distmem;
  struct genome_struct current_genome, trial_genome;
  /*  struct genome_struct *dest_genome = g2; */
  struct genome_struct Dest_genome, *dest_genome=&Dest_genome;


  int i,j;
  int d;

  int c1,c2;           /* chromosomes being examined */
  int sc1, sc2;        /* signed versions, in case need to flip them */
  int f1;              /* is chromo c1 flipped? */
  int doubflip;        /* do we flip 2 chromos to do a fusion? */

  int s1=0,e1=0,s2=0,e2=0;   /* start & end of chromosomes being examined */
  int s2_orig=0, e2_orig=0;  /* original start/end of c2, in case it's moved */
  int move_c2;         /* boolean, indicating if c2 is moved to abut c1 */


  /* breakpoint lists */
  /* bptype TRUE: ordinary breakpoints; FALSE: every elt is breakpoint */
  int bptype = TRUE;
  /*  int bptype = FALSE; */

  int *bpmem, *bp_list1, *bp_list2;  /* list of b.p. in 1 or 2 chromos */
  int bp_len1, bp_len2;              /* lengths of those lists */
  int bp_len1f, bp_len2f;
  int bp_max1;

  int *potential_revs;   /* iterators for potential reversals
			    potential_revs[i]=j
			    means the next reversal attempted from
			    b.p. block i is through b.p. block j,
			    where i <= j < bplen
			    j=bplen when no possible reversals left
			    */

  int bp_start, bp_end;  /* b.p. block #s for start/end of reversal; */
  int min_bp_start;      /* 1st block to examine on new line */
  int last_bp_end;       /* last block # in an "on" region on this line,
			    -1 if no prior block on this line */
  int gap;

  int fusionAB, fusionBA;

  tryrevparams_t tryrevparams;

  cbounds_t cb;        /* structure with chromosome boundaries */

  int *successors;     /* array of successors in final genome */




  /* we pad the genomes with a null chromosome to allow fissions */
  int num_genes_p = num_genes + 1;                /* # genes, incl. pad */
  int num_chromosomes_p = num_chromosomes + 1;    /* # chromos, incl. pad */
  int num_chromosomes_nn;                      /* # nonnull chromos in g1 */

  cb.cNum = (int *) NULL;
  cb.cBound = (int *) NULL;

  /* determine space for padded genome */
  if (num_chromosomes == 0) {
    num_genes_p = num_genes;
    num_chromosomes_p = num_chromosomes;
  } else {
    num_genes_p = num_genes + 2;
    num_chromosomes_p = num_chromosomes + 1;
  }


  init_scenario_mem(num_genes,num_chromosomes,&distmem,
		    &current_genome,&trial_genome,
		    g1,

		    /* if multichrom, add a null chromosome
		       so we can do fissions */
		    num_chromosomes > 0);

  init_scenario_bp(num_genes,num_chromosomes,
		   /*&distmem,*/
		   g2,

		   /* if multichrom, add a null chromosome
		      so we can do fissions */
		   num_chromosomes > 0,

		   &successors);


  /* copy g2 to destination genome, and pad with null chromosome at end
   to match the null added onto copy of g1 in previous call */
  alloc_simple_genome(num_genes_p, dest_genome);
  copy_genes(g2->genes,dest_genome->genes,num_genes);

  /* Need nulls at end regardless of how capping was done,
   * reformat them now */
  reformat_caps(&current_genome, num_genes, num_chromosomes);
  reformat_caps(dest_genome, num_genes, num_chromosomes);


  if (num_chromosomes > 0) {
    /* add null chromo for doing fissions */
    dest_genome->genes[num_genes_p-2] = num_genes_p - 1;
    dest_genome->genes[num_genes_p-1] = num_genes_p;

    /* get the chromosome boundaries */
    init_cbounds(num_genes_p,
		 num_chromosomes_p,
		 &current_genome,
		 &cb);

    /* get the number of nonnull chromosomes */
    /* find last nonnull chromosome in g1 */
    num_chromosomes_nn = num_chromosomes;
    i = num_genes;
    while (num_chromosomes_nn > 0 && current_genome.genes[i-2] == i-1) {
      /* ends in two caps, which is our standardized null */
      i -= 2;
      num_chromosomes_nn--;
    }
  } else {
    /* unichromosomal case: number of nonnull chromosomes considered 1 */
    num_chromosomes_nn = 1;
  }




  /* initial distance between the genomes */
  d = mcdist_noncircular(&current_genome, dest_genome,
			 num_genes, num_chromosomes,
			 &distmem,
			 FALSE);




  /* memory for breakpoints */
  /* most memory is for unichromosomal mode.
     In multichromo mode, some of this is unused;
     the most needed would be (a+b+2)*sizeof(int),
     where a,b are the sizes of the 2 largest chromos
     */
  bpmem = (int *) e_malloc((num_genes+1)*sizeof(int), "bpmem");
  bp_list1 = bpmem;

  potential_revs = (int *) e_malloc(num_genes*sizeof(int), "bpmem");


  /* initialize parameters used for every distance trial */
  init_try_operation(&tryrevparams);   /* zero counters in all categories */
  tryrevparams.trial_genome = &trial_genome;
  tryrevparams.dest_genome = dest_genome;
  tryrevparams.num_genes = num_genes_p;
  tryrevparams.num_chromosomes = num_chromosomes_p;
  tryrevparams.distmem = &distmem;



  /* Now try ALL rearrangements starting/ending at breakpoints.
     Count how many change the distance by -1, 0, +1, other (error).
     Output coded list of all the ones that decrease the distance by 1.
     */

  /* categorize by operation type */

  /************************************************************************/
  /*               Test all 1 chromo operations                           */
  /************************************************************************/

  /* start and end of the first null chromosome;
     used in fission operations */
  if (num_chromosomes>0) {
    s2 = cb.cBound[num_chromosomes_nn+1 -1];
    e2 = cb.cBound[num_chromosomes_nn+1] -1;
  }

  for (c1 = 1; c1 <= num_chromosomes_nn; c1++) {
    fprintf(outfile,"C%d:\n", c1);

    /* start & end of chromosome */
    if (num_chromosomes > 0) {
      /* start, end exlude caps */
      s1 = cb.cBound[c1-1]+1;   /* point to first gene */
      e1 = cb.cBound[c1]-1;     /* point to 1 past last gene (i.e., to rcap) */
    } else {
      s1 = 0;
      e1 = num_genes_p;
    }


    /* initialize breakpoint list */
    calc_bp_array(current_genome.genes, num_genes,
		  s1,e1,  /* exclude caps */
		  bptype,
		  successors,
		  bp_list1,
		  &bp_len1);


    /****************************************/
    /*  print chromosome, showing fissions  */
    /****************************************/

    if (num_chromosomes == 0) {
      for (i=0 ; i<=bp_len1-2 ; i++) {
	if (i>0)
	  print_bp_delim(BP_BP);
	print_bp_seg(current_genome.genes,
		     bp_list1[i], bp_list1[i+1]);
      }
      fprintf(outfile,"\n");
    } else {
      /* print first segment of genes */
      print_bp_seg(current_genome.genes,
		   bp_list1[0], bp_list1[1]);

      /* for fissions, need to use padded # genes, chromos */
      tryrevparams.num_genes = num_genes_p;
      tryrevparams.num_chromosomes = num_chromosomes_p;

      for (i=1; i<=bp_len1-2 ; i++) {
	/* test the fission starting just before breakpoint  i */

	copy_and_reverse_genes(current_genome.genes,
			       trial_genome.genes,
			       bp_list1[i],s2,
			       num_genes_p);
	if (try_operation(CFIS,d,&tryrevparams))
	  print_bp_delim(BP_FIS);
	else
	  print_bp_delim(BP_BP);

	print_bp_seg(current_genome.genes,
		     bp_list1[i], bp_list1[i+1]);
      }
      fprintf(outfile,"\n");
    }


    /************************/
    /*  internal reversals  */
    /************************/

    /* for reversals, need to use unpadded # genes, chromos */
    tryrevparams.num_genes = num_genes;
    tryrevparams.num_chromosomes = num_chromosomes;




    /* if unichromosomal, the start & end may or may not be breakpoints,
       so delete them if appropriate: */
    min_bp_start = 0;
    if (bptype && num_chromosomes==0) {
      if (current_genome.genes[0] == dest_genome->genes[0]) {
	/* chromos start w/same gene, so it's not potential start of rev */
	/* bp_list1++; bp_len1--; */
	min_bp_start = 1;
      }

      if (current_genome.genes[num_genes-1]
	  == dest_genome->genes[num_genes-1]) {
	/* chromos end w/same gene, so it's not potential end of rev */
	bp_len1--;
      }
      
    }


    bp_max1 = bp_len1 - 1;

    /* initialize reversal extent iterators
       potential_revs[i]=j
       means the next reversal attempted from
       b.p. block i is through b.p. block j,
       where i < j < bplen
       j=bplen when no possible reversals left
     */
    for (i=0 ; i<bp_max1 ; i++)
      potential_revs[i] = i;


    last_bp_end = -1;
    bp_start = min_bp_start;

    /* loop till all reversals examined */
    while (1) {

      /* find next potential starting position */
      while (bp_start < bp_max1 && potential_revs[bp_start] >= bp_max1)
	bp_start++;

      if (bp_start >= bp_max1) {
	/* terminate previous line if necessary */
	if (last_bp_end > -1)
	  fprintf(outfile,"\n");

	/* find first possible starting interval */
	while (min_bp_start < bp_max1 &&
	       potential_revs[min_bp_start] >= bp_max1)
	  min_bp_start++;

	if (min_bp_start >= bp_max1)
	  break; /* all reversals have been tested */
	
	/* first place to look for a reversal start */
	bp_start = min_bp_start;
	last_bp_end = -1; /* last bp region on this line ended @ line start */
      }

      /* examine reversals starting here */
      while (potential_revs[bp_start] < bp_max1) {
	bp_end = potential_revs[bp_start]++;

	/* chromo flips not allowed */
	if (bp_start == 0 && bp_end == bp_max1-1
	    && num_chromosomes>0) break;

	copy_and_reverse_genes(current_genome.genes,
			       trial_genome.genes,
			       bp_list1[bp_start],
			       bp_list1[bp_end+1]-1,
			       num_genes);

	if (try_operation(CREV,d,&tryrevparams)) {
	  /* print # b.p. blocks to turn off */
	  gap = bp_start - last_bp_end - 1;
#if 0
	  fprintf(outfile,"bp_start=%d bp_end=%d last_bp_end=%d gap=%d\n",bp_start,bp_end,last_bp_end,gap);
#endif
	  if (gap > 0)
	    fprintf(outfile," -%d",gap);

	  /* print # b.p. blocks to turn on */
	  fprintf(outfile," %d", bp_end - bp_start + 1);

	  /* where to start next interval */
	  last_bp_end = bp_end;
	  bp_start = bp_end + 1;
	  if (bp_start >= bp_max1) break; /* this line is over */
	}

      }
    }
    fprintf(stdout,"\n"); /* blank line after all ops on this chromo */
  } /* end for c1 */


  /************************************************************************/
  /*               Test all 2 chromo operations                           */
  /************************************************************************/


  if (num_chromosomes>0) {

    /* for fusions and translocations, need to use unpadded # genes, chromos */
    tryrevparams.num_genes = num_genes;
    tryrevparams.num_chromosomes = num_chromosomes;

    for (c1=1 ; c1 < num_chromosomes_nn ; c1++) {
      for (f1=0 ; f1<=1 ; f1++) { /* do we flip chromo f1 or not? */
	/* start & end of chromosome */
	s1 = cb.cBound[c1-1]+1;         /*  ->  1st gene */
	e1 = cb.cBound[c1]-1;           /*  ->  1 past last gene (rcap) */


	if (f1) {
	  /* when we flip this chromosome, we will never be looking
	     at any earlier chromosomes again, so there is no
	     need to restore it later */

	  /* could either reverse just the genes (s1, e1-1)
	     or could also include the caps (s1-1, e1)
	     */
	  reverse_in_place_genes(current_genome.genes, s1, e1-1);
	}

	/* get b.p. blocks for this chromo */
	calc_bp_array(current_genome.genes, num_genes,
		      s1,e1,  /* exclude caps */
		      bptype,
		      successors,
		      bp_list1,
		      &bp_len1);
	bp_list2 = bp_list1 + bp_len1;  /* b.p. for chromo 2 listed afterwards */



	for (c2=c1+1 ; c2 <= num_chromosomes_nn ; c2++) {

	  /* identify chromosomes involved */
	  /* - indicates chromo is flipped */
	  sc1 = f1 ? -c1 : c1;
	  sc2 = c2;


	  /* start & end of chromosome */
	  s2 = cb.cBound[c2-1]+1;    /* -> first gene */
	  e2 = cb.cBound[c2]-1;      /* -> 1 past last gene (i.e., rcap) */



	  /* speedup: move chromo c2 to just after c1 */
	  move_c2 = (e1+1 < s2-1);
	  if (move_c2) {
	    s2_orig = s2;
	    e2_orig = e2;

	    /* from:  xxx 1 A 2 yyy 3 B 4 zzz */
	    /* to:    xxx 1 A 2 -4 -B -3 -yyy zzz */
	    reverse_in_place_genes(current_genome.genes,
				   e1+1,e2);

	    /* new location */
	    s2 = e1+2;
	    e2 = s2+(e2_orig-s2_orig);

	    /* indicate chromo2 was flipped */
	    sc2 = -sc2;
	  }




	  /* get b.p. blocks for this chromo */
	  calc_bp_array(current_genome.genes, num_genes,
			s2,e2,  /* exclude caps */
			bptype,
			successors,
			bp_list2,
			&bp_len2);


	  /*******************************/
	  /*         fusions             */
	  /*******************************/

	  /* test fusion AB */

	  /* have  xxx 1 A 2 yyy 3 B 4 zzz
	     xxx,yyy,zzz = other whole chromos
	     1,2,3,4 = caps
	     A = chromo c1
	     B = chromo c2
	     */
	  copy_genes(current_genome.genes,
		     trial_genome.genes,
		     num_genes);

	  /* xxx 1 A -B -3 -yyy -2 4 zzz */
	  reverse_in_place_genes(trial_genome.genes, e1, e2-1);

	  /* xxx 1 A B -3 -yyy -2 4 zzz */
	  reverse_in_place_genes(trial_genome.genes,
				 e1, e1+(e2-s2)-1);


	  fusionAB = try_operation(CFUS, d, &tryrevparams);





	  /* test fusion BA = (-A)(-B) */
	  /* xxx 1 A 2 yyy 3 B 4 zzz */
	  copy_genes(current_genome.genes,
		     trial_genome.genes,
		     num_genes);

	  /* xxx 1 A -B -3 -yyy -2 4 zzz */
	  reverse_in_place_genes(trial_genome.genes, e1, e2-1);

	  /* xxx 1 -A -B -3 -yyy -2 4 zzz */
	  reverse_in_place_genes(trial_genome.genes,
				 s1, e1-1);


	  fusionBA = try_operation(CFUS, d, &tryrevparams);


	  /* +=fusion, |=no fusion */
	  /* A|B|  A+B|  A+B+  stay as is
	     A|B+ should be transformed to (-A)+(-B)|
                  because looks better aesthetically */
	  doubflip = !fusionAB  && fusionBA;
	  if (doubflip) {
	    sc1 = -sc1; sc2 = -sc2; /* flipping both chromosomes */

	    fusionAB = TRUE;
	    fusionBA = FALSE;

	    /* flip both chromos in-place in current_genome */
	    /* xxx 1 A 2 yyy 3 B 4 zzz */
	    reverse_in_place_genes(current_genome.genes, s1, e1-1); /* -A */
	    reverse_in_place_genes(current_genome.genes, s2, e2-1); /* -B */

	    /* recalculate b.p. blocks for chromo c1 */
	    calc_bp_array(current_genome.genes, num_genes,
			  s1,e1,  /* exclude caps */
			  bptype,
			  successors,
			  bp_list1,
			  &bp_len1f);


	    /* DEBUGGING.  SHOULD NEVER OCCUR. */
	    if (bp_len1 != bp_len1f) {
	      fprintf(stdout,"ERROR: bp_len1=%d, bp_len1f=%d\n",
		      bp_len1, bp_len1f);
	      bp_len1 = bp_len1f;
	      bp_list2 = bp_list1 + bp_len1;
	    }


	    /* recalculate b.p. blocks for chromo c2 */
	    calc_bp_array(current_genome.genes, num_genes,
			  s2,e2,  /* exclude caps */
			  bptype,
			  successors,
			  bp_list2,
			  &bp_len2f);


	    /* DEBUGGING.  SHOULD NEVER OCCUR. */
	    if (bp_len2 != bp_len2f) {
	      fprintf(stdout,"ERROR: bp_len2=%d, bp_len2f=%d\n",
		      bp_len2, bp_len2f);
	      bp_len2 = bp_len2f;
	    }


	  }




	  /*******************************/
	  /*   print chrom #s            */
	  /*******************************/

	  fprintf(outfile,"\nP%d,%d:\n",
		  sc1, sc2);





	  /*******************************/
	  /*   print genes in chromos    */
	  /*******************************/

	  /* chromo 1: */
	  for (i = 0 ; i <= bp_len1-2 ; i++) {
	    if (i>0)
	      print_bp_delim(BP_BP);
	    print_bp_seg(current_genome.genes,
			 bp_list1[i], bp_list1[i+1]);
	  
	  }

	  /* fusion AB or not? */
	  print_bp_delim(fusionAB ? BP_ABFUS : BP_ABNOFUS);

	  /* chromo 2: */
	  for (i = 0 ; i <= bp_len2-2 ; i++) {
	    if (i>0)
	      print_bp_delim(BP_BP);
	    print_bp_seg(current_genome.genes,
			 bp_list2[i], bp_list2[i+1]);
	  
	  }

	  /* fusion BA or not? */
	  print_bp_delim(fusionBA ? BP_BAFUS : BP_BANOFUS);

	  fprintf(outfile,"\n");




	  /*******************************/
	  /*       translocations        */
	  /*******************************/

	  for (i=0 ; i <= bp_len1-1 ; i++) {   /* start at start of block i */
	    for (j=0 ; j <= bp_len2-1 ; j++) { /* end before start of block j */

	      /* events to discard:
		 cap flip:        i=bp_len1-1, j=0
		 fusion (-A)+B:   i=0, j=0
		 fusion A+(-B):   i=bp_len1-1, j=bp_len2-1
		 chromo exchange: i=0, j=bp_len2-1
		 */
	      if ((i==0 || i == bp_len1-1) &&
		  (j==0 || j == bp_len2-1))
		continue;


	      copy_and_reverse_genes(current_genome.genes,
				     trial_genome.genes,
				     bp_list1[i],bp_list2[j]-1,
				     num_genes);

	      if (try_operation(CTRA, d, &tryrevparams)) {
		if (i>0)
		  fprintf(outfile," -%d",i);  /* # blocks off */

		/* # blocks on =
		   # blocks remaining in chrom 1
		   + # initial blocks of chrom 2
		   + 1 for crossing chrom boundary
		   */
		fprintf(outfile," %d\n",
			bp_len1-i-1 + j + 1);
	      }

	    } /* end for j (ending block within chromosome c2) */
	  } /* end for i (starting block within chromosome c1) */



	  /* if did     A B  ->  -A -B,  need to put it back */
	  if (doubflip) {
	    /* flip both chromos in-place in current_genome */
	    /* xxx 1 -A 2 yyy 3 -B 4 zzz */
	    reverse_in_place_genes(current_genome.genes, s1, e1-1); /* A */
	    reverse_in_place_genes(current_genome.genes, s2, e2-1); /* B */

	    /* recalculate b.p. blocks for chromo c1 */
	    calc_bp_array(current_genome.genes, num_genes,
			  s1,e1,  /* exclude caps */
			  bptype,
			  successors,
			  bp_list1,
			  &bp_len1);
	  }

	  /* if moved c2 for speed, move it back */
	  if (move_c2) {
	    /* from:  xxx 1 A 2 -4 -B -3 -yyy zzz */
	    /* to:    xxx 1 A 2 yyy 3 B 4 zzz */
	    reverse_in_place_genes(current_genome.genes,
				   e1+1,e2_orig);

	  }

	} /* end: c2 */
      } /* end: f1 */
    } /* end: c1 */
  } /* end: if (num_chromosomes>0) */

  /********************************************************/

  /* print summary */
  print_tryops_counts(&tryrevparams);


  /* clean up memory */

  if (num_chromosomes>0)
    free_cbounds(&cb);

  free_simple_genome(dest_genome);

  clean_scenario_bp_mem(successors);
  clean_scenario_mem(&current_genome,&trial_genome,&distmem);
}


/* calculate indices of breakpoints
   in genes[s..e-1]
   these are the genes, no caps, so caps are at s-1, e if it's capped

   position i is a breakpoint if not a successor of position i-1
   start s, end e  are considered breakpoints

   bptype = TRUE: ordinary breakpoints
   bptype = FALSE: break between every entry

   output: bplist = array of breakpoint positions, increasing order
           bplen = # breakpoints, incl. start and end
*/

void calc_bp_array(int *genes,
		   int num_genes,
		   int s, int e,
		   int bptype,
		   int *successors,
		   int *bplist,
		   int *bplen)
{
  int i;

  /* break between all entries */
  if (!bptype) {
    for (i=s; i<=e; i++)
      *bplist++ = i;
    *bplen = e-s+1;
    return;
  }


#if 0
  fprintf(outfile,"breakpoints: %d", s);
#endif

  /* break only at breakpoints */
  *bplist++ = s; *bplen=1;

  for(i=s+1; i<e; i++) {
    if (SUCCESSOR(genes[i-1]) != genes[i]) {
      *bplist++ = i; ++*bplen;
#if 0
      fprintf(outfile, " %d", i);
#endif
    }
  }

  *bplist++ = e; ++*bplen;

#if 0
  fprintf(outfile, " %d\n", e);
#endif

}
		   

void print_bp_seg(int *genes,
		  int s,
		  int e)
{
  int i;

  for(i = s ; i<e ; i++)
    fprintf(outfile," %d", genes[i]);
}

void print_bp_delim(char *delim)
{
  fprintf(outfile,delim);
}


/* try an operation
   update relevant counters
   return TRUE if distance decreased by 1, FALSE o.w.
   */
   
int try_operation(int optype,
		  int d,
		  tryrevparams_t *p)
{
  int d2;


  d2 = mcdist_noncircular(p->trial_genome, p->dest_genome,
			  p->num_genes, p->num_chromosomes,
			  p->distmem,
			  FALSE);


  if (d2 == d-1) {
    if ((p->counts)[optype][CDOWN1]++ == 0) {
      /* note when the first success was found */
      p->counts[optype][CFIRST] =
	p->counts[optype][CDOWN1] +
	p->counts[optype][CSAME] +
	p->counts[optype][CUP1]	+
	p->counts[optype][CERROR] ;
      };
    return TRUE;
  } else if (d2 == d) {
    (p->counts)[optype][CSAME]++;
  } else if (d2 == d+1) {
    (p->counts)[optype][CUP1]++;
  } else {
    (p->counts)[optype][CERROR]++;
  }
  return FALSE;
}

void init_try_operation(tryrevparams_t *p)
{
  int i,j;
  /* zero counters in all categories */
  for (i=0 ; i<CNF0 ; i++)
    for (j=0 ; j<CNF1 ; j++)
      (p->counts)[i][j]=0;

}

void print_tryops_counts(tryrevparams_t *p)
{
  int i,j;

  fprintf(outfile,"\nSummary:\n"); fflush(outfile);

  /* compute count totals per operation */
  for (i=0 ; i<CNF0 ; i++) {
    p->counts[i][CTOT1] = 0;
    for (j=0 ; j<CTOT1 ; j++) 
      p->counts[i][CTOT1] += p->counts[i][j];
  }

  /* compute count totals per delta(distance) */
  for (j=0 ; j<=CTOT0 ; j++) {
    p->counts[CTOT0][j] = 0;
    for (i=0 ; i<CTOT0 ; i++)
      p->counts[CTOT0][j] += p->counts[i][j];
  }

  fprintf(outfile,
	  "              \tDown 1\tSame\tUp 1\tError\tTotal\tFirst\n");
  for(i=0 ; i<CNF0 ; i++) {
    if (p->num_chromosomes==0 && i!=CREV) continue;
    switch(i) {
    case CREV: fprintf(outfile,"Reversal:     \t"); break;
    case CTRA: fprintf(outfile,"Translocation:\t"); break;
    case CFIS: fprintf(outfile,"Fission:      \t"); break;
    case CFUS: fprintf(outfile,"Fusion:       \t"); break;
    case CTOT0: fprintf(outfile,"Total:        \t"); break;
    }
    for (j=0 ; j<CNF1 ; j++) {
      if (i != CTOT0 || j != CFIRST)
	fprintf(outfile,"%d\t",p->counts[i][j]);
    }
    fprintf(outfile,"\n");
  }
  fflush(outfile);
}
