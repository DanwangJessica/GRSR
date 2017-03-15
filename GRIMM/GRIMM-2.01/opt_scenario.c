/* opt_scenario.c
 *    Find sequence of steps of a most parsimonious scenario.
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
#include "mcrdist.h"
/*#include "e_malloc.h"*/
#include "scenario.h"
#include "testrev.h"
#include "opt_scenario.h"
#include "write_data.h"
#include "mcread_input.h"
#include "graph_edit.h"

/*#include <time.h>*/

/* construct an optimal reversal scenario
 * go from genome g1 (=gamma)  to  genome g2  (=pi)
 * prioritize the steps:
 * 1. reversals: chromosome by chromosome, shortest to longest per chromosome
 * 2. block flip
 * 3. fusion
 * 4. translocation, shortest to longest on each pair of chromosomes
 * 5. fission
 */

/* scenario_type:
 *  current values:
 *        1: old routine in scenario.c.  Respects initial capping.
 *           output as one line perm with caps.  Doesn't prioritize
 *           by type of move.
 *    2,3,4 = produce scenario that respects the initial capping
 *        2:  output as one line perms with caps
 *        3:  output as multiline perms (one per chrom) with caps
 *            / at end of each chrom
 *        4:  output as multiline perms (one per chrom) with caps omitted
 *            $ at end of each chrom
 *        5:  produce scenario that may not respect initial capping
 *            but which allows us to avoid undesireable operations.
 *            Usually use capping for efficiency.
 *            Output style is like 4.
 *        6:  use capping while can do reversals, then switch to method 7
 *        7:  try all reversals regardless of consistency with capping
 *        TODO: all translocations, fissions, fusions, too.
 */

void print_scenario_2(struct genome_struct *g1, struct genome_struct *g2,
		      int num_genes, int num_chromosomes,
		      int scenario_type)
{
  distmem_t distmem;

  struct genome_struct *dest_genome = g2;
  struct genome_struct tmp_genome1, tmp_genome2, *tmpgenome;
  /* current_genome, trial_genome alternate between tmp_genome1, tmp_genome2 */


  cbounds_t cb;        /* structure with chromosome boundaries */

  rearrangement_t Step1, Step2;
  rearrangement_t *laststep, *thisstep, *tmpstep;

  optrevparams_t Rparams;
  optrevparams_t *rparams = &Rparams;

  int d;

#if 0
  int *successors; /* array of successors in final genome */
#endif


  /* breakpoint graph on pi=current genome, gamma=dest genome */
  distmem_t distmem0;
  graph_t G0;
  pathcounts_t pcounts_G0;
  graphstats_t graphstats0;

  int multichrom = num_chromosomes>0;


  /**********************************************************************/

  init_scenario_mem(num_genes,num_chromosomes,&distmem,
		    &tmp_genome1, &tmp_genome2,
		    g1,
		    0); /* no pad */

  /* initialize successor array for determining breakpoints */
  init_scenario_bp(num_genes,num_chromosomes,
		   g2,
		   0, /* no pad */
		   &rparams->successors);

  /* initialize memory for breakpoint graph */
  mcdist_allocmem(num_genes, num_chromosomes, &distmem0);


  /**********************************************************************/
  /*                   store vars used at all steps                     */
  /**********************************************************************/
  


  /* prior step */
  laststep = &Step1;
  clearstep(laststep);

  thisstep = &Step2;


  rparams->current_genome = &tmp_genome1;
  rparams->trial_genome = &tmp_genome2;
  rparams->dest_genome = dest_genome;
  rparams->num_genes = num_genes;
  rparams->num_chromosomes = num_chromosomes;
  rparams->distmem = &distmem;
  rparams->stepno = 1;
  rparams->cb = multichrom ? &cb : (cbounds_t *) NULL;
  /*  rparams->successors = successors; */

  /* breakpoint graph between pi=current genome, gamma=dest genome */
  rparams->distmem0 = &distmem0;
  rparams->G0 = &G0;
  rparams->pcounts_G0 = &pcounts_G0;
  rparams->graphstats0 = &graphstats0;

  /* cap style */
  if (num_chromosomes == 0) {
    rparams->showcaps = 2;
  } else switch(scenario_type) {
  case 2: rparams->showcaps = 2; break;
  case 3: rparams->showcaps = 1; break;
  case 4: rparams->showcaps = 0; break;
  case 5: rparams->showcaps = 0; break;
  case 6: rparams->showcaps = 0; break;
  case 7: rparams->showcaps = 0; break;
  }

  rparams->multiline = multichrom && (scenario_type > 2);
  rparams->can_recap = multichrom &&
    (scenario_type == 5
     || scenario_type == 6
     || scenario_type == 7);
  rparams->try_all_possibilities = multichrom &&
    (scenario_type == 6
     || scenario_type == 7);
  rparams->try_all_possibilities_now = multichrom && (scenario_type == 7);


  if (multichrom) {
    rparams->cb = &cb;
    rparams->num_flip_left =
       (g1->genes[0] == g2->genes[0] &&
	g1->genes[num_genes-1] == g2->genes[num_genes-1])
      ? 0
      : 1;

    /* Assume at most one genome has null chromos, which is our
     * standardized format but could be overridden
     */
    rparams->num_fus_left = count_num_null_chromos(rparams,g2);
    rparams->num_fis_left = count_num_null_chromos(rparams,g1);
  } else {
    rparams->cb = (cbounds_t *) NULL;
    rparams->num_flip_left = 0;
    rparams->num_fus_left = 0;
    rparams->num_fis_left = 0;
  }


  /* get the chromosome boundaries
   * Note we only use cb->cBounds, and ignore cb->cNum
   */
  if (multichrom) {
    init_cbounds(num_genes,
		 num_chromosomes,
		 rparams->current_genome,
		 &cb);
  }




  /**********************************************************************/
  /*                   print an optimal scenario                        */
  /**********************************************************************/


  fprintf(outfile,"Step 0: (Source)\n");
  /*  print_genes(rparams->current_genome->genes,num_genes); */
  print_genome_multichrom_3(rparams->current_genome,
			    rparams->num_genes,
			    rparams->num_chromosomes,
			    -1,
			    rparams->multiline,
			    rparams->showcaps,
			    0);
    

  while (1) {
    if (!rparams->try_all_possibilities_now) {
      d = invdist_noncircular_G(dest_genome,                /* gamma */
				rparams->current_genome,    /* pi */
				0, /* offset */
				rparams->num_genes,
				&distmem0,
				&G0,
				&pcounts_G0,
				&graphstats0);
    } else {
      d = mcdist_noncircular(dest_genome,                /* gamma */
			     rparams->current_genome,    /* pi */
			     rparams->num_genes,
			     rparams->num_chromosomes,
			     &distmem0,
			     &graphstats0);
    }
    if (d == 0) break;
    rparams->d = d;



    opt_getstep(rparams,
		laststep,

		/* return */
		thisstep);

    if (thisstep->optype == S_FAIL) {
      /* a recapping could cause this */
      if (rparams->d == 0) break;

      fprintf(outfile,"ERROR: reversals ended prematurely\n");
    } else {
      if (!(rparams->can_recap && thisstep->optype == S_FLIP)) {
	opt_printstep(rparams,thisstep);
      } else {
#if 0
	fprintf(outfile,"Ignoring flip...\n"); /* DEBUG */
#endif
      }
    }

    if (thisstep->optype == S_FAIL) break;

    if (multichrom) {
      opt_update_cbounds(rparams, thisstep);
    }

    /* update record of prior step */
    tmpstep = thisstep;
    thisstep = laststep;
    laststep = tmpstep;

    /* update genomes */
    tmpgenome = rparams->current_genome;
    rparams->current_genome = rparams->trial_genome;
    rparams->trial_genome = tmpgenome;

    if (!(rparams->can_recap && thisstep->optype == S_FLIP)) {
      rparams->stepno++;
    }
  }


  /* clean up memory */

  if (multichrom)
    free_cbounds(&cb);

  mcdist_freemem(&distmem0);
  clean_scenario_bp_mem(rparams->successors);
  clean_scenario_mem(&tmp_genome1, &tmp_genome2, &distmem);
}


/* return values:
 * rparams->trial_genome has the successful reversal
 * thisstep has the description of the successful step
 */
void opt_getstep(optrevparams_t *rparams,
		 rearrangement_t *laststep,

		 /* return values: */
		 rearrangement_t *thisstep)
{
  int try_optype = S_FAIL;
  int try_chrom1=0;

  int stat;

  /* depending on previous step, we can skip certain steps now */
  if (rparams->num_chromosomes == 0) {
    try_optype = S_REV;
    /* and start, end... ? */
  } else {
    switch (laststep->optype) {
      /* TODO: set start, end appropriately in each */
    case S_REV:
      try_optype = S_REV;
      try_chrom1 = laststep->chr1;
      break;
    case S_FLIP:
      try_optype = S_FUS;
      try_chrom1 = 1;
      break;
    case S_FUS:
      try_optype = S_REV;
      try_chrom1 = 1;   /* TODO: the chromosome # of the fused chromosome */
      break;
    case S_TRANSLOC:
      try_optype = S_REV;
      try_chrom1 = 1;   /* TODO: the chromosome # of the left chromosome */
      break;
    case S_FIS:
      try_optype = S_FLIP;
      try_chrom1 = 1;  /* irrelevant */
      break;
    }
  }

  /* internal reversals */
  if (try_optype <= S_REV) {
    if (try_optrev(rparams, try_chrom1,
		   thisstep)) {
      return;
    }
  }

  /* translocations, fissions, fusions, flips: only in multichromo mode */

  if (rparams->num_chromosomes > 0) {
    if (rparams->can_recap && laststep->optype != S_FLIP) {
      /* recap chromos and try again */
      opt_recap(rparams);
      clearstep(laststep);
      try_optype = S_REV;
      try_chrom1 = 1;
      if (try_optrev(rparams, try_chrom1,
		     thisstep)) {
	return;
      }
    }

    /* switch to trying all possibilities (SLOW) */
    if (rparams->try_all_possibilities
	&& !rparams->try_all_possibilities_now) {
      rparams->try_all_possibilities_now = 1;
      clearstep(laststep);
      try_optype = S_REV;
      try_chrom1 = 1;
      stat = try_optrev(rparams, try_chrom1, thisstep);
      /*      rparams->try_all_possibilities_now = 0; */
      if (stat) {
	return;
      }
    }

    if (try_optype <= S_FLIP) {
      if (rparams->num_flip_left > 0
	  &&
	  try_optflip(rparams,
		      thisstep)) {
	rparams->num_flip_left --;
	return;
      }
    }

    /*    fprintf(outfile, "num_fis_left = %d, num_fus_left = %d\n", rparams->num_fis_left, rparams->num_fus_left); */

    if (rparams->num_fus_left > 0
	&&
	try_optype <= S_FUS) {
      if (try_optfus(rparams,
		     thisstep)) {
	rparams->num_fus_left --;
	return;
      }
    }

    if (try_optype <= S_TRANSLOC) {
      if (try_opttrans(rparams,
		       thisstep)) {
	return;
      }
    }

    if (rparams->num_fis_left > 0
	&&
	try_optype <= S_FIS) {
      if (try_optfis(rparams,
		     thisstep)) {
	rparams->num_fis_left --;
	return;
      }
    }

    /* Bad!  Already used the max # of fissions/fusions necessary, but
     * it wasn't compatible with this capping.  So, do an extra fusion.
     * Note: if recapping, we'll never get here.
     */

    fprintf(outfile, "Stuck!  Trying extra fusions/fissions\n");

    if (rparams->num_fus_left == 0
	&&
	try_optype <= S_FUS) {
      if (try_optfus(rparams,
		     thisstep)) {
	rparams->num_fis_left ++;
	return;
      }
    }

    if (rparams->num_fis_left == 0
	&&
	try_optype <= S_FIS) {
      if (try_optfis(rparams,
		     thisstep)) {
	rparams->num_fus_left ++;
	return;
      }
    }



  }

  /* ERROR! */
  thisstep->optype = S_FAIL;
}


void clearstep(rearrangement_t *step)
{
  step->optype = S_REV;
  step->chr1 = step->chr2 = 1;
  step->s = 0;
  step->e = -1;
}



/* try all reversals in order by chromosome, and then by increasing length */
int try_optrev(optrevparams_t *rparams,
	       int c1,  /* smallest chromosome # to try reveral in */
	       rearrangement_t *thisstep)
{
  int max_chr;

  int s1=0;     /* start and end of chromosome (caps) */
  int e1=0;

  int s1i=0;   /* internal start, end (inside of caps) */
  int e1i=0;

  int len0, maxlen0=0;
  int i, max_i=0;
  int j;

  cbounds_t *cb = rparams->cb;

  /* params needed for ISUCCESSOR macro to work */
  int num_genes = rparams->num_genes;
  /*  int num_chromosomes = rparams->num_chromosomes; */
  int *successors = rparams->successors;
  int gene;

  int *current_genome_genes = rparams->current_genome->genes;
  int *dest_genome_genes = rparams->dest_genome->genes;


  max_chr = rparams->num_chromosomes;
  if (max_chr == 0) {
    c1 = max_chr = 1;
  }


  /* test chromosomes in order */
  while (c1 <= max_chr) {
    /* start & end of chromosome */
    if (rparams->num_chromosomes > 0) {
      s1 = cb->cBound[c1-1];   /* lcap */
      e1 = cb->cBound[c1]-1;   /* rcap */

      /* internal start, end, that omits caps */
      s1i = s1+1;
      e1i = e1-1;

      /* disallow internal chromosome flip by setting
       * length 1 less than full length of chromosome */
      maxlen0 = e1i-s1i;
    } else {
      s1i = s1 = 0;
      e1i = e1 = rparams->num_genes - 1;
      maxlen0 = rparams->num_genes;
    }

    /* reversal length = j-i+1 = len0 + 1 */
    for (len0=0 ; len0 < maxlen0 ; len0++) {
      max_i = e1i-len0;
      for (i=s1i ; i <= max_i ; i++) {
	j = i + len0;


	/***************************************
	 * only start and end at breakpoints
	 ***************************************/

	/* don't split signed 2-strip at start of reversal */
	if (i==0) {
	  /* if genomes start w/same gene, don't reverse it */
	  if (current_genome_genes[i] == dest_genome_genes[i]) continue;
	} else {
	  /* don't start here if it breaks adjacency */
	  if (!rparams->try_all_possibilities_now || i>s1i) {
	    gene = current_genome_genes[i-1];
	    if (SUCCESSOR(gene) == current_genome_genes[i])
	      continue;
	  }
	}

	/* don't split signed 2-strip at end of reversal */
	if (j==num_genes-1) {
	  /* if genomes end w/same gene, don't reverse it */
	  if (current_genome_genes[j] == dest_genome_genes[j]) continue;
	} else {
	  /* don't end here if it breaks adjacency */
	  if (!rparams->try_all_possibilities_now || j<e1i) {
	    gene = current_genome_genes[j];
	    if (SUCCESSOR(gene) == current_genome_genes[j+1])
	      continue;
	  }
	}


	/***************************************
	 * try the reversal
	 ***************************************/

	if (optrev_try_op(rparams, i, j)) {
	  thisstep->optype = S_REV;
	  thisstep->s = i;
	  thisstep->e = j;
	  thisstep->chr1 = thisstep->chr2 = c1;
	  return 1;
	}

      } /* end for i (start of reversal) */
    } /* end for len0 (length of reversal) */

    /* try next chromosome */
    c1++;
  }

  thisstep->optype = S_FAIL;
  return 0;
}


/* try all translocations in order by pairs of chromosomes,
 * and then by increasing length */
int try_opttrans(optrevparams_t *rparams,
		 rearrangement_t *thisstep)
{
  int s1=0, e1=0;     /* start and end of chromosome (caps) */
  int s1i=0, e1i=0;   /* internal start, end (inside of caps) */
  int s2=0, e2=0, s2i=0, e2i=0;
  int c1,c2;

  int off, maxoff, minoff;
  int i, max_i=0;
  int j;

  cbounds_t *cb = rparams->cb;

  /* params needed for ISUCCESSOR macro to work */
  int num_genes = rparams->num_genes;
  int num_chromosomes = rparams->num_chromosomes;
  int *successors = rparams->successors;
  int gene;

  int *current_genome_genes = rparams->current_genome->genes;
  /*  int *dest_genome_genes = rparams->dest_genome->genes; */


  /* test chromosome pairs in order */
  for (c1 = 1; c1 < num_chromosomes ; c1++) {

    /* start & end of chromosome */
    s1 = cb->cBound[c1-1];   /* lcap */
    e1 = cb->cBound[c1]-1;   /* rcap */

    /* disallow translocations with null chromosome c1
     * (they are fissions, considered separately)
     */
    /*    fprintf(outfile,"A1    c1=%d s1=%d e1=%d\n",c1,s1,e1); */
    if (e1 - s1 == 1) continue;

    /* internal start, end, that omits caps */
    s1i = s1+1;
    e1i = e1-1;

    /* 2nd chromosome */
    for (c2 = c1+1 ; c2 <= num_chromosomes ; c2++) {
      /* start & end of chromosome */
      s2 = cb->cBound[c2-1];   /* lcap */
      e2 = cb->cBound[c2]-1;   /* rcap */

      /* disallow translocations with null chromosome c2
       * (they are fissions, considered separately)
       */
      /*            fprintf(outfile,"A2    c2=%d s2=%d e2=%d\n",c2,s2,e2); */
      if (e2 - s2 == 1) continue;

      
      /* internal start, end, that omits caps */
      s2i = s2+1;
      e2i = e2-1;

      /* disallow translocation (s1i,e2i), which is a cap exchange */
      maxoff = e2i-s1i;

      /* minimum translocation shifts one gene */
      minoff = s2-e1+1;

      /*    fprintf(outfile,"A3    minoff=%d maxoff=%d\n",minoff,maxoff); */

      for (off=minoff ; off < maxoff ; off++) {
	max_i = e2i-off;
	if (max_i > e1) max_i = e1;

	/*	fprintf(outfile,"A4    off=%d max_i=%d\n",off,max_i); */
	for (i=s2-off ; i <= max_i ; i++) {
	  if (i <= s1) continue;
	  j = i + off;

	  /*	    fprintf(outfile,"A5    i=%d j=%d\n",i,j); */

	  /***************************************
	   * only start and end at breakpoints
	   ***************************************/

	  /* don't split signed 2-strip at start of reversal */
	  /* don't start here if it breaks adjacency */
	  gene = current_genome_genes[i-1];
	  if (SUCCESSOR(gene) == current_genome_genes[i])
	    continue;

	  /* don't split signed 2-strip at end of reversal */
	  /* don't end here if it breaks adjacency */
	  gene = current_genome_genes[j];
	  if (SUCCESSOR(gene) == current_genome_genes[j+1])
	    continue;


	  /**************************************************
	   * disallow translocations that are actually fusion
	   * as these should be tested for separately
	   **************************************************/

	  if ((i==s1i && j==s2) || (i==e1 && j==e2i)) continue;

	  /*	  fprintf(outfile,"A6    i=%d j=%d\n",i,j); */

	  /***************************************
	   * try the translocation
	   ***************************************/

	  if (optrev_try_op(rparams, i, j)) {
	    thisstep->optype = S_TRANSLOC;
	    thisstep->s = i;
	    thisstep->e = j;
	    thisstep->chr1 = c1;
	    thisstep->chr2 = c2;
	    return 1;
	  }

	} /* end for i (start of reversal) */
      } /* end for len0 (length of reversal) */
    } /* end for c2 */
  } /* end for c1 */

  thisstep->optype = S_FAIL;
  return 0;
}


/* try all fissions */
int try_optfis(optrevparams_t *rparams,
	       rearrangement_t *thisstep)
{
  int s1=0, e1=0;     /* start and end of chromosome (caps) */
  int s1i=0, e1i=0;   /* internal start, end (inside of caps) */
  int s2=0, e2=0;
  int c1,c2;

  int i;
  int j;

  cbounds_t *cb = rparams->cb;

  /* params needed for ISUCCESSOR macro to work */
  int num_genes = rparams->num_genes;
  int num_chromosomes = rparams->num_chromosomes;
  int *successors = rparams->successors;
  int gene;

  int *current_genome_genes = rparams->current_genome->genes;

  /* c1 is a nonnull chromosome, c2 is a null chromosome */
  for (c1 = 1; c1 <= num_chromosomes ; c1++) {
    /* start & end of chromosome */
    s1 = cb->cBound[c1-1];   /* lcap */
    e1 = cb->cBound[c1]-1;   /* rcap */

    /* make sure c1 is nonnull */
    if (e1 - s1 == 1) continue;

    /* internal start, end, that omits caps */
    s1i = s1+1;
    e1i = e1-1;

    /* 2nd chromosome */
    for (c2 = 1 ; c2 <= num_chromosomes ; c2++) {
      if (c1 == c2) continue;

      /* start & end of chromosome */
      s2 = cb->cBound[c2-1];   /* lcap */
      e2 = cb->cBound[c2]-1;   /* rcap */

      /* make sure c2 is null */
      if (e2 - s2 != 1) continue;

      /* make sure point between caps of c2 is a breakpoint */
      gene = current_genome_genes[s2];
      if (SUCCESSOR(gene) == current_genome_genes[e2]) continue;

      if (c1<c2) {
	for (i=s1i+1 ; i<=e1i ; i++) {
	  /* make sure b.p. at position i */
	  gene = current_genome_genes[i-1];
	  if (SUCCESSOR(gene) == current_genome_genes[i])
	    continue;

	  /* try the fission */
	  if (optrev_try_op(rparams, i, s2)) {
	    thisstep->optype = S_FIS;
	    thisstep->s = i;
	    thisstep->e = s2;
	    thisstep->chr1 = c1;
	    thisstep->chr2 = c2;
	    return 1;
	  }
	}
      } else {
	for (j=s1 ; j<e1i ; j++) {
	  /* make sure b.p. at position j */
	  gene = current_genome_genes[j];
	  if (SUCCESSOR(gene) == current_genome_genes[j+1])
	    continue;

	  /* try the fission */
	  if (optrev_try_op(rparams, e2, j)) {
	    thisstep->optype = S_FIS;
	    thisstep->s = e2;
	    thisstep->e = j;
	    thisstep->chr1 = c2;
	    thisstep->chr2 = c1;
	    return 1;
	  }
	}

      }

    } /* end for c2 */

  } /* end for c1 */

  thisstep->optype = S_FAIL;
  return 0;
}

/* attempt to do a block flip
 * our capping algorithm only permits prefix or suffix block flip
 */
int try_optflip(optrevparams_t *rparams,
		rearrangement_t *thisstep)
{
  int *current_genome_genes = rparams->current_genome->genes;
  int *dest_genome_genes = rparams->dest_genome->genes;

  int num_chromosomes = rparams->num_chromosomes;
  int num_genes = rparams->num_genes;

  int lowcap = num_genes - 2*num_chromosomes + 1;
  int gene0;
  int g;

  int ncaps = 0;
  int chr1 = 1;
  int chr2 = num_chromosomes;

  int i,j;

  /* prefix or suffix? */
  gene0 = dest_genome_genes[0];
  if (gene0 != current_genome_genes[0]) {
    /* look for prefix reversal */
    gene0 = -gene0;
    i = 0;
    for (j=i; j<num_genes; j++) {
      g = current_genome_genes[j];
      if (g >= lowcap || g <= -lowcap) ncaps++;
      if (g == gene0) break;
    }
    chr2 = ncaps/2;
  } else {
    /* look for suffix reversal */
    j = num_genes-1;
    gene0 = -dest_genome_genes[j];
    for (i=j; i>=0; i--) {
      g = current_genome_genes[i];
      if (g >= lowcap || g <= -lowcap) ncaps++;
      if (g == gene0) break;
    }
    chr1 = num_chromosomes - ncaps/2 + 1;
  }

  if (optrev_try_op(rparams,i,j)) {
	  thisstep->optype = S_FLIP;
	  thisstep->s = i;
	  thisstep->e = j;
	  thisstep->chr1 = chr1;
	  thisstep->chr2 = chr2;
	  return 1;
  }

  return 0;
}


/* test fusions */
int try_optfus(optrevparams_t *rparams,
	       rearrangement_t *thisstep)
{
  int num_chromosomes = rparams->num_chromosomes;
  int c1, c2;
  int s1=0,e1=0,s2=0,e2=0;
  cbounds_t *cb = rparams->cb;

  for (c1 = 1; c1 < num_chromosomes; c1++) {
    s1 = cb->cBound[c1-1];   /* lcap of c1 */
    e1 = cb->cBound[c1]-1;   /* rcap of c1 */
    
    for (c2 = c1+1 ; c2 <= num_chromosomes; c2++) {
      s2 = cb->cBound[c2-1];   /* lcap of c2 */
      e2 = cb->cBound[c2]-1;   /* rcap of c2 */


      if (optrev_try_op(rparams, s1+1, s2)) {         /* fusion (-c1) (c2) */
	thisstep->optype = S_FUS;
	thisstep->chr1 = c1;
	thisstep->chr2 = c2;
	thisstep->s = s1+1;
	thisstep->e = s2;
	return 1;
      } else if (optrev_try_op(rparams, e1, e2-1)) {  /* fusion (c1) (-c2) */
	thisstep->optype = S_FUS;
	thisstep->chr1 = c1;
	thisstep->chr2 = c2;
	thisstep->s = e1;
	thisstep->e = e2-1;
	return 1;
      }
    }
  }

  return 0;
}



/* try the reversal from positions i to j
 * input: various fields of rparams
 *      rparams->current_genome, rparams->dest_genome are the genomes
 *      rparams->d is reversal distance from current_genome to dest_genome
 * output: trial_genome is the reversal (i,j) applied to current_genome
 * return value: 0 = doesn't reduce distance, 1 = does reduce distance
 */

int optrev_try_op(optrevparams_t *rparams,
		  int i,
		  int j)
{
  int d2;
  graphstats_t *graphstats0 = rparams->graphstats0;
  /*  int *cycle = rparams->distmem0->labeled; */
  /*  int e; */

  if (rparams->try_all_possibilities_now) {
    /* brute force: do the reversal */
    copy_and_reverse_genes(rparams->current_genome->genes,
			   rparams->trial_genome->genes,
			   i,j,
			   rparams->num_genes);


    d2 = mcdist_noncircular(rparams->trial_genome,
			    rparams->dest_genome,
			    rparams->num_genes,
			    rparams->num_chromosomes,
			    rparams->distmem,
			    (graphstats_t *)0);
    /*
    if (d2 < rparams->d) {
      fprintf(outfile,"try_all_pos: d2=%d d=%d badbonds=%d\n",d2,rparams->d,rparams->graphstats0->badbonds);
    }
    */
    return (d2 == rparams->d - 1) ? 1 : 0;
  }

  /* tests on breakpoint graph that are valid w/o hurdles */
  if (graphstats0->h == 0) {
    /* the two ends of the reversal must be in the same cycle, an
     * even number of vertices apart
     */

    /* check same cycle:
     * cycle indices agree and are >=0 
     * (if agree but <0, it's actually different cycles)
     */
    /*    if ((e=cycle[2*i+1]) != cycle[2*j+2] || e<0) return 0; */

    /* check if an even number of vertices apart */

    /* TODO:
     * compare two methods:
     *  1. computing distance for each i,j as necessary;
     *  2. precomputing positions in cycles (color 0,1,0,1,...)
     */
    if (graphdist(rparams->G0, 2*i+1, 2*j+2) % 2 != 0) return 0;

  }

  /* brute force: do the reversal */
  copy_and_reverse_genes(rparams->current_genome->genes,
			 rparams->trial_genome->genes,
			 i,j,
			 rparams->num_genes);


  d2 = invdist_noncircular(rparams->trial_genome,
			   rparams->dest_genome,
			   0, /* offset */
			   rparams->num_genes,
			   rparams->distmem);

  return (d2 == rparams->d - 1) ? 1 : 0;
}


static char *stepnames[] = {
  "Reversal",
  "Trivial chromosome flip",
  "Fusion",
  "Translocation",
  "Fission",
  "Illegal reversal",
  "Cap exchange",
  "FAIL"
};

void opt_printstep(optrevparams_t *rparams,
		   rearrangement_t *step)
{
  int chr1 = step->chr1;
  int chr2 = step->chr2;
  cbounds_t *cb = rparams->cb;

  int chr1s = step->s;
  int chr2e = step->e;

  int *trial_genome_genes = rparams->trial_genome->genes;
  int startvalue, endvalue;

  /*  time_t clock; */

  startvalue = -trial_genome_genes[step->e];
  endvalue = -trial_genome_genes[step->s];

  if (cb) {
    chr1s -= cb->cBound[chr1-1];
    chr2e -= cb->cBound[chr2-1];
  } else {
    chr1 = 0;
    chr2 = 0;
  }

  /*
   time(&clock);
   fprintf(outfile, "%s", ctime(&clock));
  */

  fprintf(outfile, "Step %d: ", rparams->stepno);

  fprintf(outfile,
	  "Chrom. %d, gene %d [%d] through chrom. %d, gene %d [%d]: ",
	  chr1,
	  chr1s,
	  startvalue,
	  chr2,
	  chr2e,
	  endvalue);
  
  fprintf(outfile, "%s", stepnames[step->optype]);

  if (rparams->d <= 1)
    fprintf(outfile, " (Destination)");

  fprintf(outfile, "\n");

  /*  print_genes(rparams->trial_genome->genes, rparams->num_genes); */


  print_genome_multichrom_3(rparams->trial_genome,
			    rparams->num_genes,
			    rparams->num_chromosomes,
			    -1,
			    rparams->multiline,
			    rparams->showcaps,
			    0);
}


int count_num_null_chromos(optrevparams_t *rparams,
			   struct genome_struct *g1)
{
  int lowcap = rparams->num_genes - 2*rparams->num_chromosomes + 1;
  int prev_is_cap = 0;
  int nnulls = 0;
  int inchrom = 0;
  int g, i;

  int *genes = g1->genes;

  int num_genes = rparams->num_genes;

  for (i=0; i<num_genes; i++) {
    g = genes[i];
    if (g >= lowcap || g <= -lowcap) {
      inchrom = !inchrom;
      if (!inchrom && prev_is_cap) nnulls++;
      prev_is_cap = 1;
    } else {
      prev_is_cap = 0;
    }
  }

  /*  fprintf(stderr,"Num nulls = %d\n", nnulls); */
  return nnulls;
}

/* count number of null chromos with standard capping */
int count_num_null_chromos_std(optrevparams_t *rparams,
			       struct genome_struct *g1)
{
  int nnulls = 0;

  int *genes = g1->genes;

  int num_genes = rparams->num_genes;
  while (num_genes > 0 &&
	 genes[num_genes-2] == num_genes - 1) {
    nnulls++;
    num_genes -= 2;
  }

  return nnulls;
}

void opt_update_cbounds(optrevparams_t *rparams,
			rearrangement_t *step)
{
  int c1 = step->chr1;
  int c2 = step->chr2;
  int s  = step->s;
  int e  = step->e;

  /* mirror image points by subtracting from this
   * if the lcap was at i, it is reflected into an rcap at s+e-i;
   * thus the lcap of the next chromo is at s+e+1 - i
   */
  int m  = s+e+1;        
  int temp;

  int d1,d2;


  int *cBound = rparams->cb->cBound;


  /* update starts of chromos c1+1,c1+2,...,c2
   * note cBound[c-1] is start of chromo c
   */

  /* first swap the boundaries */
  for (d1 = c1 , d2 = c2-1 ;
       d1 < d2 ;
       d1++, d2--) {
    temp = cBound[d2];
    cBound[d2] = cBound[d1];
    cBound[d1] = temp;
  }

  /* then mirror image them */
  for (d1=c1 ; d1<c2 ; d1++) {
    cBound[d1] = m - cBound[d1];
  }
}

void opt_recap(optrevparams_t *rparams)
{
  int d;

  reformat_caps(rparams->dest_genome,
		rparams->num_genes,
		rparams->num_chromosomes);

  reformat_caps(rparams->current_genome,
		rparams->num_genes,
		rparams->num_chromosomes);

  opt_strip_endnulls(rparams);

  rparams->num_fus_left = count_num_null_chromos_std(rparams,
						     rparams->dest_genome);
  rparams->num_fis_left = count_num_null_chromos_std(rparams,
						     rparams->current_genome);

  /*  fprintf(outfile, "recap: num_fis_left = %d, num_fus_left = %d\n", rparams->num_fis_left, rparams->num_fus_left); */

  /* recap the genomes */
  d = mcdist_capgraph_nomem(rparams->dest_genome,
			    rparams->current_genome,
			    rparams->num_genes,
			    rparams->num_chromosomes,
			    rparams->graphstats0);

  /* recompute other quantities */

  /* recompute successor array for dest genome */
  /*  clean_scenario_bp_mem(rparams->successors); */
  init_scenario_bp_wmem(rparams->num_genes,
			rparams->num_chromosomes,
			rparams->dest_genome,
			0, /* no pad */
			&rparams->successors);

  /* recompute chromosome bounds for current genome */
  /*  free_cbounds(rparams->cb); */
  init_cbounds_wmem(rparams->num_genes,
		    rparams->num_chromosomes,
		    rparams->current_genome,
		    rparams->cb);


  /* recompute breakpoint graph on _permutations_ */
  if (!rparams->try_all_possibilities_now) {
    d = invdist_noncircular_G(rparams->dest_genome,       /* gamma */
			      rparams->current_genome,    /* pi */
			      0, /* offset */
			      rparams->num_genes,
			      rparams->distmem0,
			      rparams->G0,
			      rparams->pcounts_G0,
			      rparams->graphstats0);
  }
  rparams->d = d;
}


/* assume genomes are formatted with consecutive cap #s,
 * and all nulls at end.
 * Strip off common nulls at end.
 */
void opt_strip_endnulls(optrevparams_t *rparams)
{
  int num_genes = rparams->num_genes;
  int num_chromosomes = rparams->num_chromosomes;
  int *dest_genome_genes = rparams->dest_genome->genes;
  int *current_genome_genes = rparams->current_genome->genes;

  while (num_genes > 0 &&
	 dest_genome_genes[num_genes-2] == num_genes-1 &&
	 current_genome_genes[num_genes-2] == num_genes-1) {
    num_genes -= 2;
    num_chromosomes--;
  }

  rparams->num_genes = num_genes;
  rparams->num_chromosomes = num_chromosomes;
}

