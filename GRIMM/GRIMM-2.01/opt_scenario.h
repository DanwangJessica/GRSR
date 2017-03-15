/* opt_scenario.h
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

#ifndef OPT_SCENARIO_H
#define OPT_SCENARIO_H

/* prioritize the step types */
#define S_REV      0
#define S_FLIP     1
#define S_FUS      2
#define S_TRANSLOC 3
#define S_FIS      4

/* define the other operations, even though we don't allow them */
#define S_ILL_REV   5    /* illegal reversal (bad combination of caps) */
#define S_CAP_EXCH  6    /* cap exchange */ 
#define S_FAIL  7

/* describe a step */

typedef struct {
  int optype;        /* operation type */
  int s, e;          /* start, end of reversal that emulates it */
  int chr1, chr2;    /* the chromosome number(s) involved */
} rearrangement_t;


typedef struct {
  struct genome_struct *current_genome;
  struct genome_struct *trial_genome;
  struct genome_struct *dest_genome;
  int num_genes;
  int num_chromosomes;
  distmem_t *distmem;

  cbounds_t *cb;
  int *successors;

  int d;                        /* rev distance from current_genome
				 * to dest_genome
				 */
  int stepno;


  /* breakpoint graph */
  distmem_t *distmem0;
  graph_t *G0;
  pathcounts_t *pcounts_G0;
  graphstats_t *graphstats0;


  /* interchromosomal steps:
   * count how many of each type are left (when known)
   */
  int num_flip_left;
  int num_fus_left;
  int num_fis_left;
  /* int num_transloc_left; */
  /* int_num_rev_left; */

  /* cap style */
  int showcaps;
  int multiline;

  /* can we recap for better steps, or are we constrained to
   * use initial cappping/concatenation? */
  int can_recap;

  /* should we try all possibilities, or only those consistent with
   * current capping?
   */
  int try_all_possibilities;       /* yes in general... */
  int try_all_possibilities_now;   /* right now */



} optrevparams_t;


void print_scenario_2(struct genome_struct *g1, struct genome_struct *g2,
		      int num_genes, int num_chromosomes,
		      int scenario_type);


void opt_getstep(optrevparams_t *rparams,
		 rearrangement_t *laststep,

		 /* return values: */
		 rearrangement_t *thisstep);

void clearstep(rearrangement_t *step);

int try_optrev(optrevparams_t *rparams,
	       int c1,  /* smallest chromosome # to try reveral in */
	       rearrangement_t *thisstep);

int try_optflip(optrevparams_t *rparams,
		rearrangement_t *thisstep);

int try_optfus(optrevparams_t *rparams,
	       rearrangement_t *thisstep);

int optrev_try_op(optrevparams_t *rparams,
		  int i,
		  int j);

int try_opttrans(optrevparams_t *rparams,
		 rearrangement_t *thisstep);

int try_optfis(optrevparams_t *rparams,
	       rearrangement_t *thisstep);


void opt_printstep(optrevparams_t *rparams,
		   rearrangement_t *step);


int count_num_null_chromos(optrevparams_t *rparams,
			   struct genome_struct *g1);

void opt_update_cbounds(optrevparams_t *rparams,
			rearrangement_t *step);


void opt_recap(optrevparams_t *rparams);

void opt_strip_endnulls(optrevparams_t *rparams);

#endif /* OPT_SCENARIO_H */
