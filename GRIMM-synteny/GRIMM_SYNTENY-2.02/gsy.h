/* gsy.h
 *    GRIMM-Synteny parameters used globally.
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

/* Last modified on Mon Sep 4, 2006, by Glenn Tesler
 */

#ifndef GSY_H
#define GSY_H

typedef struct {
  /* 0: GRIMM-Synteny
   * 1: GRIMM-Anchors
   */
  int algorithm;

  char *outputdir;


  /* GRIMM-Synteny parameters */
  /* total gap threshold */
  int gapthresh;

  /* per-species gap threshold */
  int *gapthresh_s;
  int gapthresh_s_min;  /* minimum of those */

  /* per-species minimum block size in metric 1 */
  int *min_block_size_s;
  int min_block_size_s_max;  /* maximum of those minimums */

  /* per-species minimum block size in metric 2 */
  int *min_block_size2_s;
  int min_block_size2_s_max;  /* maximum of those minimums */

  /* per-species minimum block support */
  int *min_block_sup_s;

  /* minimum # anchors per block */
  int min_block_nanc;

  /* permutation metric or nucleotide metric? */
  int perm_metric_use;  /* 0=nuc, 1=perm */
  int perm_metric_spacing;
  int perm_metric_length;

  /*
   * 1 = strips of consecutive blocks should be condensed into single block
   * 0 = don't do that
   */
  int condense_strips;

  /*
   * 1 = analyze and repair overlaps/containments
   * 0 = don't
   */
  int remove_overlap;

  /*
   * 0 = original output format:
   *     signs are one character: -,0,+
   *     mgr_micro.txt has info on blocks, compressed anchors,
   *     microrearrangements of compressed anchors,
   *     using nested begin/end pairs
   * 1 = alternate format:
   *     signs are 1-2 characters: -1,0,1
   *     mgr_micro.txt has info on blocks & compressed anchors but not
   *     microrearrangements, and it's formatted as a matrix
   */
  int format_style;
  char **signstrs;

  int sign_err_del;    /* 1: delete blocks with ambiguous signs
			*    or low correlations
			* 0: flag them but keep them
			*/
  /* minimum correlation coefficient */
  double min_corr;

#if 0
  /* reject if RS tableau diagonal length >= tableau_d */
  int tableau_d;
#endif

  /* list of cells to avoid in RS shape */
  int RS_complex_n;
  int *RS_complex_cells;

  /* GRIMM-Anchors */
  /* When signs are ambiguous, if a fraction > sign_threshold of the
   * signs go one way, choose that sign.
   */
  double sign_threshold;

  int same_species;   /* 1: species-on-self alignment coordinates */
} GSparams;


/* E_ISEDGE: (v,w) is an edge
 * E_NOTEDGE: (v,w) is not an edge
 * E_TOOHIGH: (v,w) is not an edge AND if v was < w on input, then
 *            without any more tests there is no w'>w with (v,w') an edge
 */

#define E_ISEDGE 0
#define E_NOTEDGE 1
#define E_TOOHIGH 2


#define min2(x,y) ((x)<(y) ? (x) : (y))
#define max2(x,y) ((x)>(y) ? (x) : (y))

#ifndef MAXINT
#define MAXINT (~((unsigned)1 << (8*sizeof(int) - 1)))
#endif

#endif /* GSY_H */
