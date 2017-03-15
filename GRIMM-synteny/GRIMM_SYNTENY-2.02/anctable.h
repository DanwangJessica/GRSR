/* anctable.h
 *    Format of anchor table.
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

#ifndef ANCTABLE_H
#define ANCTABLE_H


/* maximum number of chromosomes per species */
#define MAXCHR 1000

typedef struct {
  int nanc;                   /* number of anchors */
  int nspecies;               /* number of species */
  int fields_per_species;     /* fields per species */
  int fields_per_anchor;      /* total # fields per anchor */

  /* anclist[k] = anchor k, 0<=k<nanc
   * anclist[k][j] = field j of anchor k
   * j = (species number)*fields_per_species + FIELDR_?
   *
   * anclist2: same format
   *
   * Clustering is done based on anclist
   *
   * If using nucleotide coordinates:
   *   anclist   = raw anchors
   *   anclist2  = NULL
   * If using permutation metric:
   *   anclist   = anchors in permutation metric
   *   anclist2  = raw anchors
   *
   * anclist3: same metric as anclist, but length field = support
   */
  int **anclist;
  int **anclist2;
  int **anclist3;

  /* orders[i][j] = k: species i, j th anchor is anclist[k] */
  int **orders;

  /* snames[i] = name of species i */
  char **snames;

  /* chrnames[i] = hash table converting between string & numeric chromosome
   * names for species i
   */
  HTABLE **chrnames;
  int max_chr;                /* maximum # chromosomes per species */
} ANCTABLE;

/* fields that are one per anchor are first
 * anclist[i][FIELD_x]
 *
 * fields that are per species come after that;
 * anclist[i][FIELD_base + species_number*FIELDR_n_per_species + FIELDR_x]
 */

/* in anchor table */
#define FIELD_id 0
/* #define FIELD_score 0 */

/* in block table */
#define FIELD_nanc 0

/* fields that are one per species then repeat after that */
#define FIELD_base 1
#define FIELDR_chr 0
#define FIELDR_start 1
#define FIELDR_len 2
#define FIELDR_sign 3
#define FIELDR_order 4
#define FIELDR_n_per_species 5

#endif /* ANCTABLE_H */
