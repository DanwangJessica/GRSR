/* mcstructs.h
 *    Data structures used throughout GRIMM
 *
 * Copyright (C) 2001-2006 The Regents of the University of California
 * by Glenn Tesler
 *
 * Contains code from structs.h in GRAPPA 1.02
 * Copyright (C) 2000-2001  The University of New Mexico and
 *                          The University of Texas at Austin
 * by David A. Bader, Bernard M.E. Moret, Tandy Warnow, Stacia K Wyman, Mi Yan
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

/* Modified from GRAPPA 1.02: structs.h */

#ifndef MCSTRUCTS_H
#define MCSTRUCTS_H

#include <stdio.h>
#include <stdlib.h>


#define MAX_NUM_GENES   50000
#define MAX_GENOMES     20000
#define TRUE            1
#define FALSE           0

#define MAX_STR_LEN     1024


/* for unsigned.c */
/* maximum # of singletons to exhaustively loop over
 *  (time 2^MAX_UNSIGNED_SING)
 */
/* #define MAX_UNSIGNED_SING  17 */     /* web version */
#define MAX_UNSIGNED_SING  22        /* stand alone version */


extern FILE *outfile;

#ifdef TESTING
extern double time_linear;
#endif



struct genome_struct {
  int *genes;
  int genome_num;
  char *encoding;
  char *gnamePtr; /* Used to copy the gname when adding genomes to the tree */
};



/* Structs for invdist.c */
typedef struct {
  int index;    /* index of component's root */
  int oriented; /* Boolean: Is component oriented */
  int intrachrom;
                /* Boolean: Is component intrachromosomal? */

#if 0
  int real;     /* Boolean: Is component "real"? */
#endif
  int realf;    /* Bitmask of REAL properties 
		 * Bit 0: 1=has pi-cap or gamma-tail inside, so not real
		 * Bit 1: 1=has gg path inside
		 */

  int scanstate; /* used while scanning components */

  int blocks;   /* Number of blocks in nonoriented component */
  int hurdle;   /* Bitmask of HURDLE properties.
		 * Set U (unoriented components):
		 * Bit 0 = hurdle
		 * Bit 1 = wrap-around hurdle
		 * Bit 2 = superhurdle

		 * Set IU (intrachromosomal unor. components): bits 3,4,5
		 * Set RU (real intrachrom. unor. components): bits 6,7,8
		 */

  int left;     /* Index of component to the left of my rightmost block */
  int right;    /* Index of component to the right of my rightmost block */
} component_t;

/* values of scanstate: */
#define SCAN_PENDING   0   /* haven't yet activated component */
#define SCAN_ACTIVE    1   /* activated component and have not yet
			      seen enough to determine the final answer */
#define SCAN_DONE      2   /* have final answer for this component */

/* REAL flags: */
#define NOTREAL  1
#define HASGG    (1<<1)



/*
 * sets of unoriented components:
 *         U  = all unoriented components
 *	   IU = intrachromosomal unoriented components
 *	   RU = real intrachromosomal unoriented components
 *
 * hurdles, etc.:
 *
 * set U: bit 0 hurdle     1 greatest hurdle   2 superhurdle  | fortress
 *    IU      3 knots      4 greatest knot     5 superknots   | f-of-knots
 *    RU :    6 real-knot  7 greatest real-knot8 super-real-k | f-of-real-knots
 * more bits:
 *            9 semi-knot
 *           10 simple
 */

/* hurdle flags: */
#define HURDLE          1
#define GREATHURDLE (1<<1)
#define SUPERHURDLE (1<<2)
#define HURDLE_MASK (HURDLE|GREATHURDLE|SUPERHURDLE)
#define SEMIKNOT    (1<<9)
#define SIMPLE      (1<<10)
#define U_SHIFT  0
#define IU_SHIFT 3
#define RU_SHIFT 6






/* greyEdges[i] values for non-edges */
/* Grey edge handling
 *    It was: greyEdges[i]=j with j>=0:  grey edge from i to j
 *                                j=-1:  adjacency, generally ignored
 * Multichromosome case:
 *
 *    tail pi-cap inside-verts pi-cap tail
 * and if one of the inside-verts is connected by grey edge to pi-cap,
 * the grey edge is removed and this inside-vert becomes a gamma-tail
 */
#define V_ADJ    -1    /* adjacency between 2i & 2i+1 */
#define V_TAIL   -2    /* pad at end of chromosome; delete from graph */
#define V_PICAP  -3    /* connected w/black edge to actual end of chromosome */
#define V_GTAIL  -4    /* needs connection w/grey edge to some pi-cap */

#if 0
#define V_TAILREAD -5  /* used in cap reading algorithm to re-mark V_TAIL's */
#endif


/* cycle and path types */
#define P_NONE    -1   /* vertex is NOT the leftmost one in its path/cycle */
#define P_CYCLE   -2   /* it's a cycle */
#define P_PP      -3   /* Bcap is Pi-cap,      Gcap is Pi-cap */
#define P_PG      -4   /* Bcap is Pi-cap,      Gcap is Gamma-tail */
#define P_GP      -5   /* Bcap is Gamma-tail,  Gcap is Pi-cap */
#define P_GG      -6   /* Bcap is Gamma-tail,  Gcap is Gamma-tail */



typedef struct {
  int *perm1;      /* DIST_INV: 2*num_genes + 2 */
  int *perm2;
  int *perm;
  int *done;
  int *greyEdges;

  int *ptype;      /* path/cycle type */
  int *Bcap;       /* vertex @ end of path when start @ left, follow
		    * black edge, and continue until reach the end */
  int *Gcap;       /* same, but starting @ left and following grey edge */
  int *chromNum;   /* chromosome # 1,...,N with the vertex
		    * 0 at start, N+1 at end vertex
		    */


  int *chromBd;   /* chromosome i=1,...,N (in doubled pi) occupies
		   * vertices chromBd[i],chromBd[i]+1,...,chromBd[i+1]-1
		   * size = num_chromosomes+2
		   */

  int *capMark;   /* normalized caps i=0,1,...,2*num_chromosomes-1
		   *
		   * in compute_concatenates:
		   * Array used to mark that normalized cap i has been
		   * used.
		   *
		   * in find_optimal_bond:
		   * Array used to store which caps are at opposite ends
		   * of same chromosome strings
		   */

  int *cappedp1, *cappedp2;   /* size = num_genes */
  
  int *stack;      /* DIST_INV: size */
  int *oriented;
  int *cc;
  int *labeled;
  component_t *components;
} distmem_t;


typedef struct {
  int size;             /* # vertices */
  int num_components;   /* # components */

  int components_good;  /* Boolean:
			 * TRUE if components have been computed/updated
			 * FALSE if not yet computed/updated
			 *
			 * This only describes the accuracy of the
			 * components, not their classsification flags.
			 */


  distmem_t *distmem;   /* all the structures and arrays defining the graph */
} graph_t;


typedef struct {
  int num_pp;           /* # pi-pi paths */
  int num_gg;           /* # gamma-gamma paths */
  int num_pg;           /* # pi-gamma + gamma-pi paths */
  int num_cycles;       /* # cycles other than adjacencies */
  int num_adjacencies;  /* # adjacencies = 1 black + 1 grey edge */
} pathcounts_t;


typedef struct {
  int num_hurdles;
  int num_great;
  int num_super;
  int num_fortress;

  int num_unoriented;
} hurdlecounts_t;



/* graph parameters for verbose output */
typedef struct {
  /* all */
  int n;
  int d;

  /* unichromosomal */
  int br;                /* "b" = # breakpoints */
  int c4;                /* "c" = # cycles except adjacencies */
  int h;                 /* h = # hurdles */
  int f;                 /* f = # fortresses */

  int num_components_u;  /* number of unoriented components */
  int num_components_o;  /* number of oriented components */
  /* adjacencies are not counted in either of those, but can be computed
   * as   n + 1 - br
   */

  /* multichromosomal */
  int bl;                /* "b" = # black edges */
  int cp;                /* "c" = # cycles and paths, of all lengths */
  int pgg;               /* "p" = # gamma-gamma paths */
  int s;                 /* s = # semi-knots (semi-real-knots) */
#if 0
  int rr,fr,gr;          /* parameters about G-bar */
#endif
  int r;                 /* r = # real-knots */
  int fr,gr;


  int badbonds;          /* number of bad bonds, usu. 0, sometimes 1 */
  int bp_int;            /* number of internal breakpoints */
  int bp_ext;            /* number of external breakpoints */

  /* (un)signed */

} graphstats_t;




#endif  /* #ifndef MCSTRUCTS_H */
