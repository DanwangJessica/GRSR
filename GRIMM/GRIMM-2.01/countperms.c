/* countperms.c
 *    Count how many permutations have certain distances & breakpoint stats.
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
/* Fix n
   Do unichromosomal alg on all signed perms of length n,
   counting how many achieve each distance.

   TODO:
   count how many achieve each combination of various parameters:
     distance
     # cycles
     h
     f
     # components
     # (un)oriented components
     # max components

   option t:
     # | d c2 h f c4 br ku ko | t(1) t(2) t(3) t(4) ... t(n+1)
     # = count of how many times that tuple occurs
     t(1),...,t(n+1): decreasing list of comp. sizes (size=# black edges)
     a = n+1 - br
     c2 = c4 + a = c4 + n+1 - br
     c4 = c2 + br - (n+1)
     ku = # unoriented components
     ko = # oriented components
     (adjacencies are separate from ku and ko)

   TODO:
   variants for
     unsigned
     circular
     multichromosomal

   Did unsigned BUT not with other options
*/

#include <stdio.h>
#include <stdlib.h>
#include <search.h>
#include <string.h>
#include "uniinvdist.h"
#include "mcstructs.h"
#include "scenario.h"
#include "mcrdist.h"
#include "unsigned.h"
#include "graph_edit.h"
#include "graph_components.h"
#include "e_malloc.h"
/*#include "circ_align.h"*/
#include "countperms.h"



typedef struct {
  int reclen;
  void *root;
} CTHASH;
CTHASH *cthash_create(int reclen);
int *cthash_alloc(CTHASH *h);
int *cthash_insert(CTHASH *h, int *rec);
int cthash_cmp(const void *s1, const void *s2);
void ctype_calc(int num_genes, distmem_t *distmem, int *ctype_mem);
void ctprint(FILE *stream, int n, const int *p);
void cthash_print(CTHASH *cthash);
void cthash_print_node(const void *node, VISIT order, int level);
int ctype_ptn_cmp(const void *a, const void *b);







int inc_sign_simple(int *g, int n);
int inc_sign_findadj(int *g, int n);
int inc_sign_adj(int *g, int n);
int inc_sign(int *g, int n, int adj_only);
int inc_perm(int *g, int n);



void countperms(char *options, int circular, int unsigned_dist) {
  struct genome_struct g1_mem, g2_mem, *g1=&g1_mem, *g2=&g2_mem;
  distmem_t distmem_mem, *distmem=&distmem_mem;
  int num_chromosomes=0;

  graphstats_t graphstats_mem, *graphstats=&graphstats_mem;


  /* options */
  char *s = options, c;
  int num_genes = 0;
  int adj_only = 0;

  int maxdist;
  int *results_d;
  int d;

  int i;
  int r;
  int unsigned_dist_err;

  int do_ctype = 0;
  CTHASH *cthash = (CTHASH *) NULL;
  int ctype_len = 0;
  int *ctype_mem = (int *) NULL;
  int *ctype_rec = (int *) NULL;
  int index;




  /* parse options */
  while ((c=*s++) != '\0') {
    switch (c) {
    case '0': case '1': case '2': case '3': case '4':
    case '5': case '6': case '7': case '8': case '9':
      num_genes = 10*num_genes + (c-'0');
      break;

    case 'a':
      adj_only = 1;
      break;

    case 't':
      do_ctype = 1;
      break;
    }
  }

  fprintf(outfile,
	  "Options: num_genes=%d, adj_only=%d, circular=%d, unsigned=%d, do_ctype=%d\n",
	  num_genes, adj_only, circular, unsigned_dist, do_ctype);

  maxdist = num_genes+2; /* higher than should be */



  /* allocate memory */

  /* memory for breakpoint graph computations */
  mcdist_allocmem(num_genes, num_chromosomes,
		  distmem);

  /* memory for permutations */
  alloc_simple_genome(num_genes, g1);
  alloc_simple_genome(num_genes, g2);

  /* memory for recording results: distance counts */
  /* TODO: structure to store k-tuples of joint distrib. of statistics */
  results_d = (int *) e_malloc((maxdist+1) * sizeof(int),
			       "results_d (countperms)");
  for (i=0; i<=maxdist; i++)
    results_d[i] = 0;

  /* memory for recording results: component types */
  if (do_ctype) {
    /* structure in which to store results */
    /* each record contains these items, stored as array of integers
     *    value:   count
     *    key:     d c2 h f c4 br  ku ko  t(1) ... t(n+1)
     */
    cthash = cthash_create(9 + num_genes+1);

    ctype_len = (num_genes+1) * sizeof(int);
    ctype_mem = (int *) e_malloc(ctype_len, "ctype_mem");
  }


  /* g1=gamma=identity, g2=pi is looped over
   * init both to identity
   */
  for (i=0; i<num_genes; i++) {
    g1->genes[i] = g2->genes[i] = i+1;
  }
  
  while (1) {

    /* TODO:
     * increase efficiency
     * currently, have this here and more at the end of loop;
     * and have extra calls to inc_sign_findadj
     */
    if (adj_only) {
      r = 1;
      while (inc_sign_findadj(g2->genes, num_genes)>=0 && r) {
	r = inc_sign(g2->genes, num_genes, adj_only);
	if (!r) r=inc_perm(g2->genes, num_genes);
      }
    }


    /* debugging */
#if 0
    for (i=0; i<num_genes; i++) {
      fprintf(outfile, "%d ", g2->genes[i]);
    }
    fprintf(outfile, "\n");
#endif



    if (unsigned_dist) {
      d = unsigned_mcdist(g1, g2,
			  num_genes,
			  num_chromosomes,
			  circular,
			  distmem,
			  0,
			  &unsigned_dist_err);
    } else {
      d = invdist_noncircular_v(g1, g2,
				0, /* offset */
				num_genes, distmem, graphstats);

      /* record results
	 TODO: add more statistics!
      */
      if (do_ctype) {
	/* allocate new record */
	ctype_rec = cthash_alloc(cthash);
	index = 1;
	ctype_rec[index++] = d;
	ctype_rec[index++] = graphstats->c4 + (num_genes+1) - graphstats->br; /* c2 */
	ctype_rec[index++] = graphstats->h;
	ctype_rec[index++] = graphstats->f;
	ctype_rec[index++] = graphstats->c4;
	ctype_rec[index++] = graphstats->br;

	ctype_rec[index++] = graphstats->num_components_u;
	ctype_rec[index++] = graphstats->num_components_o;

	ctype_calc(num_genes, distmem, ctype_mem);
	memcpy(&ctype_rec[index], ctype_mem, ctype_len);

	/* insert record and destroy original if necessary */
	cthash_insert(cthash, ctype_rec);
      }
    }

    /* should never happen! */
    if (d>maxdist) d=maxdist;

    results_d[d]++;



    /* advance to next permutation */
    if (unsigned_dist) {
      if (!inc_perm(g2->genes,num_genes))
	break;
    } else {
      if (!inc_sign(g2->genes,num_genes,adj_only)
	  && !inc_perm(g2->genes,num_genes))
	break;
    }
  }

  /* show results */
  for (i=0; i<=maxdist; i++)
    fprintf(outfile,"%d %d\n",i,results_d[i]);

  if (do_ctype) {
    cthash_print(cthash);
  }


  /* free mem */

  if (do_ctype) {
    /* TODO:
     * deallocate cthash
     */

    free((void *) ctype_mem);
  }

  free((void *) results_d);
  free_simple_genome(g2);
  free_simple_genome(g1);
  mcdist_freemem(distmem);
}


/* advance signs in current permutation
 * return:
 *   0: done (all are back to +)
 *   1: not done
 */

int inc_sign_simple(int *g, int n)
{
  int i, x;

#if 0
  fprintf(outfile, "+ "); /* debugging */
#endif

  /* negate g[n-1], g[n-2], ..., g[k]
   * where g[n-1], g[n-2], ..., g[k+1] are positive
   * and g[k] is negative
   */
  for (i=n-1; i>=0; i--) {
    x = g[i];
    g[i] = -x;
    if (x>0) return 1;
  }
  return 0;
}


/*
 * does g contain a pattern
 * (1)  g[p-1]=g[p]-1    for some p=1,...,n-1
 * (2)  g[0]=1           so p=0
 * (3)  g[n-1]=n         so p=n
 * If so, return p
 * If not, return -1
 *
 * OLD:
 * does g contain pattern
 * (i)    0<g[p-1]=g[p]-1 &&  g[p+1],...,g[n-1]>0   for some p>=1
 * or
 * (ii)   g[n-1]==n                                 (so p=n)
 * or
 * (iii)  0>g[p-1]=g[p]-1 &&  g[p+1],...,g[n-1]>0   for some p>=1
 * or
 * (iv)   g[0]==1 && g[1],...,g[n-1]>0              (so p=0)
 * if so, return p.
 * if not, return -1.
 */
int inc_sign_findadj(int *g, int n)
{
  int p=n-1;
  if (g[p]==n) return n;      /* pattern (3) */

  while (p>0) {
    if (g[p-1] == g[p]-1) return p;              /* pattern (1) */
#if 0
    //    if (g[p]>0 && g[p-1]==g[p]-1) return p;   /* pattern (i) */
    //    if (g[p]<0 && g[p-1]==g[p]-1) return p;   /* pattern (iii) */
    //    if (g[p]<0) return -1;
#endif /* 0 */
    p--;
  }
  if (g[p]==1) return 0;      /* pattern (2) */

  return -1;
}

int inc_sign_adj(int *g, int n)
{
  int r;
  int p=0;

  r = inc_sign_simple(g,n);
  if (!r)
    return r;

  /* is there an adjacency? */
  while (p>=0 && r) {
    p = inc_sign_findadj(g,n);
    if (p<0)
      return r; /* no adjacency! */
    r = inc_sign_simple(g,(p==n) ? n : p+1);
    /* set later signs to + */
    p++;
    while (p<n) {
      if (g[p]<0) g[p]=-g[p];
      p++;
    }
  }

  return r;
} 

/* advance signs in current permutation
 * if adj_only=1, must advance in a way that avoids signed adjacencies
 */
int inc_sign(int *g, int n, int adj_only)
{
  if (adj_only)
    return inc_sign_adj(g,n);
  return inc_sign_simple(g,n);
}


/* advance current permutation
 * return:
 *   0: done (back to identity)
 *   1: not done
 */
int inc_perm(int *g, int n)
{
  int i, j, k, r, temp;

#if 0
  fprintf(outfile, "s "); /* debugging */
#endif

  /* permutation is
   * ..., g[k], g[k+1], g[k+2], g[k+3], ..., g[k+r], ..., g[n]
   * where g[k] < g[k+1] > g[k+2] > g[k+3] > ... > g[n]
   * and r>0 is least with g[k+r]>g[k]
   * 
   * change it to
   * ..., g[k+r], g[n],g[n-1],...,g[k+r+1], g[k], g[k+r-1],...,g[k+1]
   */

  /* find rightmost ascent */
  for (k=n-2; k>=0; k--) {
    if (abs(g[k]) < abs(g[k+1])) break;
  }

#if 0
  fprintf(outfile, "s%d ", k); /* debugging */
#endif


  if (k<0)
    return 0; /* done
	       * TODO: reset perm to identity?
	       */


  /* find smallest r>k with g[r]>g[k]
   * note g[k+1]>g[k] so r>=k+1 exists
   */
  r=k+1;
  while (r<n && abs(g[r])>abs(g[k])) r++;
  r--;

  /* swap g[k], g[r] */
  temp = g[r];
  g[r] = g[k];
  g[k] = temp;


  /* g[k+1],...,g[n-1] are descending,
   * reverse them to be in ascending order
   */

#if 0
  fprintf(outfile, "cA(%d,%d) ",k+2,n-1); /* debugging */
#endif

  for (i=k+1, j=n-1;
       i<j;
       i++, j--) {
#if 0
    fprintf(outfile, "swap(%d,%d) ", i, j); /* debugging */
#endif
    temp = g[i];
    g[i] = g[j];
    g[j] = temp;
  }

  return 1;
} 


/*****************************************************************************
 * cthash
 *****************************************************************************/

static int CTHASH_n = 0;

CTHASH *cthash_create(int reclen)
{
  CTHASH *h;
  h = e_malloc(sizeof(CTHASH),"cthash_create");
  h->reclen = reclen;
  h->root = (void *) NULL;

  CTHASH_n = reclen;

  return h;
}

int *cthash_alloc(CTHASH *h)
{
  int *p = e_malloc(h->reclen * sizeof(int), "cthash_alloc");
  return p;
}

int *cthash_insert(CTHASH *h, int *rec)
{
  int ***p;

  int **node = (int **) e_malloc(sizeof(int *), "cthash_insert.node");
  *node = rec;

  p = (int ***) tfind((const void *) node, &h->root, cthash_cmp);
  if (p == (int ***) NULL) {
    /* new record */
    /* insert it into tree */
    p = (int ***) tsearch((const void *) node, &h->root, cthash_cmp);
    if (p == (int ***) NULL) {
      exit(-1);
    }
    (***p) = 1;  /* init counter */
  } else {
    /* existing record */

    (***p)++;               /* bump count by 1 */
    free((void *) rec);    /* delete new version of record */
    free((void *) node);   /* and pointer to it */
  }
  return **p;
}

int cthash_cmp(const void *s1, const void *s2)
{
  const int *p1 = *((const int **)s1);
  const int *p2 = *((const int **)s2);

  int i;
  int d=0;
  int n = CTHASH_n;

  if (s1 == s2)
    return 0;

  for (i=1; i<n; i++) {
    d = p1[i] - p2[i];
    if (d != 0) break;
  }
  return d;
}


void ctype_calc(int num_genes, distmem_t *distmem, int *ctype_mem)
{
  int i;
  int cc_ind;
  int cc_left;

  int *cc         = distmem->cc;        /* connected component w/vertex */
  component_t *components = distmem->components; /* list of components */
  
  /* clear counts */
  memset(ctype_mem, 0, (num_genes+1)*sizeof(int));

  for (i=2*num_genes; i>=0; i -= 2) {
    cc_ind = cc[i];

    cc_left =
      (cc_ind < 0) ? i             /* adjacency, i is the leftmost vertex */
      : components[cc_ind].index;  /* normal component, get leftmost vert */

    /* bump count for # black edges in component */
    ctype_mem[cc_left/2]++;
  }

  /* sort counts into order */
  qsort((void *) ctype_mem, num_genes+1, sizeof(int), ctype_ptn_cmp);
}

int ctype_ptn_cmp(const void *a, const void *b)
{
  int x = *(const int *)a;
  int y = *(const int *)b;
  return y-x;
}


void ctprint(FILE *stream, int n, const int *p)
{
  int i;
  if (n>0) {
    fprintf(stream, " %10d", p[0]);
    for (i=1; i<n; i++) {
      fprintf(stream, " %3d", p[i]);
    }
  }
  fprintf(stream, "\n");
}

void cthash_print(CTHASH *cthash)
{
  fprintf(outfile,"          #   d  c2   h   f  c4  br  ku  ko  t1...\n");
  twalk(cthash->root, cthash_print_node);
}

void cthash_print_node(const void *node, VISIT order, int level)
{
  int *rec = **(int ***)node;
  if (order == preorder || order == leaf) {
    ctprint(outfile, CTHASH_n, rec);
  }
}

