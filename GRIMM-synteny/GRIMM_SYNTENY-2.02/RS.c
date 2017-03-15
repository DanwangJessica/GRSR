/* RS.c
 *    Robinson-Schensted algorithm for unsigned permutations
 *    (signs are ignored).
 *    Experimental for filtering out complex microrearrangements.
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

#include <stdio.h>
#include <stdlib.h>
#include "ckalloc.h"
#include "RS.h"

#define min2(x,y) ((x)<(y) ? (x) : (y))


/* allocate space for a Young tableau */
tableau_t *tableau_alloc(int n)
{
  int r;
  tableau_t *T = ckalloc(1,sizeof(tableau_t));
  T->nmax = n;
  T->nr = n;
  T->lens = (int *) ckalloc(n,sizeof(int));
  T->data = (int **) ckalloc(n,sizeof(int *));
  for (r=0; r<n; r++) {
    /* allocate space for row r;
     * the largest it could be is if there are r rows
     * with lengths (k+1,k+1,...,k+1,k,k,...,k)
     * where k=floor(n/(r+1))
     */
    T->data[r] = (int *) ckalloc(n/(r+1),sizeof(int));
  }
  return T;
}


void tableau_destroy(tableau_t *T)
{
  array2d_destroy(T->nr, (void **) T->data);

  /* TODO: check the frees below... */
#if 0
  free((void *) T->lens);
  free((void *) T);
#endif
}


void tableau_clear(tableau_t *T)
{
  int n = T->nr;
  int r;
  for (r=0; r<n; r++)
    T->lens[r] = 0;

  T->ncur = 0;
}

void tableau_print(FILE *output, tableau_t *T)
{
  int r, c;
  int n = T->nr;
  int ndig = 1;
  while (n>9) {
    n /= 10;
    ndig++;
  }

  if (T->lens[0] == 0) {
    fprintf(output,"(null tableau)\n");
  } else for (r=0; r<T->nr; r++) {
    if (T->lens[r] == 0)
      break;

    for (c=0; c < T->lens[r]; c++) {
      fprintf(output, " %*d", ndig, T->data[r][c]);
    }
    fprintf(output, "\n");
  }
}

void tableau_sh_print(FILE *output, tableau_t *T)
{
  int r;

  for (r=0; r<T->nr; r++) {
    if (T->lens[r] == 0 && r>0)
      break;
    fprintf(output, " %d", T->lens[r]);
  }
  fprintf(output,"\n");
}


/* insert k */
void RS_ins(tableau_t *P, tableau_t *Q, int k)
{
  int r, c;
  int cmax = P->lens[0];
  int knext;
  int ncur = P->ncur+1;

  for (r=0;
       r < ncur;
       r++) {
    /* TODO:
     * binary search?
     */
    /* determine where on row to place k */
    cmax = min2(cmax, P->lens[r]);
    for (c=cmax; c>0; c--) {
      /* could k be placed at row r, column c? */
      if (P->data[r][c-1] < k)
	break;
    }
    /* Place k at row r, column c.
     * Is that the end of the row or are we bumping something?
     */

    if (P->lens[r] == c) {
#if 0
      printf("At end of row %d, col %d, k=%d\n", r, c, k);
#endif

      /* end of row */
      P->lens[r]++;
      P->data[r][c] = k;
      P->ncur++;

      Q->lens[r]++;
      Q->data[r][c] = ncur;
      Q->ncur++;
      return;
    }

    /* bumping something */
    /* swap it with k */
    knext = P->data[r][c];
#if 0
    printf("Bumping in row %d, col %d, k=%d, knext=%d\n", r, c, k, knext);
#endif
    P->data[r][c] = k;
    k = knext;

    /* on next row, will not go past this column */
    cmax = c;
  }
  fprintf(stderr,"\nRS_ins: shouldn't reach here. r=%d, c=%d, knext=%d\n", r, c, knext);
  fprintf(stderr,"P=\n");
  tableau_print(stderr,P);
  fprintf(stderr,"Q=\n");
  tableau_print(stderr,Q);
}

/*
 * Perform ordinary Robinson-Schensted on signed permutation perm,
 * ignoring the signs.
 * Space for P, Q must be pre-allocated.
 */
void RSalg(int n, int *perm,
	   tableau_t *P, tableau_t *Q)
{
  int i;

  tableau_clear(P);
  tableau_clear(Q);
  for (i=0; i<n; i++)
    RS_ins(P,Q,abs(perm[i]));
}


/* compute length of the diagonal of a tableau */
int tableau_diag(tableau_t *T)
{
  int r;
  int nr = T->nr;
  for (r=1; r<=nr; r++) {
    if (T->lens[r-1] < r)
      return r-1;
  }

  /* can only get here if nr=0 or 1 */
  return nr;
}

/* determine if shape contains any particular cells
 *
 * input:
 *    T = tableau
 *    c = { i(1), j(1), i(2), j(2), i(3), j(3), ..., i(m), j(m) }
 *        0-based coordinates, matrix style
 *    n = 2m
 * return:
 *    1 if sh(T) contains any cells (i(k),j(k))
 *    0 o.w.
 */

int tableau_too_complex(tableau_t *T, int n, int *cells)
{
  int i,j,k;
  for (k=0; k<n; k+=2) {
    i = cells[k];
    j = cells[k+1];
    if (T->nr > i && T->lens[i] > j)
      return 1;
  }
  return 0;
}
