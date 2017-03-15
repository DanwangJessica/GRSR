/* matrixmisc.c
 *    Routines for 2-dim matrices.
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

/* Misc. routines for matrices
 */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "matrixmisc.h"
#include "e_malloc.h"
#include "mcread_input.h"

/* calc number of digits in n base 10 */
int num_digits(int n)
{
  int k=1;

  if (n<0) {
    k++;
    n = -n;
  }

  while (n>=10) {
    n /= 10;
    k++;
  }
  return k;
}


/* load matrix from a file */
void load_matrix(FILE *input,

		 /* return: */
		 int *nrows,
		 int *ncols,
		 int ***m
		 )
{
  int nr=0;
  int nc=0;
  int nc_cur = 0;  /* number of columns in current row */
  int **mat;
  int c;

  /* pass 1: validate input, calc dimensions
   * then: allocate memory  
   * pass 2: read input
   */

  /* pass 1 */
  while ((c=getc(input)) != EOF) {
    if (c == '#') {                          /* comment */
      discard_line(input);
      c = '\n';
    }

    if (c == '\n' || c == '\r' || c == EOF) { /* end of row */
      if (nc_cur == 0) {                      /* empty row */
	if (c == EOF) break;
	continue;
      }
      if (nc == 0) {
	nc = nc_cur;
      } else if (nc != nc_cur) {
	fprintf(stderr,
		"ERROR: load_matrix(): first %d rows have %d columns but next row has %d\n",
		nr, nc, nc_cur);
	exit(-1);
      }
      nr++;
      nc_cur = 0;
      continue;
    }

    /* skip white space other than newline */
    if (isspace(c)) continue;

    if (c == '-' || isdigit(c)) {
      /* skip over number */
      while ((c = getc(input)) != EOF && isdigit(c));
      if (c != EOF) ungetc(c,input);
      nc_cur++;
      continue;
    }

    /* unknown character */
    fflush(outfile); fclose(outfile);
    fprintf(stderr,
	    "ERROR: bad character %c in matrix file at position %d.\n",
	    c,
	    (int) ftell(input));
    exit(-1);
  }

  /* finished pass 1 */

  *nrows = nr;
  *ncols = nc;

  /* allocate memory */
  *m = mat = alloc_mat_int2d(nr,nc);

  /* pass 2 */
  nr = 0;
  rewind(input);
  while ((c=getc(input)) != EOF) {
    if (c == '#') {                          /* comment */
      discard_line(input);
      c = '\n';
    }

    if (c == '\n' || c == '\r' || c == EOF) { /* end of row */
      if (nc_cur == 0) {                      /* empty row */
	if (c == EOF) break;
	continue;
      }
      if (nc == 0) {
	nc = nc_cur;
      } else if (nc != nc_cur) {
	fprintf(stderr,
		"ERROR: load_matrix(): first %d rows have %d columns but next row has %d\n",
		nr, nc, nc_cur);
	exit(-1);
      }
      nr++;
      nc_cur = 0;
      continue;
    }

    /* skip white space other than newline */
    if (isspace(c)) continue;

    if (c == '-' || isdigit(c)) {
      ungetc(c,input);    /* put lead digit back in */
      fscanf(input, "%d", &mat[nr][nc_cur++]);
      continue;
    }
  }
}

int **alloc_mat_int2d(int nr, int nc)
{
  int i;
  int **mat = (int **) e_malloc(nr * sizeof(int *), "alloc_mat_int2d: mat");
  for (i=0; i<nr; i++) {
    mat[i] = (int *) e_malloc(nc * sizeof(int), "alloc_mat_int2d: mat[i]");
  }
  return mat;
}

void destroy_mat_int2d(int **mat, int nr)
{
  int i;
  for (i=nr-1; i>=0; i--) {
    free((void *) mat[i]);
  }
  free((void *) mat);
}


/* return: 1 = symmetric, 0 = not symmetric */
int check_sym_mat_int2d(int **mat, int nr, int nc)
{
  int i, j;
  if (nr != nc || nr<0 || nc<0) {
    return 0;
  }
  for (i=1; i<nr; i++) {
    for (j=0; j<i; j++) {
      if (mat[i][j] != mat[j][i]) {
	return 0;
      }
    }
  }

  return 1;
}		 
