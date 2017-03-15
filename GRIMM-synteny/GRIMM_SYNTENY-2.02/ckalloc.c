/* ckalloc.c
 *    Memory allocation and deallocation.
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

/* Last modified on Sun Sep 3, 2006, by Glenn Tesler
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ckalloc.h"

void *ckalloc(size_t nelem, size_t elsize)
{
  void *p = calloc(nelem,elsize);
  if (p != (void *) 0)
    return(p);
  fprintf(stderr,"Allocation error: nelem=%d, elsize=%d\n", nelem,elsize);
  exit(-1);
}

int **int_array2d_create(size_t n1, size_t n2)
{
  size_t i;
  int **p = (int **) ckalloc(n1, sizeof(int *));
  for (i=0; i<n1; i++) {
    p[i] = (int *) ckalloc(n2, sizeof(int));
  }
  return p;
}

int **int_array2d_resize(int **p, size_t n1_old, size_t n2, size_t n1_new)
{
  size_t i;
  int **p_new;
  if (n1_new < n1_old) {
    fprintf(stderr,"int_array2d_resize: bad dimensions, from %d x %d to %d x %d\n",
	    n1_old,n2, n1_new,n2);
  }
  p_new = (int **) realloc((void *) p, n1_new * sizeof(int *));
  if (p_new == (void *) NULL) {
    fprintf(stderr,"Resize allocation error: from %d x %d to %d x %d\n",
	    n1_old,n2, n1_new,n2);
  }
  for (i=n1_old; i<n1_new; i++) {
    p_new[i] = (int *) ckalloc(n2, sizeof(int));
  }
  return p_new;
}


/* copy data in 2d int array p2 into p1.
 * Space for p1 must already have been allocated.
 * Dimensions n1 x n2.
 */
void int_array2d_cpy(int **p1, const int **p2, size_t n1, size_t n2)
{
  size_t i;
  size_t nbytes = n2*sizeof(int);
  for (i=0; i<n1; i++) {
    memcpy((void *) p1[i], (const void *) p2[i], nbytes);
  }
}

void array2d_destroy(size_t n1, void **p)
{
  size_t i;

  if (p == (void **) NULL)
    return;

  for (i=0; i<n1; i++) {
    if (p[i] != (void *) NULL) {
      free((void *) p[i]);
    }
  }
  free((void *) p);
}

/*
 * set array to 0,1,2,3,...
 */
void vec_012(int *p, int n)
{
  int i;
  for (i=0; i<n; i++)
    p[i] = i;
}
