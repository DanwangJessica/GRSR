/* RS.h
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

#ifndef RS_H
#define RS_H

typedef struct {
  int nmax;      /* max # elements */
  int ncur;      /* current # elements */
  int nr;        /* # rows allocated */
  int *lens;     /* lengths */
  int **data;    /* contents of tableau */
} tableau_t;


tableau_t *tableau_alloc(int n);
void tableau_destroy(tableau_t *T);
void tableau_clear(tableau_t *T);
void tableau_print(FILE *output, tableau_t *T);
void tableau_sh_print(FILE *output, tableau_t *T);
void RSalg(int n, int *perm,
	   tableau_t *P, tableau_t *Q);
int tableau_diag(tableau_t *T);
int tableau_too_complex(tableau_t *T, int n, int *cells);

#endif /* RS_H */
