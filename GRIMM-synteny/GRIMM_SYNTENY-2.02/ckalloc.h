/* ckalloc.h
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

#ifndef CKALLOC_H
#define CKALLOC_H

void *ckalloc(size_t nelem, size_t elsize);
int **int_array2d_create(size_t n1, size_t n2);
int **int_array2d_resize(int **p, size_t n1_old, size_t n2, size_t n1_new);
void int_array2d_cpy(int **p1, const int **p2, size_t n1, size_t n2);
void array2d_destroy(size_t n1, void **p);
void vec_012(int *p, int n);

#endif /* CKALLOC_H */
