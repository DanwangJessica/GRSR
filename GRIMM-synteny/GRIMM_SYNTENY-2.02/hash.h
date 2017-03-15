/* hash.h
 *    Use a hash to convert chromosome names (strings) to numbers
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

#ifndef HASH_H
#define HASH_H

typedef struct {
  int len;         /* number of currently allocated entries */
  int maxlen;      /* number of entries allowed */
  char **strings;  /* strings in hash table */
  void *root;      /* tree to traverse strings */
} HTABLE;

HTABLE *hash_create(int maxlen);
void hash_destroy(HTABLE *h);
int hash_str2num(HTABLE *h, char *s);
int hash_cmp(const void *s1, const void *s2);
char *hash_num2str(HTABLE *h, int n);
void hash_sort(HTABLE *h, int (*compar)(const void *, const void *));

#endif /* HASH_H */
