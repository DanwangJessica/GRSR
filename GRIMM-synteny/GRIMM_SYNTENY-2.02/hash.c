/* hash.c
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

#include <stdio.h>
#include <search.h>
#include <stdlib.h>
#include <string.h>
#include "hash.h"
#include "ckalloc.h"

/*****************************************************************************
 * hash table
 *****************************************************************************/

HTABLE *hash_create(int maxlen)
{
  HTABLE *h;
  h = ckalloc(1, sizeof(HTABLE));
  h->len = 0;
  h->maxlen = maxlen;
  h->strings = (char **) ckalloc(maxlen, sizeof(char *));
  h->root = (void *) NULL;
  return h;
}

void hash_destroy(HTABLE *h)
{
  int i;

  for (i=0; i < h->len; i++) {
    tdelete((const void *) &h->strings[i],
	    &h->root, hash_cmp);
    free((void *) h->strings[i]);
  }
  free((void *) h->strings);
  free((void *) h);
}

/* convert a string to a number in a hash table */
int hash_str2num(HTABLE *h, char *s)
{
  char *buf;
  char ***p;

  p = (char ***) tfind((const void *) &s, &h->root, hash_cmp);
  if (p != (char ***) 0) {
    return *p - &h->strings[0];
  }

  if (h->len >= h->maxlen) {
    fprintf(stderr,"Error: hash table overflow\n");
    exit(-1);
  }

  buf = (char *) ckalloc(1,sizeof(char) * (strlen(s)+1));
  strcpy(buf, s);
  h->strings[h->len] = buf;

  /* insert key into table */
  p = (char ***) tsearch((const void *) &h->strings[h->len], &h->root, hash_cmp);
  h->len++;
  if (p != (char ***) 0) {
    return *p - &h->strings[0];   /* TODO: long unsigned int */
  }

  fprintf(stderr,"Error inserting %s into table %p\n", s, (const void *) h);
  exit(-1);
}

int hash_cmp(const void *s1, const void *s2)
{
  return strcmp(*((const char **)s1),*((const char **)s2));
}

char *hash_num2str(HTABLE *h, int n)
{
  if (n >= 0
      && n < h->len)
    return h->strings[n];

  fprintf(stderr, "Hash %p entry # %d out of bounds\n", (const void *) h, n);
  return "";
}



void hash_sort(HTABLE *h, int (*compar)(const void *, const void *))
{
  /*
   * 1. empty the tree structure pointing at the strings
   * 2. sort strings by user specified sort
   * 3. reinsert strings into tree structure in new order
   */

  int i;

#if 0
  fprintf(stderr,"Order before:\n");
#endif

  /* 1. empty the tree structure pointing at the strings */
  for (i=0; i < h->len; i++) {
    tdelete((const void *) &h->strings[i],
	    &h->root, hash_cmp);
#if 0
    fprintf(stderr,"\t%d: %s\n", i, h->strings[i]);
#endif
  }

  /* 2. sort strings by user specified sort */
  qsort(h->strings, h->len, sizeof(char *), compar);

  /* 3. reinsert strings into tree structure in new order */
#if 0
  fprintf(stderr,"\nOrder after:\n");
#endif
  for (i=0; i < h->len; i++) {
    tsearch((const void *) &h->strings[i], &h->root, hash_cmp);
#if 0
    fprintf(stderr,"\t%d: %s\n", i, h->strings[i]);
#endif
  }
#if 0
  fprintf(stderr,"\n");
#endif
}
