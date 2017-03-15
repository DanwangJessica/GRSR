/* gsfile.c
 *    Open a file in a certain directory.
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
#include <string.h>
#include "gsfile.h"
#include "ckalloc.h"


/*****************************************************************************
 * Open a file in a directory
 *****************************************************************************/

FILE *fopen_dir(const char *dir, const char *fname, const char *mode)
{
  FILE *f;
  char *s;
  int len_dir, len_fname, len;

  /* make sure directory name and file name are non-null */
  if (dir == (char *) NULL) {
    fprintf(stderr,"fopen_dir: directory name is null!\n");
    exit(-1);
  }
  if (fname == (char *) NULL) {
    fprintf(stderr,"fopen_dir: filename is null!\n");
    exit(-1);
  }

  /* form string "dir/fname" */
  len_dir = strlen(dir);
  len_fname = strlen(fname);

  len = len_dir + len_fname + 2;
  s = (char *) ckalloc(1, len * sizeof(char));

  strcpy(s,dir);
  strcat(s,"/");
  strcat(s,fname);


  /* open file */
  f = fopen(s, mode);

  /* check for errors */
  if (f == (FILE *) NULL) {
    fprintf(stderr,"fopen_dir: could not open file %s\n", s);
    perror("");
    free((void *) s);
    exit(-1);
  }
  free((void *) s);
  return f;
}
