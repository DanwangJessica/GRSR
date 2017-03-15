/* time_stamp.c
 *    Print a time stamp and message.
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
#include <time.h>
#include "time_stamp.h"

void time_stamp(FILE *output, char *msg)
{
  time_t clock;
  char *s;

  if (output == (FILE *) NULL)
    output = stderr;

  /* get time */
  time(&clock);

  /* represent as string, and remove trailing \n */
  s = ctime(&clock);
  s[24] = '\0';

  fprintf(output, "%s", s);
  if (msg != (char *) NULL) {
    fprintf(output,": %s", msg);
  }
  fprintf(output,"\n");
  fflush(output);
}
