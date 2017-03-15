/* grimm_anch.h
 *    GRIMM-Anchors algorithm: filter out orthologs with colliding coordinates.
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

/* GRIMM-Anchor graph components.
 */

#ifndef GRIMM_ANCH_H
#define GRIMM_ANCH_H

typedef struct {
  int n;                    /* nodes 0,...,n-1 */
  int nspecies;
  int *nodelist;            /* nodelist[i] = k
			     * means node i is externally called anchor k
			     */
  int **comp;               /* comp[j+1][i] =
			     * root of component with node i in species j,
			     * for j=-1,0,...,nspecies-1
			     */
} GANCTREE;

void grimm_anch(FILE *report_file, ANCTABLE *anctable, GSparams *gs_params);

#endif /* GRIMM_ANCH_H */
