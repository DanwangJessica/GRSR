/* ext_function.h
 *    Template for new feature in GRIMM.  For testing only.
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

#ifndef EXT_FUNCTION_H
#define EXT_FUNCTION_H

void do_fn_external(int NUM_GENES, int NUM_CHROMOSOMES, int NUM_GENOMES,
		    int circular,
		    int unsigned_dist,
		    int verbose,
		    struct genome_struct *genome_list,
		    char *ext_args);

#endif /* EXT_FUNCTION_H */
