/* gscomp.h
 *    Form the components on the anchors and analyze them.
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


/* Anchor graph components.
 * We do not actually store the whole graph, we just store a spanning tree
 * to get its components.
 */

#ifndef GSCOMP_H
#define GSCOMP_H


#if 0
typedef struct {
  int *data;
  int len;                  /* current length */
  int maxlen;               /* maximum length */
} int_stack;
#endif

typedef struct {
  int n;                    /* nodes 0,...,n-1 */
  int *nodelist;            /* nodelist[i] = k
			     * means node i is externally called anchor k
			     */
  int *comp;                /* comp[i] = root of component with node i
			     * i=-1,0,1,2,...,n-1
			     * component -1 signifies deleted component
			     */
  int *comp_mem;            /* comp[i] is comp_mem[i+1] */

#if 0
  int *parent;              /* parent[i] = j
			     * when j is parent of i;
			     * parent[i] = i
			     * when i is root
			     */
  int *dtime;               /* dtime[i] = discovery time of node i */
  int *ftime;               /* ftime[i] = finishing time of node i */
  int dtime_cur, ftime_cur; /* current times */
  

  int *bridge;              /* bridge[i] =
			     *   1: if edge from i to parent is a bridge
			     *   0: o.w.
			     */

  int *degree;              /* degree[i] = degree of node i in spanning tree */

  int_stack *roots_l;       /* list of roots */
  int_stack *preorder_l;    /* list of nodes in preorder */
  int_stack *postorder_l;   /* list of nodes in postorder */
#endif
} ANCTREE;

ANCTREE *gstree_alloc(int n);
void gstree_destroy(ANCTREE *anctree);
void gstree(ANCTABLE *anctable, GSparams *gs_params,
	    ANCTREE **anctree_ret
	    );
int is_edge(ANCTABLE *anctable, ANCTREE *anctree, GSparams *gs_params,
	    int v, int w);

void comp_join(ANCTREE *anctree, int v, int w);
int comp_get(ANCTREE *anctree, int v);
void comp_optimize(ANCTREE *anctree);

void display_components_byanchor(FILE *output,
				 ANCTABLE *anctable, ANCTREE *anctree);


void list_anc_by_block(ANCTABLE *anctable, ANCTREE *anctree,
		       int ***blockanc_ret, int **node2blockpos_ret
		       );

void gs_filter_small(ANCTABLE *anctable,
		     GSparams *gs_params,
		     ANCTREE *anctree,
		     ANCTABLE *blocktable);

void gs_comp_coords(ANCTABLE *anctable, GSparams *gs_params, ANCTREE *anctree,
		    ANCTABLE **blocktable_ret,
		    int *max_nanc_ret
		    );

int gs_signs(FILE *output,
	     GSparams *gs_params,
	     ANCTABLE *anctable, ANCTREE *anctree,
	     ANCTABLE *blocktable,
	     int max_nanc,
	     int **blockanc, int *node2blockpos
	     );

int gs_filter_perms(FILE *output,
		    GSparams *gs_params,
		    ANCTABLE *anctable, ANCTREE *anctree,
		    ANCTABLE *blocktable,
		    int max_nanc,
		    int **blockanc, int *node2blockpos
		    );

void gs_micro(FILE *output,
	      ANCTABLE *anctable, ANCTREE *anctree,
	      ANCTABLE *blocktable,
	      int max_nanc,
	      int **blockanc, int *node2blockpos,
	      int RSon
	      );

void gs_macro(FILE *output,
	      GSparams *gs_params,
	      ANCTABLE *anctable, ANCTREE *anctree,
	      ANCTABLE *blocktable,
	      int **blockanc, int *node2blockpos
	      );


int gs_overlap(FILE *output,
	       GSparams *gs_params,
	       ANCTABLE *anctable, ANCTREE *anctree,
	       ANCTABLE *blocktable,
	       int **blockanc, int *node2blockpos
	       );

void gs_micro_mgr(GSparams *gs_params,
		  ANCTABLE *anctable, ANCTREE *anctree,
		  ANCTABLE *blocktable,
		  int max_nanc,
		  int **blockanc, int *node2blockpos
		  );

void gs_combine_consec_blocks(
			      ANCTABLE *anctable, ANCTREE *anctree,
			      ANCTABLE *blocktable,
			      int **blockanc, int *node2blockpos
			      );

void print_block(GSparams *gs_params, FILE *output,
		 ANCTABLE *anctable, int *anc);
void print_block_nocr(GSparams *gs_params, FILE *output,
		      ANCTABLE *anctable, int *anc);
void print_block2m(GSparams *gs_params, FILE *output,
		   ANCTABLE *anctable, int i);
void print_block_header(FILE *output, ANCTABLE *anctable);

/* E_ISEDGE: (v,w) is an edge
 * E_NOTEDGE: (v,w) is not an edge
 * E_TOOHIGH: (v,w) is not an edge AND if v was < w on input, then
 *            without any more tests there is no w'>w with (v,w') an edge
 */

#define E_ISEDGE 0
#define E_NOTEDGE 1
#define E_TOOHIGH 2

#endif /* GSCOMP_H */
