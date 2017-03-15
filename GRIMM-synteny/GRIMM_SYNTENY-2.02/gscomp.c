/* gscomp.c
 *    Form the components on the anchors and analyze them.
 *
 * Copyright (C) 2001-2008 The Regents of the University of California
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

/* Last modified on Sun Mar 2, 2008, by Glenn Tesler
 */

/*
 * TODO:
 *  incrementally update the anchors-by-component table
 *  instead of recomputing it from scratch every time it changes
 */


#include <stdio.h>
#include <stddef.h>
/*#include <values.h>*/
#include <math.h>
#include <string.h>

#include "ckalloc.h"
#include "hash.h"
#include "anctable.h"
#include "gsy.h"
#include "gscomp.h"
#include "mcstructs.h"
#include "mcrdist.h"
#include "gs_grimm.h"
#include "gsfile.h"
#include "RS.h"

#define DB_RS
/* #undef DB_RS */

/* private */

int SNBC_cmp(const void *a, const void *b);

int calc_block_sign(int *blockanc_species, int base, int len,
		    int *node2blockpos, char *mem);


void gs_micro_1block(FILE *output,
		     int unsigned_flag,
		     int b_name,
		     int *block,   /* coords of this block */
		     int **blockanc,
		     int base,
		     int nanc_in_block,
		     int *node2blockpos,
		     ANCTABLE *anctable, ANCTREE *anctree,
		     distmem_t *distmem,
		     graphstats_t **statsmatrix,
		     struct genome_struct *genome_list);


void lmatrix_accum2(
				int field_offset,
				graphstats_t **graphstats_matrix,
				int num_genomes,
				struct genome_struct *genome_list);


void grimm_encode_mc(
		      int nblocks,
		      int max_chr,
		      int base_field,
		      struct genome_struct *g,
		      int *blockorder_species,
		      int *blockinds,
		      ANCTABLE *blocktable);

void grimm_encode_block(
		     int *block,   /* coords of this block */
		     int **blockanc,
		     int base,
		     int nanc_in_block,
		     int *node2blockpos,
		     ANCTABLE *anctable, ANCTREE *anctree,
		     struct genome_struct *genome_list);

void mgr_macro_mc(
		  FILE *output,
		      int nblocks,
		      int base_field,
		      struct genome_struct *g,
		      HTABLE *chrnames_table,
		      int *blockorder_species,
		      int *blockinds,
		      ANCTABLE *blocktable);



void gs_micro_mgr_1block(
			 FILE *output,
			 GSparams *gs_params,
			 int unsigned_flag,
			 int b_name,
			 int *block,   /* coords of this block */
			 int **blockanc,
			 int base,
			 int nanc_in_block,
			 int *node2blockpos,
			 ANCTABLE *anctable, ANCTREE *anctree,
			 ANCTABLE *blocktable,
			 int *strip_mem,
			 struct genome_struct *genome_list,
			 int **compressed_blocks
			 );

void compare_blocks_report(GSparams *gs_params,
			   int nspecies,
			   int nblocks,
			   int **compressed_blocks);
int cbr_cmp(const void *a0, const void *b0);
int cbr_cmp_stable(const void *a0, const void *b0);

int split_block(FILE *report_file,
		ANCTABLE *anctable, ANCTREE *anctree,
		ANCTABLE *blocktable,
		int *b_0,
		int j, /* species # */
		int c_A, int c_B,
		int A_startanc
		);

void sort_anc_by_chr(ANCTABLE *anctable,
		     int **nodelist_ret);
int SABC_cmp(const void *a0, const void *b0);
int calc_chr_window_end(ANCTABLE *anctable, ANCTREE *anctree, int v);
int calc_comp_window_end(ANCTREE *anctree, int v);



/* block list
 * format:
 *    blockinds = {v0,k0,l0,v1,k1,l1,v2,k2,l2,...}
 *    blockinds[i*BL_tot + BL_off_v] = v_i = component
 *    blockinds[i*BL_tot + BL_off_k] = k_i
 *                       = offset into sorted component/anchor table
 *    blockinds[i*BL_tot + BL_off_l] = l_i
 *                       = 1-based label used to name the block,
 *                         based on order in 1st genome
 */

#define BL_off_v 0
#define BL_off_k 1
#define BL_off_l 2
#define BL_tot   3


void list_blocks(
		 ANCTREE *anctree,
		 ANCTABLE *blocktable,
		 int **blockanc,

		 int *nblocks_ret,
		 int ***blockorders_ret,
		 int **blockinds_ret
		 );
int count_blocks(ANCTREE *anctree,
		 ANCTABLE *blocktable,
		 int **blockanc);



void print_anc_correlation(
			   FILE *output,			   
			   int *block,   /* coords of this block */
			   int **blockanc,
			   int base_anchor,
			   int nanc_in_block,
			   int *node2blockpos,
			   ANCTABLE *anctable, ANCTREE *anctree
		     );
double calc_anc_corr(
		     int *block,   /* coords of this block */
		     int **blockanc,
		     int base_anchor,
		     int nanc_in_block,
		     int *node2blockpos,
		     ANCTABLE *anctable, ANCTREE *anctree,
		     int species_x, int species_y,
		     int sign_b
		     );


void print_anc_support(
		       FILE *output,
		       int *block,   /* coords of this block */
		       int **blockanc,
		       int base_anchor,
		       int nanc_in_block,
		       int *node2blockpos,
		       ANCTABLE *anctable, ANCTREE *anctree
		       );


/*
 * Allocate space for a spanning tree with n nodes
 * node i is anchor # nodelist[i]
 * where these are sorted in an order suitable to finding edges
 *
 * For now, don't need whole tree, just the components, identified
 * by the lowest numbered anchor in the component
 */
ANCTREE *gstree_alloc(int n)
{
  ANCTREE *anctree = (ANCTREE *) ckalloc(1, sizeof(ANCTREE));

  anctree->n = n;
  anctree->nodelist = (int *) NULL;
  anctree->comp_mem = (int *) ckalloc(n+1, sizeof(int));
  anctree->comp     = &anctree->comp_mem[1];
  anctree->comp[-1] = -1;

  vec_012(anctree->comp, n);  /* set to {0,1,2,...,n-1} */

  return anctree;
}

void gstree_destroy(ANCTREE *anctree)
{
  free((void *) anctree->nodelist);
  free((void *) anctree->comp_mem);
  free((void *) anctree);
}


void gstree(ANCTABLE *anctable, GSparams *gs_params,
	    ANCTREE **anctree_ret
	    )
{
  int nanc = anctable->nanc;
  int edge_status;
  int v,w;


  int chr_window_end = -1;  /* highest node in current chromosome window */
  int comp_end = -1;        /* nodes (v,v+1,v+2,...,comp_end)
			     * are all in the same component.
			     * Higher ones may or may not be in the component.
			     */
  int has_gap;


  /* form a forest on all the anchors */
  ANCTREE *anctree = gstree_alloc(nanc);

  /* order them by chromosome window, then 1st genome */
  sort_anc_by_chr(anctable, &anctree->nodelist);
#if 0
  anctree->nodelist = anctable->orders[0];
#endif

  *anctree_ret = anctree;

  /* loop over potential edges */

  for (v=0; v<nanc; v++) {
    if (v > chr_window_end) {
      /* find new chromosome window end */
      chr_window_end = calc_chr_window_end(anctable, anctree, v);
    }
    if (v > comp_end) {
      comp_end = calc_comp_window_end(anctree, v);
    }

    has_gap = 0;   /* are there anchors > v intervening with block
		    * that are not in the block?
		    */
    for (w=comp_end+1; w<=chr_window_end; w++) {
      edge_status = is_edge(anctable,anctree,gs_params,v,w);
      /* will there possibly be ANY more edges (v,w') with w'>=w? */
      if (edge_status == E_TOOHIGH)
	break;  /* no */

      if (edge_status == E_NOTEDGE) {
	has_gap = 1;
	continue;
      }

      /* it's an edge! */
      comp_join(anctree,v,w);

      if (!has_gap)
	comp_end = w;
    }
  }

  comp_optimize(anctree);
}

/* find largest w s.t. nodes v and w have same chromosome window */
int calc_chr_window_end(ANCTABLE *anctable, ANCTREE *anctree, int v)
{
  int nanc = anctable->nanc;
  int nspecies = anctable->nspecies;
  int fields_per_species = anctable->fields_per_species;
  int **anclist = anctable->anclist;
  int *nodelist = anctree->nodelist;
  int *anc_v = anclist[nodelist[v]];
  int *anc_w;
  int w;
  int j;
  int col;

  for (w=v+1; w<nanc; w++) {
    anc_w = anclist[nodelist[w]];
    for (j=0; j<nspecies; j++) {
      col = j*fields_per_species + FIELD_base + FIELDR_chr;
      if (anc_v[col] != anc_w[col])
	return w-1; /* chromosome window has changed */
    }
  }

  /* window runs up to final node */
  return nanc-1;
}

/* find largest w s.t. nodes v,v+1,...,w are already in same component */
int calc_comp_window_end(ANCTREE *anctree, int v)
{
  int n = anctree->n;
  int w;
  int c_v = comp_get(anctree,v);

  for (w=v+1; w<n; w++) {
    if (c_v != comp_get(anctree,w)) {
      return w-1;  /* component has changed */
    }
  }

  /* window runs up to final node */
  return n-1;
}


/*
 * Determine if (v,w) is an edge
 *
 * Return:
 *  E_ISEDGE: it's an edge
 *  E_TOOHIGH: it's not an edge, and the distance in first org is
 *             so high that there will not be any more edges (v,w')
 *             with w' > w
 *  E_NOTEDGE: it's not an edge
 *
 *
 * TODO: circular chromosomes; unsigned anchors
 */
#undef is_edge_TEST_CHR
int is_edge(ANCTABLE *anctable, ANCTREE *anctree, GSparams *gs_params,
	    int v, int w)
{

  int nspecies = anctable->nspecies;
  int fields_per_species = anctable->fields_per_species;
  int *nodelist = anctree->nodelist;
  int base_field;
  int j;
  int gapthresh = gs_params->gapthresh;
  int *gapthresh_s = gs_params->gapthresh_s;  /* per-species */

  /* anchor coordinates */
  int start_v, len_v, end_v, sign_v;
  int start_w, len_w, end_w, sign_w;
#ifdef is_edge_TESTCHR
  int chr_v, chr_w;
#endif

  /* internal distance between anchors, ignoring their signs */
  int d_int = 0;
  int d_tot;

  /* adjustments for signs of anchors */
  int pen_v_pos=0, pen_v_neg=0, pen_w_pos=0, pen_w_neg=0;
  int v_sign, w_sign;

  /* data for the anchors */
  int *anc_v = anctable->anclist[nodelist[v]];
  int *anc_w = anctable->anclist[nodelist[w]];

  for (j=0; j<nspecies; j++) {
    base_field = j*fields_per_species + FIELD_base;

    /* chromosome, sign, start, end for each anchor */
#ifdef is_edge_TESTCHR
    chr_v = anc_v[base_field + FIELDR_chr];
    chr_w = anc_w[base_field + FIELDR_chr];

    /* different chromosomes, not an edge */
    if (chr_v != chr_w) {
      return (j>0) ? E_NOTEDGE : E_TOOHIGH;
    }
#endif

    start_v = anc_v[base_field + FIELDR_start];
    len_v = anc_v[base_field + FIELDR_len] - 1;
    end_v = start_v + len_v;
    sign_v = anc_v[base_field + FIELDR_sign];

    start_w = anc_w[base_field + FIELDR_start];
    len_w = anc_w[base_field + FIELDR_len] - 1;
    end_w = start_w + len_w;
    sign_w = anc_w[base_field + FIELDR_sign];

    /* ASSUME LINEAR CHROMOSOMES AND ANCHORS DO NOT OVERLAP */
    /* unsigned distance between their closest ends */
    d_int += min2(abs(start_w-end_v),abs(start_v-end_w));

    /* additional distance to +/- terminals
     * ASSUMES linear chromosomes and anchors do not overlap
     */
    if (end_v < start_w) {
      if (sign_v > 0) {
	pen_v_neg += len_v;
      } else {
	pen_v_pos += len_v;
      }
      if (sign_w > 0) {
	pen_w_pos += len_w;
      } else {
	pen_w_neg += len_w;
      }
    } else {
      if (sign_w > 0) {
	pen_w_neg += len_w;
      } else {
	pen_w_pos += len_w;
      }
      if (sign_v > 0) {
	pen_v_pos += len_v;
      } else {
	pen_v_neg += len_v;
      }
    }

    /* TODO: adjustments if circular */

    /* fail if partial distance already too high
     *
     * TODO: could use a better bound than gapthresh in determining E_TOOHIGH
     * but this one is valid, just not optimal
     */

    d_tot = d_int + min2(pen_v_pos,pen_v_neg) + min2(pen_w_pos,pen_w_neg);
    if (d_tot >= gapthresh) {
      return (j>0) ? E_NOTEDGE : E_TOOHIGH;
    }
  }

  /* 2nd test:
   * per-species threshold
   */

  /* total is small enough that each component had to be small enough too */
  if (d_tot < gs_params->gapthresh_s_min)
    return E_ISEDGE;


  /* + if pen_v_pos smaller, - if pen_v_neg smaller, 0 if tie */
  v_sign = pen_v_neg - pen_v_pos;
  w_sign = pen_w_neg - pen_w_pos;

  /* secondary test for it being an edge, species-by-species */
  for (j=0; j<nspecies; j++) {
    base_field = j*fields_per_species + FIELD_base;

    start_v = anc_v[base_field + FIELDR_start];
    len_v = anc_v[base_field + FIELDR_len] - 1;
    end_v = start_v + len_v;
    sign_v = anc_v[base_field + FIELDR_sign];

    start_w = anc_w[base_field + FIELDR_start];
    len_w = anc_w[base_field + FIELDR_len] - 1;
    end_w = start_w + len_w;
    sign_w = anc_w[base_field + FIELDR_sign];

    /* ASSUME LINEAR CHROMOSOMES AND ANCHORS DO NOT OVERLAP */
    /* unsigned distance between their closest ends */
    d_int = min2(abs(start_w-end_v),abs(start_v-end_w));

    /* additional distance to +/- terminals
     * ASSUMES linear chromosomes and anchors do not overlap
     */
    if (end_v < start_w) {
      if (sign_v > 0) {
	pen_v_neg = len_v;
	pen_v_pos = 0;
      } else {
	pen_v_pos = len_v;
	pen_v_neg = 0;
      }
      if (sign_w > 0) {
	pen_w_pos = len_w;
	pen_w_neg = 0;
      } else {
	pen_w_neg = len_w;
	pen_w_pos = 0;
      }
    } else {
      if (sign_w > 0) {
	pen_w_neg = len_w;
	pen_w_pos = 0;
      } else {
	pen_w_pos = len_w;
	pen_w_neg = 0;
      }
      if (sign_v > 0) {
	pen_v_pos = len_v;
	pen_v_neg = 0;
      } else {
	pen_v_neg = len_v;
	pen_v_pos = 0;
      }
    }
    if (v_sign < 0)
      pen_v_pos = pen_v_neg;
    if (w_sign < 0)
      pen_w_pos = pen_w_neg;
    d_tot = d_int + pen_v_pos + pen_w_pos;

    if (d_tot >= gapthresh_s[j]) {
      if (j==0 && d_int >= gapthresh_s[j])
	return E_TOOHIGH;
      return E_NOTEDGE;
    }
  }

  /* success, edge */
  return E_ISEDGE;
}


/* join the components with nodes v and w */
void comp_join(ANCTREE *anctree, int v, int w)
{
  int root_v = comp_get(anctree,v);
  int root_w = comp_get(anctree,w);
  int root = min2(root_v,root_w);
  anctree->comp[root_v] = anctree->comp[root_w] = root;
}

/* delete component with node v
 * by joining it to component -1
 */
void comp_del(ANCTREE *anctree, int v)
{
  anctree->comp[comp_get(anctree,v)] = -1;
}

/*
 * comp_get: get the component # with node v
 * updates intermediate nodes for efficiency
 */
int comp_get(ANCTREE *anctree, int v)
{
  int v0 = v;
  int w;
  int *comp = anctree->comp;

  /* follow v, comp[v], comp[comp[v]], ... till it stabilizes
   */
  while (v != (w=comp[v])) {
    v = w;
  }

  /* then update all the pointers to point to the new component root
   */
  while (v0 != (w=comp[v0])) {
    comp[v0] = v;
    v0 = w;
  }

  return v;
}

/*
 * optimize all component pointers to point to component root
 */
void comp_optimize(ANCTREE *anctree)
{
  int n = anctree->n;
  int *comp = anctree->comp;

  int i;

  /* note comp[-1] = -1 already */
  for (i=0; i<n; i++)
    comp[i] = comp[comp[i]];
}

#if 0
/* display all anchors with their component numbers */
void display_components_byanchor(FILE *output,
				 ANCTABLE *anctable, ANCTREE *anctree)
{
  int n = anctree->n;
  int nspecies = anctable->nspecies;
  int fields_per_species = anctable->fields_per_species;
  int *nodelist = anctree->nodelist;
  int **anclist = anctable->anclist, *anc;
  int *comp = anctree->comp;
  HTABLE **chrnames = anctable->chrnames;
  int base_field;
  int i,j;
  int chr_num, start, len, sign, comp_no;
  int v;

  for (v=0; v<n; v++) {
    i = nodelist[v];      /* anchor number */
    comp_no = comp[v];    /* component number */
    anc = anclist[i];     /* data for the anchor */

    fprintf(output, "%d %d ", i, comp_no);
    fprintf(output, "%d.0", anc[FIELD_score]);
    for (j=0; j<nspecies; j++) {
      base_field = j*fields_per_species + FIELD_base;
      chr_num = anc[base_field + FIELDR_chr];
      start = anc[base_field + FIELDR_start];
      len = anc[base_field + FIELDR_len];
      sign = anc[base_field + FIELDR_sign];

      fprintf(output, " %s %d %d %s",
	      hash_num2str(chrnames[j],chr_num),start,len,
	      sign > 0 ? "+" : "-");
    }
    fprintf(output, "\n");
  }
}
#endif

/*
 * Convert components to blocks
 */

void gs_comp_coords(ANCTABLE *anctable, GSparams *gs_params, ANCTREE *anctree,

		    /* table of blocks */
		    ANCTABLE **blocktable_ret,

		    /* max # anchors per block */
		    int *max_nanc_ret
		    )
{
  int n = anctree->n;
  int nspecies = anctable->nspecies;
  int fields_per_anchor = anctable->fields_per_anchor;
  int fields_per_species = anctable->fields_per_species;
  int **anclist;
  int *anc, *block;
  ANCTABLE *blocktable;
  int *nodelist = anctree->nodelist;
  int *comp = anctree->comp;
  int **anclist_b;
  int **anclist_b_sup;
  int v, c_v;
  int j;
  int base_field;

  int start_a,end_a,len_a;
  int start_b,end_b,len_b;

  int max_nanc = 1;

  int pass;

  /* will be compressed later */
  blocktable = (ANCTABLE *) ckalloc(1, sizeof(ANCTABLE));
  *blocktable_ret = blocktable;
  blocktable->nanc = n;
  blocktable->nspecies = nspecies;
  blocktable->fields_per_species = FIELDR_n_per_species;
  blocktable->fields_per_anchor = anctable->fields_per_anchor;

  /* will be compressed later;
   * unused block #s have NULL pointers
   */
  blocktable->chrnames = anctable->chrnames;
  blocktable->max_chr = anctable->max_chr;

  /*
   * Form a table of all components, in similar format to anchors
   * blocktable->anclist[v] != NULL:
   *         info for block whose root is node v
   * blocktable->anclist[v] == NULL:
   *         no block has root v
   *         OR the block for root v was deleted
   *
   * may have two passes:
   *  anclist_b
   *  anclist2_b
   *
   */


  for (pass=0; pass<=gs_params->perm_metric_use; pass++) {
    if (pass == 0) {
      anclist = anctable->anclist;
      anclist_b =
	blocktable->anclist = (int **) ckalloc(n, sizeof(int *));
      anclist_b_sup =
	blocktable->anclist3 = (int **) ckalloc(n, sizeof(int *));
    } else {
      anclist = anctable->anclist2;
      anclist_b =
	blocktable->anclist2 = (int **) ckalloc(n, sizeof(int *));
      anclist_b_sup = (int **) NULL;
    }

    for (v=0; v<n; v++) {
      /* anchor corresponding to this node */
      anc = anclist[nodelist[v]];
      c_v = comp[v];       /* get component # with node v */
      if (c_v == v) {
	/* new component! */
	anclist_b[c_v] = (int *) ckalloc(fields_per_anchor, sizeof(int));

	if (anclist_b_sup != (int **) NULL) {
	  anclist_b_sup[c_v] = (int *) ckalloc(fields_per_anchor, sizeof(int));
	}

	/* initialize fields to agree with the root anchor;
	 * many of the fields will be changed
	 */
	memcpy((void *) anclist_b[c_v],
	       (const void *) anc,
	       fields_per_anchor * sizeof(int));
	anclist_b[c_v][FIELD_nanc] = 1;
      } else if (c_v >= 0) {
	/* determine coordinate extents in each species */
	block = anclist_b[c_v];
	block[FIELD_nanc]++;

	max_nanc = max2(max_nanc, block[FIELD_nanc]);

	for (j=0; j<nspecies; j++) {
	  base_field = j*fields_per_species + FIELD_base;

	  /* anchor coords for this organism */
	  start_a = anc[base_field + FIELDR_start];
	  len_a = anc[base_field + FIELDR_len];
	  end_a = start_a + len_a;

	  if (anclist_b_sup != (int **) NULL) {
	    /* accumulate support */
	    anclist_b_sup[c_v][base_field + FIELDR_len] += len_a;
	  }

	  /* cumulative block coords for this organism */
	  start_b = block[base_field + FIELDR_start];
	  len_b = block[base_field + FIELDR_len];
	  end_b = start_b + len_b;

	  /* combine */
	  start_b = min2(start_a,start_b);
	  end_b = max2(end_a,end_b);
	  block[base_field + FIELDR_start] = start_b;
	  block[base_field + FIELDR_len] = end_b - start_b;
	}
      }
    }
  }
  *max_nanc_ret = max_nanc;
}


/*
 * filter out small blocks
 */

void gs_filter_small(ANCTABLE *anctable,
		     GSparams *gs_params,
		     ANCTREE *anctree,
		     ANCTABLE *blocktable)
{
  int *min_block_size_s = gs_params->min_block_size_s;
  int *min_block_size2_s = gs_params->min_block_size2_s;
  int *min_block_sup_s  = gs_params->min_block_sup_s;
  int min_block_nanc = gs_params->min_block_nanc;
  int perm_metric = gs_params->perm_metric_use;

  int n = anctree->n;
  int *comp = anctree->comp;
  int nspecies = anctable->nspecies;
  int fields_per_species = anctable->fields_per_species;
  int **anclist_b = blocktable->anclist;
  int **anclist2_b = blocktable->anclist2;
  int **anclist_b_sup = blocktable->anclist3;
  int *block, *block_sup, *block_nuc;
  int v, del, j;
  int del_any = 0;
  int base_field;

  if (gs_params->min_block_size_s_max <= 0
      && (!perm_metric || gs_params->min_block_size2_s_max <= 0)
      && min_block_nanc <= 1
      )
    return;

  /* loop over roots of components */
  for (v=0; v<n; v++) {
    if (v != comp[v])
      continue;  /* not root of component */

    /* block coordinates */
    block = anclist_b[v];
    block_sup = anclist_b_sup[v];
    if (perm_metric)
      block_nuc = anclist2_b[v];
    del = 0;
    for (j=0; j<nspecies; j++) {
      base_field = j*fields_per_species + FIELD_base;
	if (block[base_field + FIELDR_len] < min_block_size_s[j]
	    ||
	    block_sup[base_field + FIELDR_len] < min_block_sup_s[j]
	    ||
	    block[FIELD_nanc] < min_block_nanc
	    || (perm_metric &&
		block_nuc[base_field + FIELDR_len] < min_block_size2_s[j])
	    ) {

#if 0
	  /* DEBUGGING */
	  if (block[base_field + FIELDR_len] >= min_block_size_s[j]) {
	    printf("small block %d: species %d, len=%d (vs %d), sup=%d (vs %d), nanc=%d (vs %d)\n",
		   v, j+1,
		   block[base_field + FIELDR_len], min_block_size_s[j],
		   block_sup[base_field+FIELDR_len], min_block_sup_s[j],
		   block[FIELD_nanc], min_block_nanc);
	    print_block(gs_params, stdout, blocktable, block);
	  }
#endif

	  /* block too small, mark for deletion */
	  del = 1;
	  break;
	}
    }

    if (del) {
      /* delete block */
      comp_del(anctree, v);
      del_any = 1;
    }
  }
  if (del_any)
    comp_optimize(anctree);
}

/*
 * Compute signs of each block
 * return value:
 *   0: good
 *   1: some blocks had bad signs or bad correlation coefficients and were
 *      marked for deletion
 */

int gs_signs(FILE *output,
	     GSparams *gs_params,
	     ANCTABLE *anctable, ANCTREE *anctree,
	     ANCTABLE *blocktable,
	     int max_nanc,
	     int **blockanc, int *node2blockpos
	     )
{
  int nspecies = anctable->nspecies;
  int fields_per_species = anctable->fields_per_species;
  int **anclist_b = blocktable->anclist;
  int perm_metric = gs_params->perm_metric_use;
  int **anclist_b2 = perm_metric ? blocktable->anclist2 : (int **) NULL;
  int i, k, c_v, nanc_in_block, m;
  int j;
  int sign;
#if 0
  int dbg_i;
#endif

  /* a bit array would work too */
  char *signmem = (char *) ckalloc(max_nanc, sizeof(char));

  int FIELD0_sign = 0*fields_per_species + FIELDR_sign + FIELD_base;
  int *bo0;
  int b_name;

  int nblocks;
  int **blockorders;
  int *blockinds;
  int *block;
  int *block_b2;
  int err;

  int any_err=0;
  double r, abs_r;

  int sign_r;



  /* get list of all blocks, sorted by order */
  list_blocks(anctree, blocktable, blockanc,
	      &nblocks, &blockorders, &blockinds);

  /* examine all anchors by component */
  bo0 = blockorders[0];
  for (i=0; i<nblocks; i++) {
    m = bo0[i];
    b_name = blockinds[m*BL_tot + BL_off_l];
    c_v = blockinds[m*BL_tot + BL_off_v];
    k = blockinds[m*BL_tot + BL_off_k];
    block = anclist_b[c_v];                      /* block coords */
    block_b2 = perm_metric ? anclist_b2[c_v] : (int *) NULL; /* 2nd coord sys */
    nanc_in_block = block[FIELD_nanc];           /* # anchors in block */

    /* if block has just 1 anchor, the sign is already determined */
    if (nanc_in_block < 2)
      continue;

    /* signs in first species are +1,
     * signs in other species are relative to that
     */
    block[FIELD0_sign] = +1;
    if (perm_metric) {
      block_b2[FIELD0_sign] = +1;
    }

    err = 0;

    for (j=1; j<nspecies; j++) {
      sign = calc_block_sign(blockanc[j], k, nanc_in_block, node2blockpos,
			     signmem);

      /* correlation coefficient too */
      r = calc_anc_corr(block, blockanc, k,
			nanc_in_block, node2blockpos,
			anctable, anctree,
			0,j,1);

      sign_r = (r > 0) ? 1 : (r<0) ? -1 : 0;

      if (sign != 0 && sign != sign_r) {
	fprintf(output,
		"Block %d, species %d: r=%f but combinatorial sign=%d\n",
		b_name, j+1, r, sign);
	err=1;
      }

      abs_r = (r>0) ? r : -r;
      if (abs_r < gs_params->min_corr) {
	if (gs_params->sign_err_del) {
	  fprintf(output,
		  "Block %d, species %d: r=%f, discarding block\n",
		  b_name, j+1, r);
	} else {
	  fprintf(output,
		  "Block %d, species %d: r=%f, low but keeping block\n",
		  b_name, j+1, r);
	}
	err=1;
      }
      
      if (sign == 0) {
	fprintf(output,
		"Block %d: combinatorial sign unresolveable in species %d\n",
		b_name, j+1);
#if 0
	for (dbg_i=0; dbg_i<nanc_in_block; dbg_i++) {
	  fprintf(output,
		  "%d ", 1+node2blockpos[blockanc[j][k+dbg_i]]);
	}
	fprintf(output,
		"\n");
#endif
	sign = (r>=0) ? 1 : -1;
	err = 1;
      }

      /*
       * TODO: compare correlation coeff & signs in all other pairs of species
       */

      block[j*fields_per_species + FIELDR_sign + FIELD_base] = sign;
      if (perm_metric) {
	block_b2[j*fields_per_species + FIELDR_sign + FIELD_base] = sign;
      }

      if (err) {
	if (gs_params->sign_err_del) {
	  /* delete block */
	  comp_del(anctree,c_v);
	  any_err = 1;
	  fprintf(output,"Deleting block\n");
	}

	print_block2m(gs_params, output, blocktable, c_v);
      }
    }
  }

  if (any_err)
    comp_optimize(anctree);

  free((void *) blockinds);
  array2d_destroy(nspecies, (void **) blockorders);

  free((void *) signmem);
  return any_err;
}

/*
 * calculate sign of a block relative to reference order
 * input:
 *   blockanc[base],blockanc[base+1],...,blockanc[base+len-1]
 *   are the node #s
 *
 *   node2blockpos: array s.t.
 *    node2blockpos[blockanc[base]],...,node2blockpos[blockanc[base+len-1]]
 *   is a permutation of 0,1,2,...,base+len-1
 *
 *   mem: a character array of length at least len
 *        which we will use for workspace
 *        TODO: use bit array?
 *
 * return value:
 *   +1,-1,0; 0 means can't determine sign
 */
int calc_block_sign(int *blockanc_species, int base, int len,
		    int *node2blockpos, char *mem)
{
  int size_prefix = 0;
  int size_suffix = 0;
  int len1 = len-1;

  int i, ic, f_i;

  for (i=0; i<len; i++) {
    mem[i] = 0;
  }

  /* Define f(i) = node2blockpos[blockanc_species[base + i]]
   * If f(0),f(1),...,f(k-1) is a perm of 0,...,k-1
   * for any k < len, then return sign +1;
   * If f(0),f(1),...,f(k-1) is a perm of len-k,...,len-1
   * for any k < len, then return sign -1.
   */
  for (i=0; i<len1; i++) {
    /* forwards strip:
     * has i already been encountered?
     * has f_i already been encountered?
     */
    if (mem[i] & 1) {
      size_prefix++;
    } else {
      mem[i] |= 1;
    }
    f_i = node2blockpos[blockanc_species[base+i]];
    if (mem[f_i] & 1) {
      size_prefix++;
    } else {
      mem[f_i] |= 1;
    }
    if (size_prefix == i+1) {
      /* f_0,...,f_i is a permutation of 0,...,i */
      return 1;
    }

    /* reverse strip:
     * has len-1-i already been encountered?
     * has f_i already been encountered?
     */
    ic = len1-i;
    if (mem[ic] & 2) {
      size_suffix++;
    } else {
      mem[ic] |= 2;
    }
    if (mem[f_i] & 2) {
      size_suffix++;
    } else {
      mem[f_i] |= 2;
    }
    if (size_suffix == i+1) {
      /* f_0,...,f_i is a perm of len1,len1-1,...,len1-i
       */
      return -1;
    }
  }

  /* unable to determine sign */
  return 0;
}

/*****************************************************************************
 * Compute Robinson-Schensted tableau for the perms given by each block,
 * and filter out those that are deemed to be too complex
 *****************************************************************************/

int gs_filter_perms(FILE *output,
		    GSparams *gs_params,
		    ANCTABLE *anctable, ANCTREE *anctree,
		    ANCTABLE *blocktable,
		    int max_nanc,
		    int **blockanc, int *node2blockpos
		    )
{
  int nspecies = anctable->nspecies;
  int **anclist_b = blocktable->anclist;
  int i, k, c_v, nanc_in_block, m;
  int j;

  struct genome_struct *genome_list;

  int *bo0;

  int nblocks;
  int **blockorders;
  int *blockinds;
  int *block;
  int err;

  int any_err=0;

  tableau_t *P, *Q;

  /* allocate maximum amount of space needed by GRIMM for any block;
   * the space can be reused for all blocks
   */
  gs_grimm_alloc(0,          /* unichromosomal */
		 max_nanc,   /* max # "genes" */
		 nspecies,   /* # species */

		 /* various GRIMM structures */
		 (distmem_t **) NULL, (graphstats_t ***) NULL, &genome_list);


  /* allocate space for Robinson-Schensted */
  P = tableau_alloc(max_nanc);
  Q = tableau_alloc(max_nanc);


  /* get list of all blocks, sorted by order */
  list_blocks(anctree, blocktable, blockanc,
	      &nblocks, &blockorders, &blockinds);

  /* examine all anchors by component */
  bo0 = blockorders[0];
  for (i=0; i<nblocks; i++) {
    m = bo0[i];
    c_v = blockinds[m*BL_tot + BL_off_v];
    k = blockinds[m*BL_tot + BL_off_k];
    block = anclist_b[c_v];                      /* block coords */
    nanc_in_block = block[FIELD_nanc];           /* # anchors in block */

    err = 0;

    /* Form GRIMM encoding of the anchors in the block for all species */
    grimm_encode_block(block, blockanc, k, nanc_in_block, node2blockpos,
		       anctable, anctree, genome_list);

    for (j=1; j<nspecies; j++) {
      RSalg(nanc_in_block, &genome_list[j].genes[0], P, Q);
      /* test complexity of shape of P */
      if (tableau_too_complex(P,
			      gs_params->RS_complex_n,
			      gs_params->RS_complex_cells)) {
	fprintf(output,"Permutation too complex, species 0 vs %d; RS shape ",
		j);
	tableau_sh_print(output, P);
	err = 1;
      }
    }

    if (err) {
      comp_del(anctree,c_v);
      any_err = 1;
      fprintf(output,"Deleting block\n");
      print_block2m(gs_params, output, blocktable, c_v);
    }
  }

  /* free memory */

  free((void *) blockinds);
  array2d_destroy(nspecies, (void **) blockorders);

  tableau_destroy(Q);
  tableau_destroy(P);
  gs_grimm_destroy(0,          /* unichromosomal */
		   max_nanc,   /* max # "genes" */
		   nspecies,   /* # species */

		   /* various GRIMM structures */
		   (distmem_t *) NULL, (graphstats_t **) NULL,
		   genome_list);


  if (any_err)
    comp_optimize(anctree);

  return any_err;
}


/*****************************************************************************
 * Combine consecutive blocks into one
 *****************************************************************************/

void gs_combine_consec_blocks(
			      ANCTABLE *anctable, ANCTREE *anctree,
			      ANCTABLE *blocktable,
			      int **blockanc, int *node2blockpos
			      )
{
  int nspecies = anctable->nspecies;
  int fields_per_species = anctable->fields_per_species;
  int i,j;

  int nblocks, nblocks1;
  int **blockorders;
  int *blockinds;
  int *strip_mem;

  int col_sign, col_chr;
  int m1,v1,sign1,*block1;
  int i2,m2,v2,sign2,*block2;
  int b_name1,b_name2;

  /* get list of all blocks, sorted by order */
  list_blocks(anctree, blocktable, blockanc,
	      &nblocks, &blockorders, &blockinds);

  strip_mem = (int *) ckalloc(nblocks, sizeof(int));


  /* compute strips
   * Let the blocks be numbered i=0,1,2,...,n-1 according to order in
   * species 0. 
   * 1. set strip_mem[i]=0 for all i, 0<=i<n.
   * 2. go through all blocks in all species j>0;
   *    if (i,i+1) is a breakpoint in any species j>0 relative to species 0,
   *    or if (i,i+1) is a chromosome break in species j,
   *    set strip_mem[i]=1
   */

#if 0
  /* Form GRIMM encoding of the anchors in the block for all species */
  grimm_encode_block(block, blockanc, base, nanc_in_block, node2blockpos,
		     anctable, anctree, genome_list);
#endif


  /* determine strips */
  nblocks1 = nblocks-1;
  for (j=0; j<nspecies; j++) {
    /* find strip boundaries in species j vs. species 0,
     * and chromosome boundaries in species j
     */
    col_sign = FIELD_base + j*fields_per_species + FIELDR_sign;
    col_chr = FIELD_base + j*fields_per_species + FIELDR_chr;
    for (i=0; i<nblocks; i++) {
      m1 = blockorders[j][i];
      b_name1 = blockinds[m1*BL_tot + BL_off_l];
      v1 = blockinds[m1*BL_tot + BL_off_v];
      block1 = blocktable->anclist[v1];
      sign1 = block1[col_sign];

      i2 = i+sign1;

      if (i2<0 || i2>=nblocks) {
	strip_mem[m1] = 1;
	continue;
      }

      m2 = blockorders[j][i2];
      b_name2 = blockinds[m2*BL_tot + BL_off_l];
      v2 = blockinds[m2*BL_tot + BL_off_v];
      block2 = blocktable->anclist[v2];
      sign2 = block2[col_sign];

      if (sign1 == sign2
	  && block1[col_chr] == block2[col_chr]
	  && b_name1 + 1 == b_name2
	  ) {
	continue;
      }
      strip_mem[m1] = 1;
    }
  }

  /* merge consec blocks into single blocks */
  for (i=0; i<nblocks1; i++) {
    if (!strip_mem[i]) {
      /* merge i-th and (i+1)-st block into one */
      comp_join(anctree,
		blockinds[i*BL_tot + BL_off_v],
		blockinds[(i+1)*BL_tot + BL_off_v]);
    }
  }

  /* recompute root pointers */
  comp_optimize(anctree);

  /* free block list */
  free((void *) blockinds);
  array2d_destroy(nspecies, (void **) blockorders);

  free((void *) strip_mem);
}


/*****************************************************************************
 * Microrearrangements
 *****************************************************************************/

#define print_lmatrix_field_asrow(display,var) \
    print_lmatrix_offset_asrow(output, display, \
                    offsetof(graphstats_t, var), \
                    statsmatrix, \
                    nspecies, \
                    genome_list)

#define print_umatrix_field_asrow(display,var) \
    print_umatrix_offset_asrow(output, display, \
                    offsetof(graphstats_t, var), \
                    statsmatrix, \
                    nspecies, \
                    genome_list)


/* print a matrix of the values of a particular field in the distance
 * structures */
void print_lmatrix_offset_asrow(
				FILE *output,
				char *matrix_name,
				int field_offset,
				graphstats_t **graphstats_matrix,
				int num_genomes,
				struct genome_struct *genome_list) {
  int i, j;     /* loop over matrix entries */

  fprintf(output,"  %-18s ",matrix_name);

  for (i=1 ; i<num_genomes ; i++) {
    fprintf(output,"| ");
    for (j=0 ; j<i ; j++) {
      fprintf(output,"%4d ",
              *(int *)((char *)(&graphstats_matrix[i][j]) + field_offset));
    }
  }
  fprintf(output,"\n");
}
void print_umatrix_offset_asrow(
				FILE *output,
				char *matrix_name,
				int field_offset,
				graphstats_t **graphstats_matrix,
				int num_genomes,
				struct genome_struct *genome_list) {
  int i, j;     /* loop over matrix entries */

  fprintf(output,"  %-18s ",matrix_name);

  for (i=1 ; i<num_genomes ; i++) {
    fprintf(output,"| ");
    for (j=0 ; j<i ; j++) {
      fprintf(output,"%4d ",
              *(int *)((char *)(&graphstats_matrix[j][i]) + field_offset));
    }
  }
  fprintf(output,"\n");
}


#define print_smatrix_field(name,var) \
    print_smatrix_offset(output, name, \
                    offsetof(graphstats_t, var), \
                    statsmatrix, \
                    nspecies, \
                    genome_list)
/* print a matrix of the values of a particular field in the distance
 * structures */
void print_smatrix_offset(
			  FILE *output,
			  char *matrix_name,
			  int field_offset,
			  graphstats_t **graphstats_matrix,
			  int num_genomes,
			  struct genome_struct *genome_list) {
  int i, j;

  fprintf(output,"%s Matrix:\n", matrix_name);
  for (i=0 ; i<num_genomes ; i++) {
    for (j=0 ; j<num_genomes ; j++) {
      fprintf(output,
	      " %4d",
	     *(int *)((char *)(&graphstats_matrix[i][j]) + field_offset));
    }
    fprintf(output,"\n");
  }
  fprintf(output,"\n");
}




#define lmatrix_accum(var) \
    lmatrix_accum2( \
                    offsetof(graphstats_t, var), \
                    statsmatrix, \
                    nspecies, \
                    genome_list)

/* accumulate running totals in the upper half of the matrix
 */
void lmatrix_accum2(
				int field_offset,
				graphstats_t **graphstats_matrix,
				int num_genomes,
				struct genome_struct *genome_list) {
  int i, j;     /* loop over matrix entries */

  for (i=1 ; i<num_genomes ; i++) {
    for (j=0 ; j<i ; j++) {
      *(int *)((char *)(&graphstats_matrix[j][i]) + field_offset)
	+= *(int *)((char *)(&graphstats_matrix[i][j]) + field_offset);
    }
  }
}


/* add something for cumulative statistics? */
void gs_micro_1block(
		     FILE *output,
		     int unsigned_flag,
		     int b_name,
		     int *block,   /* coords of this block */
		     int **blockanc,
		     int base,
		     int nanc_in_block,
		     int *node2blockpos,
		     ANCTABLE *anctable, ANCTREE *anctree,
		     distmem_t *distmem,
		     graphstats_t **statsmatrix,
		     struct genome_struct *genome_list)
{
  int nspecies = anctable->nspecies;
  int i,j;

  /* encode block for GRIMM */
  grimm_encode_block(block, blockanc, base, nanc_in_block, node2blockpos,
		     anctable, anctree, genome_list);


  /*
   * The distance parameters are symmetric and the diagonal is predictable,
   * so only compute below-diagonal portion of matrix
   */
  for (i=1; i<nspecies; i++) {
    for (j=0; j<i; j++) {
      /* TODO: unsigned */
      mcdist_noncircular(&genome_list[i],
			 &genome_list[j],
			 nanc_in_block,
			 0, /* # chromosomes: 0 means unichromosomal */
			 distmem,
			 &statsmatrix[i][j]);
      /* KLUDGE: store reuses in n field */
      statsmatrix[i][j].n =
	2*statsmatrix[i][j].d - statsmatrix[i][j].br;
    }
  }

  /*
   * display results
   */

  fprintf(output,"\nblock %d: %d anchors\n", b_name, nanc_in_block);

  /* report on support (#/% nucleotides) of anchor coordinates */
  print_anc_support(output,
		    block, blockanc, base, nanc_in_block, node2blockpos,
		    anctable, anctree);


  print_lmatrix_field_asrow("distance", d);
  print_lmatrix_field_asrow("# breakpoints", br);
  print_lmatrix_field_asrow("# long cycles", c4);
  print_lmatrix_field_asrow("# hurdles", h);
  print_lmatrix_field_asrow("# fortresses", f);

  /* KLUDGE */
  print_lmatrix_field_asrow("# reuses", n);

  /* report on correlation coefficient of anchor coordinates */
  print_anc_correlation(output,
			block, blockanc, base, nanc_in_block, node2blockpos,
			anctable, anctree);


  /* running totals */
  lmatrix_accum(d);
  lmatrix_accum(br);
  lmatrix_accum(c4);
  lmatrix_accum(h);
  lmatrix_accum(f);
}


void grimm_encode_block(
		     int *block,   /* coords of this block */
		     int **blockanc,
		     int base,
		     int nanc_in_block,
		     int *node2blockpos,
		     ANCTABLE *anctable, ANCTREE *anctree,
		     struct genome_struct *genome_list)
{
  int fields_per_species = anctable->fields_per_species;
  int *nodelist = anctree->nodelist;
  int nspecies = anctable->nspecies;
  int i,j,v,f_i,sign;
  int *anc_v;
  int sign_b;

  /* encode genes orders for GRIMM */
  for (j=0; j<nspecies; j++) {
    /* block sign in species j */
    sign_b = block[j*fields_per_species + FIELDR_sign + FIELD_base];

    for (i=0; i<nanc_in_block; i++) {
      v = blockanc[j][base+i];        /* node # */

      /* Get order of this anchor w/in block in species j rel. to species 0.
       * It's 0-based, need 1-based.
       */
      f_i = node2blockpos[v] + 1;

      /* anchor */
      anc_v = anctable->anclist[nodelist[v]];

      /* anchor sign */
      sign = anc_v[j*fields_per_species + FIELDR_sign + FIELD_base];

      /* store in block order, reverse block if block is - for this genome */
      if (sign_b > 0) 
	genome_list[j].genes[i] = sign>0 ? f_i : -f_i;
      else
	genome_list[j].genes[nanc_in_block-1-i] = sign<0 ? f_i : -f_i;
    }

#if 0
    /* debugging */
    printf("species %d:", j+1);
    for (i=0; i<nanc_in_block; i++) {
      printf(" %d", genome_list[j].genes[i]);
    }
    printf("\n");
#endif
  }
}


/*
 * Analyze and print report on microrearrangements, block-by-block
 */

void gs_micro(FILE *output,
	      ANCTABLE *anctable, ANCTREE *anctree,
	      ANCTABLE *blocktable,
	      int max_nanc,
	      int **blockanc, int *node2blockpos,
	      int RSon
	      )
{
  int nspecies = anctable->nspecies;
  int **anclist_b = blocktable->anclist;
  int *bo0;
  int *block;
  int k, c_v, nanc_in_block;
  int b_name;
  int i;
  int m;

  distmem_t *distmem;
  graphstats_t **statsmatrix;
  struct genome_struct *genome_list;

  int tot_nanc = 0;

  int nblocks;
  int **blockorders;
  int *blockinds;

  int unsigned_flag = 0;

#ifdef DB_RS
  tableau_t *P = (tableau_t *) NULL, *Q = (tableau_t *) NULL;
  int s;
  int nRS;
#endif

  /* allocate maximum amount of space needed by GRIMM for any block;
   * the space can be reused for all blocks
   */
  gs_grimm_alloc(0,          /* unichromosomal */
		 max_nanc,   /* max # "genes" */
		 nspecies,   /* # species */

		 /* various GRIMM structures */
		 &distmem, &statsmatrix, &genome_list);


#ifdef DB_RS
  /* allocate space for Robinson-Schensted */
  if (RSon) {
    P = tableau_alloc(max_nanc);
    Q = tableau_alloc(max_nanc);
  }
#endif

  /* get list of all blocks, sorted by order */
  list_blocks(anctree, blocktable, blockanc,
	      &nblocks, &blockorders, &blockinds);

  /* examine all anchors by component */
  bo0 = blockorders[0];
  for (i=0; i<nblocks; i++) {
    m = bo0[i];
    b_name = blockinds[m*BL_tot + BL_off_l];
    c_v = blockinds[m*BL_tot + BL_off_v];
    k = blockinds[m*BL_tot + BL_off_k];
    block = anclist_b[c_v];                      /* block coords */
    nanc_in_block = block[FIELD_nanc];           /* # anchors in block */

    /* analyze block */
    tot_nanc += nanc_in_block;
    gs_micro_1block(output,
		    unsigned_flag,
		    b_name,
		    block,
		    blockanc, k, nanc_in_block, node2blockpos,
		    anctable, anctree,
		    distmem, statsmatrix, genome_list);

#ifdef DB_RS
    /* TODO: debugging, move elsewhere later */
    if (RSon) {
      for (s=1; s<nspecies; s++) {

      nRS = nanc_in_block;
      RSalg(nRS, &genome_list[s].genes[0], P, Q);
      fprintf(output, "RS species 0 vs. %d\n", s);

#if 0
      fprintf(output, "perm: ");
      for (m=0; m<nRS; m++)
	fprintf(output, "%d ", genome_list[s].genes[m]);
      fprintf(output, "\n");
#endif

      fprintf(output, "RS species 0 vs %d shape: ", s);
      tableau_sh_print(output, P);

#if 0
      fprintf(output, "P=\n");
      tableau_print(output, P);

      fprintf(output, "Q=\n");
      tableau_print(output, Q);
#endif
      }
    }
#endif /* DB_RS */
  }

  /*
   * print totals
   */

  fprintf(output,"\nTotals: %d anchors in %d blocks\n", tot_nanc, nblocks);
  print_umatrix_field_asrow("distance", d);
  print_umatrix_field_asrow("# breakpoints", br);
  print_umatrix_field_asrow("# long cycles", c4);
  print_umatrix_field_asrow("# hurdles", h);
  print_umatrix_field_asrow("# fortresses", f);
  print_umatrix_field_asrow("# reuses", n);  /* KLUDGE */
  fprintf(output,"\n");

  /* free memory */

  free((void *) blockinds);
  array2d_destroy(nspecies, (void **) blockorders);

#ifdef DB_RS
  if (RSon) {
    tableau_destroy(Q);
    tableau_destroy(P);
  }
#endif

  gs_grimm_destroy(0,          /* unichromosomal */
		   max_nanc,   /* max # "genes" */
		   nspecies,   /* # species */

		   /* various GRIMM structures */
		   distmem, statsmatrix, genome_list);

}


void print_anc_correlation(
			   FILE *output,
			   int *block,   /* coords of this block */
			   int **blockanc,
			   int base_anchor,
			   int nanc_in_block,
			   int *node2blockpos,
			   ANCTABLE *anctable, ANCTREE *anctree
			   )
{
  int i,j;
  int nspecies = anctable->nspecies;
  double r;

  fprintf(output,"  %-18s ","Anchor correlation");
  for (i=1; i<nspecies; i++) {
    fprintf(output,"| ");
    for (j=0 ; j<i ; j++) {
      r = calc_anc_corr(block, blockanc, base_anchor,
			nanc_in_block, node2blockpos,
			anctable, anctree,
			i, j, 0);

      fprintf(output,"%4.2f ",r);
    }
  }
  fprintf(output,"\n");
}


/*
 * calculate correlation coefficient between coordinates of
 * anchor in species_x and species_y
 *
 * sign_b=0: lookup block signs
 *      !=0: use this sign as relative sign of block in the 2 species
 */
double calc_anc_corr(
		     int *block,   /* coords of this block */
		     int **blockanc,
		     int base_anchor,
		     int nanc_in_block,
		     int *node2blockpos,
		     ANCTABLE *anctable, ANCTREE *anctree,
		     int species_x, int species_y,
		     int sign_b
		     )
{
  double Sxy = 0.0;
  double Sx  = 0.0;
  double Sx2 = 0.0;
  double Sy  = 0.0;
  double Sy2 = 0.0;
  double r, m;

  int fields_per_species = anctable->fields_per_species;
  int *nodelist = anctree->nodelist;

  int base_x = species_x * fields_per_species + FIELD_base;
  int base_y = species_y * fields_per_species + FIELD_base;

  double start_x, len_x, end_x, start_y, len_y, end_y;
  double mid_x, mid_y;
  double len_m, len_m1;   /* len_m=min(len_x, len_y); len_m1 = len_m - 1 */ 
  int sign_x, sign_y;

#if 0
  int n = nanc_in_block;
#endif
  int n = 0; /* # points in correlation, = sum of anchor lengths */
  int *anc;

  int i, v;

  int *b_0 = blockanc[0];

  /* adjustments to improve precision */
  double min_x,min_y,max_x,max_y,adj_x,adj_y;


  /* shouldn't happen */
  if (nanc_in_block == 0)
    return 0;

  /* reduce overflows by determining an amount to downshift all coordinates */
  /* block 0 */
    v = b_0[base_anchor];        /* node # for anchor 0 */

    /* anchor */
    anc = anctable->anclist[nodelist[v]];

    /* coords */
    start_x = anc[base_x + FIELDR_start];
    len_x = anc[base_x + FIELDR_len];
    mid_x = start_x + (len_x - 1)/2;
    min_x = max_x = mid_x;

    start_y = anc[base_y + FIELDR_start];
    len_y = anc[base_y + FIELDR_len];
    mid_y = start_y + (len_y - 1)/2;
    min_y = max_y = mid_y;

  /* blocks > 0 */
  for (i=1; i<nanc_in_block; i++) {
    v = b_0[base_anchor+i];        /* node # */

    /* anchor */
    anc = anctable->anclist[nodelist[v]];

    /* coords */
    start_x = anc[base_x + FIELDR_start];
    len_x = anc[base_x + FIELDR_len];
    mid_x = start_x + (len_x - 1)/2;
    if (min_x > mid_x) { min_x = mid_x; }
    if (max_x < mid_x) { max_x = mid_x; }

    start_y = anc[base_y + FIELDR_start];
    len_y = anc[base_y + FIELDR_len];
    mid_y = start_y + (len_y - 1)/2;
    if (min_y > mid_y) { min_y = mid_y; }
    if (max_y < mid_y) { max_y = mid_y; }
  }
  adj_x = (min_x + max_x)/2;
  adj_y = (min_y + max_y)/2;
  
  

  for (i=0; i<nanc_in_block; i++) {
    v = b_0[base_anchor+i];        /* node # */

    /* anchor */
    anc = anctable->anclist[nodelist[v]];

    /* coords */
    start_x = anc[base_x + FIELDR_start] - adj_x;
    len_x = anc[base_x + FIELDR_len];
    end_x = start_x + len_x - 1;
    mid_x = (start_x + end_x)/2;
    sign_x = anc[base_x + FIELDR_sign];

    start_y = anc[base_y + FIELDR_start] - adj_y;
    len_y = anc[base_y + FIELDR_len];
    end_y = start_y + len_y - 1;
    mid_y = (start_y + end_y)/2;
    sign_y = anc[base_y + FIELDR_sign];

    /* 
     * Correlation using whole anchors.
     * Assumes length in two species is similar.
     *
     * Generate len_m=min(len_x,len_y) points.
     * For len_x<len_y, use m=(len_y-1)/(len_x-1):
     *    if sign_x==sign_y:
     *       (start_x,start_y),
     *       (start_x+1,start_y+1*m),
     *       (start_x+2,start_y+2*m), ...
     *       (start_x+len_x-1,start_y+(len_x-1)*m) = (end_x,end_y)
     *   
     *    if sign_x != sign_y:
     *       (start_x,end_y),
     *       (start_x+1,end_y-1*m),
     *       (start_x+2,end_y-2*m), ...
     *
     * For len_y>len_x, use m=(len_x-1)/(len_y-1):
     *    if sign_x==sign_y:
     *       (start_x,start_y),
     *       (start_x+1*m,start_y+1),
     *       (start_x+2*m,start_y+2), ...
     *       (start_x+(len_m-1)*m,start_y+(len_m-1)) = (end_x,end_y)
     *   
     *    if sign_x != sign_y:
     *       (end_x,start_y),
     *       (end_x-1*m,start_y+1), ...
     *
     * Formulas for Sx2 = sum of x^2 values, etc., computed in maple
     */

    if (len_x < len_y) {
      len_m = len_x;
      m = (len_y - 1)/(len_x - 1);
    } else {
      len_m = len_y;
      m = (len_x - 1)/(len_y - 1);
    }
    len_m1 = len_m - 1;

    if (len_m == 0) continue;

    if (len_m == 1) {
      /* contribution from single point (mid_x, mid_y) */
      Sx += mid_x;
      Sy += mid_y;
      Sxy += mid_x*mid_y;
      Sx2 += mid_x*mid_x;
      Sy2 += mid_y*mid_y;
      n++;
      continue;
    }


    m = (double)(len_y - 1) / (double)(len_x - 1);

    if (sign_x != sign_y) {
      /* line segment with opposite slope */
      m = -m;
    }


    Sx  += len_m * mid_x;
    Sy  += len_m * mid_y;
    Sxy += len_m * (mid_x*mid_y + m*(len_m-1)*(len_m+1)/12);
    /* actually, m*(len_m-1)*(len_m+1) = s * (len_max - 1) * (len_min + 1)
     * where len_max=max(len_x,len_y), len_min=min(len_x,len_y)
     * and s = sign_x * sign_y
     */

    Sx2 +=
      len_m / (3*len_m1) *(
			   (4*len_m - 2)* mid_x*mid_x
			   - (len_m+1) * start_x * end_x
			   );
    Sy2 +=
      len_m / (3*len_m1) *(
			   (4*len_m - 2)* mid_y*mid_y
			   - (len_m+1) * start_y * end_y
			   );
    n += len_x;
  }

  r = (Sxy - Sx*Sy/n) / sqrt((Sx2 - Sx*Sx/n) * (Sy2 - Sy*Sy/n));


  /* combined block sign */
  if (sign_b == 0)
    sign_b = block[base_x + FIELDR_sign] * block[base_y + FIELDR_sign];

  /* adjust for overall sign of block */
  r *= sign_b;

  return r;
}


/*
 * report on support (# nucleotides in alignments vs. in whole block)
 */ 

void print_anc_support(
		       FILE *output,
		       int *block,   /* coords of this block */
		       int **blockanc,
		       int base_anchor,
		       int nanc_in_block,
		       int *node2blockpos,
		       ANCTABLE *anctable, ANCTREE *anctree
		       )
{
  int i,j,p;
  int nspecies = anctable->nspecies;
  int fields_per_species = anctable->fields_per_species;
  int col_len;
  int block_len;
  int support;
  int v;
  int *anc;
  int *nodelist = anctree->nodelist;
  int *b_0 = blockanc[0];

  for (p=0; p<3; p++) {
    fprintf(output,"  %-18s ",
	    (p == 0) ? "Block lengths" :
	    (p == 1) ? "Support (nuc)" :
	               "Support (%)");

    for (j=0; j<nspecies; j++) {
      support = 0;
      col_len = j*fields_per_species + FIELD_base + FIELDR_len;
      block_len = block[col_len];

      if (p>0) {
	for (i=0; i<nanc_in_block; i++) {
	  v = b_0[base_anchor+i];        /* node # */

	  /* anchor */
	  anc = anctable->anclist[nodelist[v]];

	  support += anc[col_len];
	}
      }

      switch (p) {
      case 0:
	fprintf(output, " %9d", block_len);
	break;
      case 1:
	fprintf(output, " %9d", support);
	break;
      case 2:
	fprintf(output, " %8.5f%%", 100*(double)support / (double)block_len);
	break;
      }
    }
    fprintf(output,"\n");
  }
}

/*****************************************************************************
 * Sorts
 *****************************************************************************/

/*
 * SNBC = sort_nodes_by_comp
 * sort anchors for each species by (component,chromosome,start,len,node #)
 */

static int **SNBC_anclist;
static int *SNBC_nodelist;
static int SNBC_base;
static int *SNBC_comp;

/* list all anchors by block, sorted for each species
 * return:
 *   blockanc[j][0], blockanc[j][1], ..., blockanc[j][nanc-1]
 *   are the node #s sorted in order for species j
 *   (the anchor number of node v is anctree->nodelist[v])
 *
 *   node2anclabel[v] = k
 *     means in the block containing node v, the nodes are
 *     assigned numbers 0,1,2,...,nanc_in_block-1
 *     and k is one of these numbers
 *
 */
void list_anc_by_block(ANCTABLE *anctable, ANCTREE *anctree,
		       int ***blockanc_ret, int **node2blockpos_ret
		       )
{
  int **blockanc, *b_0, *b_j;
  int *node2blockpos;
  int j,v;
  int nspecies = anctable->nspecies;
  int fields_per_species = anctable->fields_per_species;
  int n = anctree->n;

  int cur_comp, cur_pos_in_comp, c_v, k;
  int *comp = anctree->comp;

  blockanc = int_array2d_create(nspecies, n);
  *blockanc_ret = blockanc;

  /* node lists */
  b_0 = blockanc[0];
  for (j=0; j<nspecies; j++) {
    b_j = blockanc[j];

    if (j == 0) {
      /* species 0: sort from scratch */
      vec_012(b_j, n);  /* set to {0,1,2,...,n-1} */
    } else {
      /* other species: order will be close to the order in species 0 */
      memcpy((void *) b_j, (const void *) b_0, n*sizeof(int));
    }

    /* sort it in order for species j by
     * (block #, chromosome #, start, len, node #)
     */
    SNBC_anclist = anctable->anclist;
    SNBC_nodelist = anctree->nodelist;
    SNBC_comp = comp;
    SNBC_base = j*fields_per_species + FIELD_base;
    qsort((void *) b_j, n, sizeof(int),
	  SNBC_cmp);
  }

  /*
   * compute position # within each block
   */

  node2blockpos = (int *) ckalloc(n, sizeof(int));
  *node2blockpos_ret = node2blockpos;

  cur_comp = -2; /* root of current component */
  cur_pos_in_comp = 0;

  /* examine anchors in order for species 0,
   * sorted first by block, second by coordinates
   */
  for (k=0; k<n; k++) {
    v = b_0[k];              /* node # */
    c_v = comp[v];           /* component # */
    if (c_v != cur_comp) {
      /* starting a new component, restart labels */
      cur_comp = c_v;
      cur_pos_in_comp = 0;
    }
    node2blockpos[v] = cur_pos_in_comp++;
  }
}

int SNBC_cmp(const void *a, const void *b)
{
  int *anc1 =
    SNBC_anclist[SNBC_nodelist[*(int *)a]];
  int *anc2 =
    SNBC_anclist[SNBC_nodelist[*(int *)b]];

  /* first compare component numbers */
  int d =
    SNBC_comp[*(int *)a] - SNBC_comp[*(int *)b];

  if (d == 0) {
    /* break tie with chromosomes */
    d = anc1[SNBC_base + FIELDR_chr] - anc2[SNBC_base + FIELDR_chr];

    if (d == 0) {
      /* break tie with start coordinates */
      d = anc1[SNBC_base + FIELDR_start] - anc2[SNBC_base + FIELDR_start];

      if (d == 0) {
	/* break tie with length */
	d = anc1[SNBC_base + FIELDR_len] - anc2[SNBC_base + FIELDR_len];

	if (d == 0) {
	  /* break tie with node number */
	  d = *(int *)a - *(int *)b;
	}
      }
    }
  }


  return d;
}

/*
 * sort blocks for each species by block coordinates
 */
int SB_cmp(const void *a, const void *b)
{
  int *anc1 =
    SNBC_anclist[SNBC_nodelist[(*(int *)a) * BL_tot + BL_off_v]];
  int *anc2 =
    SNBC_anclist[SNBC_nodelist[(*(int *)b) * BL_tot + BL_off_v]];

  /* first compare chromosomes */
  int d = anc1[SNBC_base + FIELDR_chr] - anc2[SNBC_base + FIELDR_chr];

  if (d == 0) {
    /* break tie with start coordinates */
    d = anc1[SNBC_base + FIELDR_start] - anc2[SNBC_base + FIELDR_start];

    if (d == 0) {
      /* break tie with length */
      d = anc1[SNBC_base + FIELDR_len] - anc2[SNBC_base + FIELDR_len];

      if (d == 0) {
	/* break tie with node number */
	d = *(int *)a - *(int *)b;
      }
    }
  }

  return d;
}

/*
 * SABC = sort anchors by chromosome window
 *
 * list all anchors in order by
 *   (chr1,chr2,chr3,...,start1,len1,start2,len2,....,node #)
 *
 * return nodelist = {a1,a2,...}
 *                 = anchor #'s in order as above
 *
 */

static int *SABC_nodelist;
static ANCTABLE *SABC_anctable;

void sort_anc_by_chr(ANCTABLE *anctable,
		     int **nodelist_ret)
{
  int nanc = anctable->nanc;
  int *nodelist = (int *) ckalloc(nanc, sizeof(int));
  vec_012(nodelist, nanc);  /* set to {0,1,2,...} */
  SABC_nodelist = nodelist;
  SABC_anctable = anctable;
  qsort((void *) nodelist, nanc, sizeof(int), SABC_cmp);
  *nodelist_ret = nodelist;
}

/*
 * compare anchors by chromosome # in each species;
 * if same in all species, then by position in 1st species, 2nd species, etc.
 */
int SABC_cmp(const void *a0, const void *b0)
{
  int nspecies = SABC_anctable->nspecies;
  int fields_per_species = SABC_anctable->fields_per_species;
  int **anclist = SABC_anctable->anclist;
  int *anc_a = anclist[*(int *)a0];
  int *anc_b = anclist[*(int *)b0];
  int j;
  int col;
  int d;

  /* compare chromosome window */
  for (j=0; j<nspecies; j++) {
    col = j*fields_per_species + FIELD_base + FIELDR_chr;
    d = anc_a[col] - anc_b[col];
    if (d != 0) return d;
  }

  /* compare coordinates, should resolve on species 0 */
  for (j=0; j<nspecies; j++) {
    /* compare by start coordinate */
    col = j*fields_per_species + FIELD_base + FIELDR_start;
    d = anc_a[col] - anc_b[col];
    if (d != 0) return d;

    /* compare by end coordinate */
    col += FIELDR_len - FIELDR_start;
    d = anc_a[col] - anc_b[col];
    if (d != 0) return d;
  }

  /* Shouldn't reach here with valid anchors.
   * If do, just make it a stable sort.
   */
  d = *(int *)a0 - *(int *)b0;
  return d;
}





/*****************************************************************************
 * 
 *****************************************************************************/


/* list all blocks, sorted for each species
 * return:
 *   nblocks = # blocks
 *
 *   blockorders[j][0], ..., blockorders[j][nblocks-1]
 *   are block indices sorted in order for species j
 *
 *   Let r = blockorders[j][i]
 *       v = blockinds[r*BL_tot + BL_OFF_v]
 *       k = blockinds[r*BL_tot + BL_OFF_k]
 *
 *   Let v = blockinds[blockorders[j][i]]
 *   then the block info is at blocktable->anclist[v]
 *   and anchor list is
 *     blockanc[j][k], blockanc[j][k+1], ..., blockanc[j][k+nanc_in_block-1]
 */

void list_blocks(
		 ANCTREE *anctree,
		 ANCTABLE *blocktable,
		 int **blockanc,

		 int *nblocks_ret,
		 int ***blockorders_ret,
		 int **blockinds_ret
		 )
{
  int n = blocktable->nanc;
  int nspecies = blocktable->nspecies;
  int fields_per_species = blocktable->fields_per_species;
  int **anclist_b = blocktable->anclist;
  int *b_0 = blockanc[0];
  int *block;
  int k, v, c_v, nanc_in_block, k_next;
  int j;
  int *comp = anctree->comp;

  int bno, nblocks;
  int *blockinds;
  int **blockorders;
  int *b_j;
  int i,m;

  nblocks = count_blocks(anctree, blocktable, blockanc);
  *nblocks_ret = nblocks;

  /* list of roots of all the blocks */
  blockinds = (int *) ckalloc(nblocks*BL_tot, sizeof(int));

  /* skip deleted nodes
   * TODO: store number of deleted nodes somewhere?
   */
  bno = 0;

  for (k=0; k<n; k = k_next) {
    v = b_0[k];                                  /* node # */
    c_v = comp[v];                               /* component # */
    if (c_v == -1) {
      k_next = k+1;                              /* component deleted */
      continue;
    }
      

    block = anclist_b[c_v];                      /* block coords */
    nanc_in_block = block[FIELD_nanc];           /* # anchors in block */
    k_next = k + nanc_in_block;
    blockinds[bno*BL_tot + BL_off_v] = c_v;
    blockinds[bno*BL_tot + BL_off_k] = k;
    bno++;
  }
  *blockinds_ret = blockinds;

  /* get block order in each species */
  blockorders = int_array2d_create(nspecies, nblocks);
  for (j=0; j<nspecies; j++) {
    b_j = blockorders[j];
    vec_012(b_j,nblocks); /* set to {0,1,2,...} */

    /* sort it in order for species j by
     * (chromosome #, start, len, node #)
     */
    SNBC_anclist = anclist_b;
    SNBC_nodelist = blockinds;
    SNBC_comp = (int *) NULL; /* not used */
    SNBC_base = j*fields_per_species + FIELD_base;
    qsort((void *) b_j, nblocks, sizeof(int),
	  SB_cmp);

  }
  *blockorders_ret = blockorders;

  /* assign block labels in order by 1st genome */
  b_j = blockorders[0];
  for (i=0; i<nblocks; i++) {
    m = b_j[i];
    blockinds[m*BL_tot + BL_off_l] = i+1;
  }

}

/*
 * count # remaining blocks
 */

int count_blocks(ANCTREE *anctree,
		 ANCTABLE *blocktable,
		 int **blockanc)
{
  int n = blocktable->nanc;
  int **anclist_b = blocktable->anclist;
  int *b_0 = blockanc[0];
  int *comp = anctree->comp;

  int k, k_next;
  int v, c_v;
  int nanc_in_block;
  int *block;
  int nblocks = 0;

  for (k=0; k<n; k = k_next) {
    v = b_0[k];                                  /* node # */
    c_v = comp[v];                               /* component # */
    if (c_v == -1) {
      k_next = k+1;                              /* component deleted */
      continue;
    }
    block = anclist_b[c_v];                      /* block coords */
    nanc_in_block = block[FIELD_nanc];           /* # anchors in block */
    k_next = k + nanc_in_block;
    nblocks++;
  }
  return nblocks;
}

/*****************************************************************************
 * Output anchors in block for MGR to analyze
 * Also output block table
 *****************************************************************************/

void gs_micro_mgr(GSparams *gs_params,
		  ANCTABLE *anctable, ANCTREE *anctree,
		  ANCTABLE *blocktable,
		  int max_nanc,
		  int **blockanc, int *node2blockpos
		  )
{
  int nspecies = anctable->nspecies;

  int **anclist_b =
    gs_params->perm_metric_use
    ? blocktable->anclist2
    : blocktable->anclist;

  int *bo0;
  int *block;
  int k, c_v, nanc_in_block;
  int i,m;
  int b_name;

  int **compressed_blocks;

  int *strip_mem;
  struct genome_struct *genome_list;

  int nblocks;
  int **blockorders;
  int *blockinds;

  int unsigned_flag = 0;



  FILE *mgr_out = fopen_dir(gs_params->outputdir,"mgr_micro.txt","w");
  FILE *block_out = fopen_dir(gs_params->outputdir,"blocks.txt","w");
  print_block_header(block_out,anctable);

  /* Allocate maximum amount of space needed for analyzing strips.
   * The space can be reused for all blocks.
   */

  strip_mem = (int *) ckalloc(max_nanc, sizeof(int));

  gs_grimm_alloc(0,          /* unichromosomal */
		 max_nanc,   /* max # "genes" */
		 nspecies,   /* # species */

		 /* various GRIMM structures */
		 (distmem_t **) NULL, (graphstats_t ***) NULL,
		 &genome_list);

  /* store all compressed blocks in a table
   * and later determine which compressed blocks are the same perms
   */
  nblocks = count_blocks(anctree, blocktable, blockanc);
  compressed_blocks = (int **) ckalloc(nblocks, sizeof(int *));

  /* get list of all blocks, sorted by order */
  list_blocks(anctree, blocktable, blockanc,
	      &nblocks, &blockorders, &blockinds);

  /* examine all anchors by component */
  bo0 = blockorders[0];
  for (i=0; i<nblocks; i++) {
    m = bo0[i];
    b_name = blockinds[m*BL_tot + BL_off_l];
    c_v = blockinds[m*BL_tot + BL_off_v];
    k = blockinds[m*BL_tot + BL_off_k];
    block = anclist_b[c_v];                      /* block coords */
    nanc_in_block = block[FIELD_nanc];           /* # anchors in block */

    /* analyze block */
    gs_micro_mgr_1block(mgr_out,
			gs_params,
			unsigned_flag,
			b_name,
			block,
			blockanc, k, nanc_in_block, node2blockpos,
			anctable, anctree,
			blocktable,
			strip_mem, genome_list,
			compressed_blocks);


    /* write out info on block to block file */
    fprintf(block_out,"%d", b_name);
    print_block(gs_params, block_out, blocktable, block);
  }
  fclose(block_out);
  fclose(mgr_out);

  /* compare compressed blocks */
  compare_blocks_report(gs_params, nspecies, nblocks, compressed_blocks);


  /* free space */
  free((void *) blockinds);
  array2d_destroy(nspecies, (void **) blockorders);
  array2d_destroy(nblocks, (void **) compressed_blocks);
  gs_grimm_destroy(0,          /* unichromosomal */
		   max_nanc,   /* max # "genes" */
		   nspecies,   /* # species */

		   /* various GRIMM structures */
		   (distmem_t *) NULL, (graphstats_t **) NULL,
		   genome_list);
  free((void *) strip_mem);
}

/*
 *
 * ...
 * compressed_blocks[b_name-1] =
 *    { length, genome1order, genome2order, genome3order,... }
 *    as a 1-dimentsonal array
 */
void gs_micro_mgr_1block(
			 FILE *output,
			 GSparams *gs_params,
			 int unsigned_flag,
			 int b_name,
			 int *block,   /* coords of this block */
			 int **blockanc,
			 int base,
			 int nanc_in_block,
			 int *node2blockpos,
			 ANCTABLE *anctable, ANCTREE *anctree,
			 ANCTABLE *blocktable,
			 int *strip_mem,
			 struct genome_struct *genome_list,
			 int **compressed_blocks
			 )
{
  int *nodelist = anctree->nodelist;
  int nspecies = anctable->nspecies;
  int i,j;
  int m,s,g;
  int nanc1 = nanc_in_block - 1;
  int bpnum;

  int compressed_block_len;
  int *compressed_block;

  int **anclist_n =
    gs_params->perm_metric_use
    ? anctable->anclist2
    : anctable->anclist;

  char **snames = anctable->snames;

  /* 0: original
   * 1: modifications Qian wanted for matlab
   */
  int format_style = gs_params->format_style;


  /* compute strips
   * Let the anchors be numbered i=0,1,2,...,n-1 according to order in
   * species 0. 
   * 1. set strip_mem[i]=0 for all i, 0<=i<n.
   * 2. go through all markers in all species j>0;
   *    if (i,i+1) is a breakpoint in any species j>0 relative to species 0,
   *    set strip_mem[i]=1
   * 3. if strip_mem[i] != 0, set strip_mem[i] = # {k<=i: strip_mem[k] != 0}
   *    if strip_mem[i] == 0, set strip_mem[i] = -(next highest one of
   *                                               above form)
   * 4. In marker order for species j, delete anchor i if strip_mem[i]>0,
   *    else replace it by +/- strip_mem[i].
   */

  /* Form GRIMM encoding of the anchors in the block for all species */
  grimm_encode_block(block, blockanc, base, nanc_in_block, node2blockpos,
		     anctable, anctree, genome_list);


  /* determine strips */
  /* initialize: nothing is end of strip */
  memset((void *) strip_mem, 0, nanc_in_block * sizeof(int));
  strip_mem[nanc1] = 1; /* final anchor is end of its strip */

  for (j=1; j<nspecies; j++) {
    /* find strip boundaries in species j vs. species 0 */
    for (i=0; i<nanc_in_block; i++) {
      m = genome_list[j].genes[i];
      if (m > 0) {
	/* m is breakpoint unless next marker is m+1 */
	if (i == nanc1 || genome_list[j].genes[i+1] != m+1) {
	  strip_mem[m-1] = 1;
	}
      } else {
	/* (-m) is breakpoint unless previous marker is m-1 */
	if (i==0 || genome_list[j].genes[i-1] != m-1) {
	  strip_mem[-m-1] = 1;
	}
      }
    }
  }

  /* At each b.p., put cumulative count of b.p. up to that one */
  bpnum = 1;
  for (i=0; i<nanc_in_block; i++) {
    if (strip_mem[i]) {
      strip_mem[i] = bpnum++;
    } else {
      strip_mem[i] = -bpnum;
    }
  }


  /* write out info on block */
  if (!format_style) {
    fprintf(output,"# begin_block %d", b_name);
  } else {
    fprintf(output,"-%d", b_name);
  }
  print_block(gs_params, output, blocktable, block);


  if (!format_style) {
    fprintf(output,"# begin_anchors\n");
  }
  /* Write out coords of compressed anchors */
  for (i=0; i<nanc_in_block; i++) {
    /* microsegment # of which this anchor is a part */
    fprintf(output,"%d", abs(strip_mem[i]));

    /* print remainder of its coordinates */
    print_block(gs_params, output, anctable,
		anclist_n[nodelist[blockanc[0][base+i]]]);
  }
  if (!format_style) {
    fprintf(output,"# end_anchors\n");
  }

  /* store all compressed strips in a table */
  compressed_block_len = nspecies * (bpnum-1);
  compressed_block = compressed_blocks[b_name-1]
    = (int *) ckalloc(compressed_block_len+1, sizeof(int));
  *compressed_block++ = compressed_block_len;

  if (!format_style) {
    /* Write out MGR order */
    fprintf(output,"# begin_mgr\n");

    for (j=0; j<nspecies; j++) {
      /* fprintf(output,">genome%d\n", j+1); */
      fprintf(output,">%s\n", snames[j]);
      for (i=0; i<nanc_in_block; i++) {
	m = genome_list[j].genes[i];
	s = strip_mem[abs(m)-1];
	if (s>0) {
	  g = (m>0) ? s : -s;
	  fprintf(output,"%d ", g);
	  *compressed_block++ = g;
	}
      }
      fprintf(output,"\n");
    }
    fprintf(output,"# end_mgr\n");
    fprintf(output,"# end_block\n");
  }
}



/*****************************************************************************
 * Report on which compressed blocks have same MGR inputs
 *****************************************************************************/

static int **cbr_blocks;

void compare_blocks_report(GSparams *gs_params,
			   int nspecies,
			   int nblocks,
			   int **compressed_blocks)
{
  FILE *brep_out = fopen_dir(gs_params->outputdir,"mgr_micro_equiv.txt","w");

  int *blocklist = ckalloc(nblocks, sizeof(int));
  int i, i_next, j;
  int blen;

  fprintf(brep_out,
	  "# Groups of blocks with identical compressed anchor permutations\n");

  /* 0-based list of blocks */
  vec_012(blocklist, nblocks);  /* set to {0,1,2,...} */

  /* sort it in an order s.t.
   *   blocks with fewer markers come first;
   *   among blocks with same # of markers, all equal perms are adjacent
   */

  cbr_blocks = compressed_blocks;
  qsort((void *) blocklist, nblocks, sizeof(int), cbr_cmp_stable);

  /* write out equivalent blocks */
  /* i runs over 1st block in the group of equivalent blocks */
  for (i=0; i<nblocks; i=i_next) {

    /* i_next points to 1st block after the group of equivalent blocks */
    for (i_next = i+1; i_next < nblocks; i_next++) {
      if (cbr_cmp(&blocklist[i], &blocklist[i_next]))
	break;
    }

    blen = compressed_blocks[blocklist[i]][0] / nspecies;
    fprintf(brep_out,
	    "group with %d block%s of length %d:",
	    i_next-i, (i_next - i == 1) ? "" : "s", blen);

    
    for (j=i; j<i_next; j++) {
      fprintf(brep_out, " %d", blocklist[j]+1);
    }
    fprintf(brep_out, "\n");
  }

  fclose(brep_out);
}

/*
 * compare two compressed blocks;
 * for equal blocks, return 0 (so not stable sort)
 */
int cbr_cmp(const void *a0, const void *b0)
{
  int *a = cbr_blocks[*(int *) a0];
  int *b = cbr_blocks[*(int *) b0];

  int len_a = *a++;
  int len_b = *b++;

  int d = len_a - len_b;

  /* if different lengths, shorter one is smaller */
  if (d != 0)
    return d;

  /* same length, compare marker by marker */
  while (len_a-- > 0) {
    d = *a++ - *b++;
    if (d != 0)
      return d;
  }

  /* equal */
  return 0;
}

/*
 * compare two compressed blocks
 * for equal blocks, keep in same order (so stable sort)
 */
int cbr_cmp_stable(const void *a0, const void *b0)
{
  int d = cbr_cmp(a0,b0);
  if (d==0)
    d = *(int *)a0 - *(int *)b0;
  return d;
}



/*****************************************************************************
 * Report on macrorearrangements
 *****************************************************************************/

void gs_macro(
	      FILE *output,
	      GSparams *gs_params,
	      ANCTABLE *anctable, ANCTREE *anctree,
	      ANCTABLE *blocktable,
	      int **blockanc, int *node2blockpos
	      )
{
  int nspecies = anctable->nspecies;
  int fields_per_species = anctable->fields_per_species;
  int i,j;

  int max_chr = blocktable->max_chr;
  int nmarkers;

  distmem_t *distmem;
  graphstats_t **statsmatrix;
  struct genome_struct *genome_list;

  int nblocks;
  int **blockorders;
  int *blockinds;

  FILE *mgr_out = fopen_dir(gs_params->outputdir,"mgr_macro.txt","w");

  HTABLE **chrnames = anctable->chrnames;

  char **snames = anctable->snames;

  /* get list of all blocks, sorted by order */
  list_blocks(anctree, blocktable, blockanc,
	      &nblocks, &blockorders, &blockinds);

  /* TODO: unichromosomal, linear and circular */

  /* allocate space needed by GRIMM */
  gs_grimm_alloc(max_chr,
		 nblocks,
		 nspecies,

		 /* various GRIMM structures */
		 &distmem, &statsmatrix, &genome_list);

  nmarkers = nblocks + 2*max_chr;  /* # markers, incl. blocks + caps */

  /* encode block table into GRIMM format */
  for (j=0; j<nspecies; j++) {
    grimm_encode_mc(
		    nblocks,
		    max_chr,
		    j*fields_per_species + FIELD_base,
		    &genome_list[j],
		    blockorders[j],
		    blockinds,
		    blocktable);
  }


  /* compute distance matrix */
  for (i=0; i<nspecies; i++) {
    for (j=0; j<nspecies; j++) {
      mcdist_noncircular(&genome_list[i],
			 &genome_list[j],
			 nmarkers,
			 max_chr,
			 distmem,
			 &statsmatrix[i][j]);
      /* KLUDGE: # reuses in field "n" */
      statsmatrix[i][j].n =
	2*statsmatrix[i][j].d - (statsmatrix[i][j].bp_int +
				 statsmatrix[i][j].bp_ext);
    }
  }

  /* print distance matrices */
  fprintf(output,"Macrorearrangement analysis\n");
  fprintf(output,"# species:              %d\n", nspecies);
  fprintf(output,"# blocks:               %d\n", nblocks);
  fprintf(output,"# chromosomes/species:  ");
  for (j=0; j<nspecies; j++) {
    fprintf(output,"%d ", chrnames[j]->len);
  }
  fprintf(output,"\n");

  print_smatrix_field("Distance", d);
  print_smatrix_field("Number of Black Edges", bl);
  print_smatrix_field("Number of Cycles and Paths", cp);
  print_smatrix_field("Number of Gamma-Gamma Paths", pgg);
  print_smatrix_field("Number of Semi-knots", s);
  print_smatrix_field("Parameter r", r);
  print_smatrix_field("Parameter fr", fr);
  print_smatrix_field("Parameter gr", gr);

#if 0
  /* can't do this w/o forming capping/concat, which we didn't do */
  print_smatrix_field("Number of bad bonds", badbonds);
#endif

  print_smatrix_field("Number of internal breakpoints", bp_int);
  print_smatrix_field("Number of external breakpoints", bp_ext);

  print_smatrix_field("Number of breakpoint reuses", n); /* KLUDGE */


  /* encode block table into MGR file format */
  /* TODO: unichromosomal linear/circular; unsigned */
  fprintf(mgr_out,"# Options:\n");  /* TODO */
  fprintf(mgr_out,"# Macrorearrangements\n");
  for (j=0; j<nspecies; j++) {
    /* fprintf(mgr_out,">genome%d\n", j+1); */
    fprintf(mgr_out,">%s\n", snames[j]);
    mgr_macro_mc(
		    mgr_out,
		    nblocks,
		    j*fields_per_species + FIELD_base,
		    &genome_list[j],
		    chrnames[j],
		    blockorders[j],
		    blockinds,
		    blocktable);
  }

  free((void *) blockinds);
  array2d_destroy(nspecies, (void **) blockorders);
  fclose(mgr_out);
}

void grimm_encode_mc(
		      int nblocks,
		      int max_chr,
		      int base_field,
		      struct genome_struct *g,
		      int *blockorder_species,
		      int *blockinds,
		      ANCTABLE *blocktable)
{
  int cur_chr = -1;
  int cur_cap = nblocks+1;
  int max_cap = nblocks + 2*max_chr;
  int *gene_list = g->genes;
  int **anclist_b = blocktable->anclist;
  int *block;
  int chr, sign;
  int m;
  int b_name;
  int i;

  for (i=0; i<nblocks; i++) {
    m = blockorder_species[i];

    block = anclist_b[blockinds[m*BL_tot + BL_off_v]];
    b_name = blockinds[m*BL_tot + BL_off_l];
    chr = block[base_field + FIELDR_chr];
    sign = block[base_field + FIELDR_sign];

    if (chr != cur_chr) {
      /* put in caps */
      if (cur_chr >= 0) {
	/* close previous chromosome */
	*gene_list++ = cur_cap++;
      }
      /* open new chromosome */
      *gene_list++ = cur_cap++;

      cur_chr = chr;
    }

    /* store marker */
    *gene_list++ = sign * b_name;
  }

  /* close current chromosome,
   * insert null chromosomes out to required number of chromosomes
   */
  while (cur_cap <= max_cap)
    *gene_list++ = cur_cap++;
}


void mgr_macro_mc(
		  FILE *output,
		      int nblocks,
		      int base_field,
		      struct genome_struct *g,
		      HTABLE *chrnames_table,
		      int *blockorder_species,
		      int *blockinds,
		      ANCTABLE *blocktable)
{
  int cur_chr = -1;
  int **anclist_b = blocktable->anclist;
  int *block;
  int chr, sign;
  int m;
  int b_name;
  int i;

  for (i=0; i<nblocks; i++) {
    m = blockorder_species[i];

    block = anclist_b[blockinds[m*BL_tot + BL_off_v]];
    b_name = blockinds[m*BL_tot + BL_off_l];
    chr = block[base_field + FIELDR_chr];
    sign = block[base_field + FIELDR_sign];

    if (chr != cur_chr) {
      if (cur_chr != -1) {
	fprintf(output,"$\n");
      }
      cur_chr = chr;
      fprintf(output,"# Chromosome %s\n", hash_num2str(chrnames_table,cur_chr));
    }

    /* print marker */
    fprintf(output,"%d ", sign*b_name);
  }

  /* close current chromosome */

  fprintf(output,"$\n");
}

/*****************************************************************************
 * Overlap detection
 *****************************************************************************/

int gs_overlap(FILE *output,
	       GSparams *gs_params,
	       ANCTABLE *anctable, ANCTREE *anctree,
	       ANCTABLE *blocktable,
	       int **blockanc, int *node2blockpos
	       )
{
  int nspecies = anctable->nspecies;
  int fields_per_species = anctable->fields_per_species;

  int nblocks;
  int **blockorders;
  int *blockinds;
  int *b_0 = blockanc[0];

  int *block_r, chr_r, start_r, end_r, order_r, bnum_r;
  int *block_s, start_s, end_s, order_s, bnum_s;
  int b_name_r, b_name_s;
  int j, r, s;
  int base_field;

  int fix=0;  /* set to 1 if "fix" any blocks */
  int fixed_block;

  /* get list of all blocks, sorted by order */
  list_blocks(anctree, blocktable, blockanc,
	      &nblocks, &blockorders, &blockinds);

  for (j=0; j<nspecies; j++) {
    /* detect overlaps in species j */
    base_field = j*fields_per_species + FIELD_base;

    for (r=0; r<nblocks; r++) {
      order_r = blockorders[j][r];
      bnum_r = blockinds[order_r*BL_tot + BL_off_v];
      block_r = blocktable->anclist[bnum_r];
      b_name_r = blockinds[order_r*BL_tot + BL_off_l];
      chr_r = block_r[base_field + FIELDR_chr];
      start_r = block_r[base_field + FIELDR_start];
      end_r = start_r + block_r[base_field + FIELDR_len];


      /* if any overlaps, they will be in all of blocks r+1,r+2,...,r+m
       * for some m
       */
      for (s=r+1; s<nblocks; s++) {
	order_s = blockorders[j][s];
	bnum_s = blockinds[order_s*BL_tot + BL_off_v];
	block_s = blocktable->anclist[bnum_s];

	if (chr_r != block_s[base_field + FIELDR_chr])
	  break; /* different chromosomes */

	start_s = block_s[base_field + FIELDR_start];
	if (start_s >= end_r)
	  break;

	b_name_s = blockinds[order_s*BL_tot + BL_off_l];

	/* overlap (acceptable) or containment (bad)? */
	end_s = start_s + block_s[base_field + FIELDR_len];
	if (end_s >= end_r) {
	  fprintf(output,
		  "Overlap: Block %d overlaps %d in species %d\n",
		  b_name_r, b_name_s, j+1);
	} else {
	  if (start_s == start_r && end_s == end_r) {
	    fprintf(output,
		    "Containment: Block %d endpoints equal block %d in species %d, cannot fix\n",
		 b_name_r, b_name_s, j+1);
	  } else {
	    fprintf(output,
		    "Containment: Block %d contains block %d in species %d,",
		    b_name_r, b_name_s, j+1);
	    fixed_block = 
	      split_block(output,
			  anctable, anctree, blocktable,
			  b_0,
			  j,
			  blockinds[order_r*BL_tot + BL_off_v],
			  blockinds[order_s*BL_tot + BL_off_v],
			  blockinds[order_r*BL_tot + BL_off_k]
			  );
	    if (fixed_block) {
	      fprintf(output," split block\n");
	      fix = 1;
	    } else {
	      fprintf(output," unable to split block\n");
	    }
	  }
	}

#if 0
	print_block(gs_params, output, blocktable, block_r);
	print_block(gs_params, output, blocktable, block_s);
	fprintf(output,"\n");
#endif

	print_block2m(gs_params, output, blocktable, bnum_r);
	print_block2m(gs_params, output, blocktable, bnum_s);
	fprintf(output,"\n");
      }
    }
  }

  free((void *) blockinds);
  array2d_destroy(nspecies, (void **) blockorders);

  return fix;
}

/*****************************************************************************/

void print_species_coords(GSparams *gs_params, FILE *output,
			  ANCTABLE *anctable, int *coords_spec,
			  int species_no)
{
  HTABLE **chrnames = anctable->chrnames;
  char **signstrs = gs_params->signstrs;
  int chr_num = coords_spec[FIELDR_chr];
  int start = coords_spec[FIELDR_start];
  int len = coords_spec[FIELDR_len];
  int sign = coords_spec[FIELDR_sign];

  fprintf(output,
	  "%s %d %d %s",
	  hash_num2str(chrnames[species_no],chr_num),start,len,
	  signstrs[sign+1]);
}

void print_block_nocr(GSparams *gs_params, FILE *output,
		      ANCTABLE *anctable, int *anc)
{
  int j, base_field, chr_num, start, len, sign;
  int nspecies = anctable->nspecies;
  int fields_per_species = anctable->fields_per_species;
  HTABLE **chrnames = anctable->chrnames;
  char **signstrs = gs_params->signstrs;

  if (anc == (int *) NULL)
    return;

  for (j=0; j<nspecies; j++) {
    base_field = j*fields_per_species + FIELD_base;
    chr_num = anc[base_field + FIELDR_chr];
    start = anc[base_field + FIELDR_start];
    len = anc[base_field + FIELDR_len];
    sign = anc[base_field + FIELDR_sign];

    fprintf(output,
	    " %s %d %d %s",
	    hash_num2str(chrnames[j],chr_num),start,len,
	    signstrs[sign+1]);
  }
}

void print_block(GSparams *gs_params, FILE *output,
		 ANCTABLE *anctable, int *anc)
{
  print_block_nocr(gs_params, output, anctable, anc);
  fprintf(output,"\n");
}

/* print block i using both metrics */
void print_block2m(GSparams *gs_params, FILE *output,
		   ANCTABLE *anctable, int i)
{
  int **anclist = anctable->anclist;
  int **anclist2 = anctable->anclist2;

  print_block(gs_params, output, anctable, anclist[i]);
  if (anclist2 != (int **) NULL
      && anclist2[i] != (int *) NULL) {
    fprintf(output,"[");
    print_block_nocr(gs_params, output, anctable, anclist2[i]);
    fprintf(output,"]\n");
  }
}

/* print comments with species names */
void print_block_header(FILE *output, ANCTABLE *anctable)
{
  int nspecies = anctable->nspecies;
  char **snames = anctable->snames;
  int i;

  for (i=0; i<nspecies; i++) {
    fprintf(output, "# genome%d: %s\n", i+1, snames[i]);
  }

  fprintf(output,"# block_id");
  for (i=1; i<=nspecies; i++) {
    fprintf(output,
	    " genome%d_chr genome%d_start genome%d_len genome%d_sign",
	    i, i, i, i);
  }
  fprintf(output, "\n");
}

/*****************************************************************************/

/* Block A contains block B in species j.
 * Split block A into 2-3 parts:
 * part before B; intersecting B (may be null); after B
 *
 * TODO: if anchors of A overlap the boundaries of B, may need to
 * handle it differently.
 *
 * return:
 *   1 = successfully split block
 *   0 = unable to split block
 */
int split_block(FILE *report_file,
		ANCTABLE *anctable, ANCTREE *anctree,
		ANCTABLE *blocktable,
		int *b_0,
		int j, /* species # */
		int c_A, int c_B,
		int A_startanc
		)
{
  int *nodelist = anctree->nodelist;
  int fields_per_species = anctable->fields_per_species;
  int **anclist = anctable->anclist;

  /* for anchors */
  int base_field = j*fields_per_species + FIELD_base;
  int col_start = base_field + FIELDR_start;

  /* for blocks */
  int col_start_b = col_start;

  int *comp = anctree->comp;
  int *block_A = blocktable->anclist[c_A];
  int *block_B = blocktable->anclist[c_B];

  /* coords of block B */
  int start_B = block_B[col_start_b];
  int end_B = start_B + block_B[base_field + FIELDR_len];
  

  int root1 = MAXINT;  /* part before B */
  int root2 = MAXINT;  /* part intersecting B */
  int root3 = MAXINT;  /* part after B */

  int nroot1=0, nroot2=0, nroot3=0;  /* # verts with each root */

  /* v loops over nodes in block A
   * c_A is root of block A
   */
  int v;
  int *anc_v;
  int start_v;

  int nanc_in_block = block_A[FIELD_nanc];
  int A_endanc = A_startanc + nanc_in_block;
  int k;

  /* loop over anchors in block
   * to determine least vertices in each of the three sections
   */
  for (k=A_startanc; k<A_endanc; k++) {
    v = b_0[k];
    /* if this block already underwent splitting during this pass,
     * abort splitting
     */
    if (comp[v] != c_A) {
      fprintf(report_file, 
	      " (v=%d, c[v]=%d, c_A=%d, c[c_A]=%d, c_B=%d, c[c_B]=%d; block already split?)",
	      v, comp[v], c_A, comp[c_A], c_B, comp[c_B]);
      return 0;
    }

    anc_v = anclist[nodelist[v]];
    start_v = anc_v[col_start];
    if (start_v < start_B) {
      /* part before block B */
      root1 = min2(root1,v);
      nroot1++;
    } else if (start_v >= end_B) {
      /* part after block B */
      root3 = min2(root3,v);
      nroot3++;
    } else {
      /* part intersecting block B */
      root2 = min2(root2,v);
      nroot2++;
    }
  }

  /* DEBUGGING */
  fprintf(report_file," anchor counts (%d,%d,%d),", nroot1, nroot2, nroot3);

  if (nroot1 == nanc_in_block
      || nroot2 == nanc_in_block
      || nroot3 == nanc_in_block) {
    /* unable to split block */
    return 0;
  }

  /* now reassign the clusters */
  for (k=A_startanc; k<A_endanc; k++) {
    v = b_0[k];
    anc_v = anclist[nodelist[v]];
    start_v = anc_v[col_start];
    if (start_v < start_B) {
      /* part before block B */
      comp[v] = root1;
    } else if (start_v >= end_B) {
      /* part after block B */
      comp[v] = root3;
    } else {
      /* part intersecting block B */
      comp[v] = root2;
    }
  }

  /* successfully split block */
  return 1;
}
