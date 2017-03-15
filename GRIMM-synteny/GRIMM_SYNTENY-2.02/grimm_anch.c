/* grimm_anch.c
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

/* Last modified on Wed Sep 6, 2006, by Glenn Tesler
 */

/* GRIMM-Anchors */
/*
 * The elements of anctable are not valid anchors,
 * must filter and combine them to get non-conflicting coordinates
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ckalloc.h"
#include "time_stamp.h"
#include "hash.h"
#include "anctable.h"
#include "gsy.h"
#include "grimm_anch.h"
#include "gsfile.h"
#include "gscomp.h"


/* Initial array size to allocate for repeat families */
#define REPEAT_LEN 500
/* If repeat families have more intervals than initial allocation size,
 * increase allocation size by this increment as many times as needed */
#define REPEAT_INC 500

/* In case of repeats with inconsistent signs,
 * this is max # trials to compute signs.
 * If signs are consistent, they will be found on first try.
 */
#define REPEAT_SIGN_TRIES 5

/* should be FIELD_base + 0*fields_per_species + FIELDR_order */
#define FIELD_mate (FIELD_base + FIELDR_order)




void ga_tree(FILE *report_file,
	     ANCTABLE *anctable, GSparams *gs_params,
	     int **orders,
	     GANCTREE **ganctree_ret
	     );

GANCTREE *ga_tree_alloc(int nspecies, int nanc);
void ga_tree_destroy(GANCTREE *ganctree);
int ga_is_edge(ANCTABLE *anctable, GANCTREE *ganctree, GSparams *gs_params,
	       int j, int v, int w);
void ga_comp_join(GANCTREE *ganctree, int j, int v, int w);
int ga_comp_get(GANCTREE *ganctree, int j, int v);
void ga_comp_optimize(GANCTREE *ganctree);


void ga_sort(FILE *report_file, ANCTABLE *anctable,
	     int ***orders_ret);
int SGA_cmp(const void *a, const void *b);

void ga_out_anchors(FILE *report_file,
		    ANCTABLE *anctable,
		    GSparams *gs_params,
		    GANCTREE *ganctree);

/* parameters that will stay the same through numerous calls to
 * ga_out_repeat and its further subroutines
 */
typedef struct {
  FILE *report_file;
  FILE *repeat_out;
  ANCTABLE *anctable;
  GSparams *gs_params;
  GANCTREE *ganctree;

  /* list of nodes in the component */
  int comp_nodes_len;    /* # elts in list */
  int *comp_nodes;       /* the list */

  /* space to construct coordinates of intervals/anchors */
  int *scratch_anc;

  int repeat_id;

  /* repeat interval info table
   * Let S(m) be all nodes with given component in particular species
   *
   *    intervals[m] = {
   *        component,
   *        species,
   *        chr,
   *        start,       // minimum of start over S(m)
   *        end,         // maximum of end over S(m)
   *        sign,        // =0 if haven't guessed a sign yet
   *                     // =+1 or -1 if have guessed a sign
   *                     // When sign !=0,
   *        count_plus,  //   # nodes in S(m) where sign in this species
   *                     //   is the same as the guessed sign
   *        count_minus  //      ... opposite ...
   *    }
   *    for m = 0,1,2,...,(len0+len1+... - 1)
   *    where there are len_j intervals in species j
   */
  
  int intervals_alloc;   /* # rows allocated in intervals array */
  int intervals_len;     /* # rows used in intervals array */
  int **intervals;
} REPEAT_PARAMS_t;
#define RP_component 0
#define RP_species 1
#define RP_chr 2
#define RP_start 3
#define RP_end 4
#define RP_sign 5
#define RP_ct_p 6
#define RP_ct_m 7
#define RP_size 8

void ga_out_repeat(REPEAT_PARAMS_t *repeat_params,
		   int merged
		   );

void ga_repeat_signs(REPEAT_PARAMS_t *repeat_params);

int ga_repeat_signs_interval(REPEAT_PARAMS_t *repeat_params,
			     int seed_m);

int interval_sign(REPEAT_PARAMS_t *repeat_params,
		  int *anc_i,
		  int w,
		  int flag_ambiguous);


int interval_cmp(const void *a, const void *b);

int *get_interval(REPEAT_PARAMS_t *repeat_params,
		  int component,
		  int species);


int same_species_mate(REPEAT_PARAMS_t *repeat_params,
		      int v, int s);
int same_species_mate_comp(REPEAT_PARAMS_t *repeat_params,
			   int v, int s);

/*****************************************************************************/





void grimm_anch(FILE *report_file, ANCTABLE *anctable, GSparams *gs_params)
{
  GANCTREE *ganctree;
  int **orders;


  /* sort nodelist for each species */
  time_stamp(report_file,"Sort for each species");
  ga_sort(report_file, anctable,
	  &orders);


  /* form components */
  time_stamp(report_file,"Determine conflicting alignments");
  ga_tree(report_file,anctable,gs_params,orders,
	  &ganctree);

  time_stamp(report_file, "Output anchors and repeats");
  ga_out_anchors(report_file, anctable, gs_params, ganctree);

  time_stamp(report_file,"Free space");
  ga_tree_destroy(ganctree);
  array2d_destroy(anctable->nspecies, (void **) orders);
  time_stamp(report_file,"Done");
}

void ga_tree(FILE *report_file,
	     ANCTABLE *anctable, GSparams *gs_params,
	     int **orders,
	     GANCTREE **ganctree_ret
	     )
{
  int nanc = anctable->nanc;
  int nspecies = anctable->nspecies;
  int **anclist = anctable->anclist;
  int fields_per_species = anctable->fields_per_species;
  
  int same_species = gs_params->same_species;

  int v,w;
  int k_v,k_w;
  int j;
  int *o_j;
  int num_edges;
  int *nodelist = anctable->orders[0];

  int *anc_v, *anc_w;

  /* coords of interval being built up */
  int chr1,start1,len1,end1;

  /* coords of next anchor */
  int chr2,start2,len2,end2;

  int base_field;

  /* form a forest on all the anchors, in order by first genome */
  GANCTREE *ganctree = ga_tree_alloc(nspecies, nanc);
  ganctree->nodelist = nodelist;

  *ganctree_ret = ganctree;

  /* if species-on-self, join each node with its mates */
  if (same_species) {
    for (k_v = 0; k_v < nanc; k_v += nspecies) {
      /* in original input order, mates are consecutive */
      v = anclist[k_v][FIELD_mate];      /* note: nodelist[k_v] = v   */
      for (k_w = k_v+1; k_w < k_v + nspecies; k_w++) {
	w = anclist[k_w][FIELD_mate];    /*       nodelist[k_w] = w */

	ga_comp_join(ganctree,-1,v,w);  /* join components globally */
      }
    }
  }


  /* loop over potential edges */
  for (j=0; j<nspecies; j++) {
    base_field = j*fields_per_species + FIELD_base;

    o_j = orders[j];
    num_edges = 0;

    for (k_v=0; k_v<nanc; ) {
      v = o_j[k_v];
      anc_v = anclist[nodelist[v]];
      chr1 = anc_v[base_field + FIELDR_chr];
      start1 = anc_v[base_field + FIELDR_start];
      len1 = anc_v[base_field + FIELDR_len];
      end1 = start1 + len1;

      for (k_w = k_v+1; k_w<nanc; k_w++) {
	w = o_j[k_w];
	anc_w = anclist[nodelist[w]];

	chr2 = anc_w[base_field + FIELDR_chr];
	if (chr1 != chr2) break;

	start2 = anc_w[base_field + FIELDR_start];
	if (start2 >= end1) break;

	len2 = anc_w[base_field + FIELDR_len];
	end2 = start2 + len2;
	end1 = max2(end1,end2);

	/* it's an edge!
	 * join the anchors in species j (adjacent anchors)
	 * join the anchors in species "-1" (global repeat family)
	 */
	ga_comp_join(ganctree,j,v,w);
	ga_comp_join(ganctree,-1,v,w);
	num_edges++;
      }
      k_v = k_w;

#if 0
      if (k_w - k_v > 100) {
	fprintf(report_file,"Warning: not pruning? j=%d v=%d k_v=%d k_w=%d\n", j, v, k_v, k_w);
      }
#endif
    }

    fprintf(report_file,"%d conflicts in species %d\n", num_edges, j+1);
  }

  ga_comp_optimize(ganctree);
}

void ga_tree_old(FILE *report_file,
	     ANCTABLE *anctable, GSparams *gs_params,
	     int **orders,
	     GANCTREE **ganctree_ret
	     )
{
  int nanc = anctable->nanc;
  int nspecies = anctable->nspecies;
  int edge_status;
  int v,w;
  int k_v,k_w;
  int j;
  int *o_j;
  int num_edges;

  /* form a forest on all the anchors, in order by first genome */
  GANCTREE *ganctree = ga_tree_alloc(nspecies, nanc);
  ganctree->nodelist = anctable->orders[0];

  *ganctree_ret = ganctree;

  /* loop over potential edges */
  for (j=0; j<nspecies; j++) {
    o_j = orders[j];
    num_edges = 0;

    for (k_v=0; k_v<nanc; k_v++) {
      v = o_j[k_v];

      for (k_w = k_v+1; k_w<nanc; k_w++) {
	w = o_j[k_w];

	edge_status = ga_is_edge(anctable,ganctree,gs_params,j,v,w);
	/* will there possibly be ANY more edges (v,w') with
	 * w' later than w in species j?
	 */

	if (edge_status == E_TOOHIGH)
	  break;  /* no */

	if (edge_status == E_NOTEDGE)
	  continue; /* maybe */

	/* it's an edge!
	 * join the anchors in species j (adjacent anchors)
	 * join the anchors in species "-1" (global repeat family)
	 */
	ga_comp_join(ganctree,j,v,w);
	ga_comp_join(ganctree,-1,v,w);
	num_edges++;
      }
#if 0
      if (k_w - k_v > 100) {
	fprintf(report_file,"Warning: not pruning? j=%d v=%d k_v=%d k_w=%d\n", j, v, k_v, k_w);
      }
#endif
    }

    fprintf(report_file,"%d conflicts in species %d\n", num_edges, j+1);
  }

  ga_comp_optimize(ganctree);
}

/*
 * Allocate space for nspecies+1 spanning forests with n nodes
 * node i is anchor # nodelist[i]
 * where nodelist = anctable->orders[0]
 *
 * Don't need whole forests, just the components, identified
 * by the lowest numbered anchor in the component
 *
 * ganctree->comp[species+1][v] = w
 *         = component number of node v in given species
 *           (v is the root if v == w)
 * species = -1,0,1,...,nspecies-1
 * v = 0,1,2,...,nanc-1
 *
 * In a given species>=0:
 *     A component is transitive closure of overlapping repeats.
 * In "species=-1":
 *     Merge components from all species to get repeat families.
 */

GANCTREE *ga_tree_alloc(int nspecies, int nanc)
{
  int j;
  GANCTREE *ganctree = (GANCTREE *) ckalloc(1, sizeof(GANCTREE));

  ganctree->n = nanc;
  ganctree->nspecies = nspecies;
  ganctree->nodelist = (int *) NULL;
  ganctree->comp = int_array2d_create(nspecies+1, nanc);

  for (j=0; j<=nspecies; j++) {
    /* set to list 0,1,2,3,...,nanc-1 */
    vec_012(ganctree->comp[j], nanc);
  }

  return ganctree;
}

void ga_tree_destroy(GANCTREE *ganctree)
{
  array2d_destroy(ganctree->nspecies+1, (void **) ganctree->comp);
  free((void *) ganctree);
}


/*
 * Determine if (v,w) is an edge in species j
 *
 * Return:
 *  E_ISEDGE: it's an edge
 *  E_TOOHIGH: it's not an edge, and the separation in species j is
 *             so high that there will not be any more edges (v,w')
 *             with anchor w' later than anchor w in species j
 *             (assuming anchor v preceeds anchor w in species j)
 *  E_NOTEDGE: it's not an edge
 *
 * TODO: circular chromosomes; unsigned anchors
 */
int ga_is_edge(ANCTABLE *anctable, GANCTREE *ganctree, GSparams *gs_params,
	       int j, int v, int w)
{
  int fields_per_species = anctable->fields_per_species;
  int *nodelist = ganctree->nodelist;
  int base_field = j*fields_per_species + FIELD_base;

  /* anchor coordinates */
  int chr_v, start_v, len_v, end_v; /*, sign_v;*/
  int chr_w, start_w, len_w, end_w; /*, sign_w;*/

  /* data for the anchors */
  int *anc_v = anctable->anclist[nodelist[v]];
  int *anc_w = anctable->anclist[nodelist[w]];

  /* chromosome, sign, start, end for each anchor */
  chr_v = anc_v[base_field + FIELDR_chr];
  start_v = anc_v[base_field + FIELDR_start];
  len_v = anc_v[base_field + FIELDR_len] - 1;
  end_v = start_v + len_v;
#if 0
  sign_v = anc_v[base_field + FIELDR_sign];
#endif

  chr_w = anc_w[base_field + FIELDR_chr];
  start_w = anc_w[base_field + FIELDR_start];
  len_w = anc_w[base_field + FIELDR_len] - 1;
  end_w = start_w + len_w;
#if 0
  sign_w = anc_w[base_field + FIELDR_sign];
#endif

  /* different chromosomes, not an edge */
  if (chr_v != chr_w) {
    return E_TOOHIGH;
  }

  /* one anchor ends before other starts, not an edge */
  if (end_v <= start_w
      || end_w <= start_v)
    return E_TOOHIGH;

  /* edges definitely overlap */
  return E_ISEDGE;

#if 0
  /* if option chosen "edge is defined by any overlap", then it's an edge */
  if (gs_params->ga_edge_olap)
    return E_ISEDGE;

  /* option chosen is "edge is defined by containment" */
  if ((start_w <= start_v && end_v <= end_w)
      || (start_v <= start_w && end_w <= end_v))
    return E_ISEDGE;
  else
    return E_NOTEDGE;
#endif
}


/* join the components with nodes v and w in species j */
void ga_comp_join(GANCTREE *ganctree, int j, int v, int w)
{
  int root_v = ga_comp_get(ganctree,j,v);
  int root_w = ga_comp_get(ganctree,j,w);
  int root = min2(root_v,root_w);
  ganctree->comp[j+1][root_v] = ganctree->comp[j+1][root_w] = root;
}

/*
 * comp_get: get the component # with node v for species j
 * updates intermediate nodes for efficiency
 */
int ga_comp_get(GANCTREE *ganctree, int j, int v)
{
  int v0 = v;
  int w;
  int **comp = ganctree->comp;
  int j1=j+1;

  /* follow v, comp[v], comp[comp[v]], ... till it stabilizes
   */
  while (v != (w=comp[j1][v])) {
    v = w;
  }

  /* then update all the pointers to point to the new component root
   */
  while (v0 != (w=comp[j1][v0])) {
    comp[j1][v0] = v;
    v0 = w;
  }

  return v;
}

/*
 * optimize all component pointers to point to component root
 */
void ga_comp_optimize(GANCTREE *ganctree)
{
  int n = ganctree->n;
  int nspecies = ganctree->nspecies;
  int **comp = ganctree->comp;
  int *comp_j;

  int i,j;
  for (j=0; j<=nspecies; j++) {
#if 0
    fprintf(stderr,"species %d:",j);
#endif
    comp_j = comp[j];
    for (i=0; i<n; i++) {
      comp_j[i] = comp_j[comp_j[i]];
#if 0
      fprintf(stderr," %d", comp_j[i]);
#endif
    }
#if 0
    fprintf(stderr,"\n");
#endif
  }
}



/*
 * SGA = sort for GRIMM-Anchors
 */

static int **SGA_anclist;
static int *SGA_nodelist;
static int SGA_base;

/*
 * Sort nodelist by order in each species
 */

void ga_sort(FILE *report_file, ANCTABLE *anctable,
	     int ***orders_ret)
{
  int nspecies = anctable->nspecies;
  int fields_per_species = anctable->fields_per_species;
  int n = anctable->nanc;
  int i,j;
  int *nodelist = anctable->orders[0];
  int *b_j;

  /*
   * TODO: anctable->orders[j] isn't suitable, so don't compute it
   */
  int **orders = int_array2d_create(nspecies, n);
  *orders_ret = orders;
  

  for (j=0; j<nspecies; j++) {
    b_j = orders[j];

    for (i=0; i<n; i++)
      b_j[i] = i;

    /* species 0: nodelist is in order already */
    if (j==0) continue;

    /* sort it in order for species j by
     *   (chromosome #, start, len, node #)
     */

    SGA_anclist = anctable->anclist;
    SGA_nodelist = nodelist;
    SGA_base = j*fields_per_species + FIELD_base;
    qsort((void *) b_j, n, sizeof(int), SGA_cmp);
  }
}

/* TODO:
 * this is same as SB_cmp in gscomp.c,
 * separate out sorting into a separate file
 */

int SGA_cmp(const void *a, const void *b)
{
  int *anc1 =
    SGA_anclist[SGA_nodelist[*(int *)a]];
  int *anc2 =
    SGA_anclist[SGA_nodelist[*(int *)b]];

  /* compare chromosomes */
  int d = anc1[SGA_base + FIELDR_chr] - anc2[SGA_base + FIELDR_chr];

  if (d == 0) {
    /* break tie with start coordinates */
    d = anc1[SGA_base + FIELDR_start] - anc2[SGA_base + FIELDR_start];

    if (d == 0) {
      /* break tie with length */
      d = anc1[SGA_base + FIELDR_len] - anc2[SGA_base + FIELDR_len];

      if (d == 0) {
	/* break tie with node number */
	d = *(int *)a - *(int *)b;
      }
    }
  }

  return d;
}

/* sort in order by (global component, node #) */
static int *COMPLIST_comp;
int COMPLIST_cmp(const void *a, const void *b)
{
  int node1 = COMPLIST_comp[*((int *)a)];
  int node2 = COMPLIST_comp[*((int *)b)];

  /* compare chromosomes */
  int d = node1 - node2;

  if (d == 0) {
    /* break tie with node number */
    d = (*((int *)a)) - (*((int *)b));
  }

  return d;
}


void ga_out_anchors(FILE *report_file,
		    ANCTABLE *anctable,
		    GSparams *gs_params,
		    GANCTREE *ganctree)
{
  /* strategy for anchors only, no analysis of repeats */

  int nanc = anctable->nanc;
  int nspecies = anctable->nspecies;
  int **anclist = anctable->anclist;
  int fields_per_anchor = anctable->fields_per_anchor;
  int bytes_per_anchor = fields_per_anchor * sizeof(int);
  int fields_per_species = anctable->fields_per_species;
  int condense_strips = gs_params->condense_strips;
  int same_species = gs_params->same_species;

  int *nodelist = ganctree->nodelist;

  int *c_r;  /* merged components of all species */
  int v, len, i, j;
  int start0, len0, end0, start_i, len_i, end_i;
  int good;
  int *anc_i;
  int base_field;

  int q, q_next, w;



#if 0
  /* count # anchors per component in species "-1" */
  int *comp_cts = (int *) ckalloc(nanc, sizeof(int));
#endif
  /* sorted list of nodes by component */
  int *comp_list = (int *) ckalloc(nanc, sizeof(int));

  int **comp = ganctree->comp;
  int *comp_j;

  /* merge coords of different alignments into one anchor */
  int *merged_anch = (int *) ckalloc(fields_per_anchor, sizeof(int));



  int count_single = 0;
  int count_merged = 0;
  int count_comp = 0;
  int count_used = 0;

  FILE *anc_out = fopen_dir(gs_params->outputdir,"unique_coords.txt","w");
  FILE *repeat_out = fopen_dir(gs_params->outputdir,"repeat_coords.txt","w");

  REPEAT_PARAMS_t repeat_params_mem, *repeat_params=&repeat_params_mem;

  /* initialize repeat_params fields */
  repeat_params->report_file = report_file;
  repeat_params->repeat_out = repeat_out;
  repeat_params->anctable = anctable;
  repeat_params->gs_params = gs_params;
  repeat_params->ganctree = ganctree;

  repeat_params->intervals_alloc = REPEAT_LEN;
  repeat_params->intervals_len = 0;
  repeat_params->intervals =
    int_array2d_create(REPEAT_LEN,RP_size);

  /* repeat_params->comp_nodes = ??? */

  repeat_params->scratch_anc =
    (int *) ckalloc(fields_per_anchor, sizeof(int));

  repeat_params->repeat_id = 1;

  /* Unbuffer the output streams */
  setvbuf(anc_out, (char *) NULL, _IONBF, (size_t) 0);
  setvbuf(repeat_out, (char *) NULL, _IONBF, (size_t) 0);


  /* print header with species names */
  print_block_header(anc_out,anctable);


#if 0
  /* count # anchors per component in species "-1" */
  c_r = ganctree->comp[0];
  for (v=0; v<nanc; v++) {
    comp_cts[c_r[v]]++;
  }
#endif

  /* initialize list of nodes: 0,1,2,... */
  vec_012(comp_list, nanc);
  /* sort it in order by (global component, node #) */
  c_r = COMPLIST_comp = ganctree->comp[0];
  qsort((void *) comp_list, nanc, sizeof(int), COMPLIST_cmp);

  /* Go through components in order by species "-1".
   * If each component is contiguous in all species,
   * and if signs are consistent,
   * then output it.
   */

  for (q_next=0, q=0; q<nanc; q=q_next) {
    count_comp++;

    /* global component number */
    v = c_r[comp_list[q]];

    /* find all entries in component list for this same global component */
    q_next = q+1;
    while (q_next < nanc
	   && c_r[comp_list[q_next]] == v)
      q_next++;
    len = q_next - q;

#if 0
    //  for (v=0; v<nanc; v++) {
    /* ^^^ */
    /* only do species "-1" component roots */
    if (c_r[v] != v) {
      continue;
    }
    count_comp++;
    len = comp_cts[v];
#endif

    /* build merged anchor, starting with 1st alignment */
    memcpy((void *) merged_anch,
	   (const void *) anclist[nodelist[v]],
	   bytes_per_anchor);

    /*
     * If all the nodes are all "same component"
     * and have consistent signs, it's mergeable.
     *
     * "same component":
     *   with condense strips option: component in species "-1"
     *   w/o: per-species component
     */

    /* good=0: not mergeable
     *      1: mergeable
     *      2: all alignments are identical (duplicate entries in input)
     */
    good = 2;

    /* can only be mergeable if all nodes in component are consecutive */
    if (comp_list[q+len-1] - comp_list[q] != len-1) {
      good = 0;
    } else for (i=0; i<len; i++) {
      w = comp_list[q+i];
#if 0
      anc_i = anclist[nodelist[v+i]];
#endif
      anc_i = anclist[nodelist[w]];

      for (j=0; j<nspecies; j++) {
	base_field = j*fields_per_species + FIELD_base;
	comp_j = condense_strips ? comp[0] : comp[j+1];

	/* bad component or wrong sign? */
#if 0
	if (comp_j[v+i] != v
#endif
	if (comp_j[w] != v
	    || (anc_i[base_field + FIELDR_sign]
		!= merged_anch[base_field + FIELDR_sign])) {
	  good = 0; /* not mergeable */
	  break;
	}

	start0 = merged_anch[base_field + FIELDR_start];
	len0 = merged_anch[base_field + FIELDR_len];
	end0 = start0 + len0;

	start_i = anc_i[base_field + FIELDR_start];
	len_i = anc_i[base_field + FIELDR_len];
	end_i = start_i + len_i;

	if (start0 != start_i || len0 != len_i)
	  good=1;

	merged_anch[base_field + FIELDR_start] = start0 = min2(start0,start_i);
	end0 = max2(end0,end_i);
	merged_anch[base_field + FIELDR_len] = end0 - start0;
      }

      if (!good) break;
    }

    /* for repeats, and for merged anchors, output analysis to
     * the repeats file.
     */

    if (len>1 && good<2) {
      repeat_params->comp_nodes = &comp_list[q];
      repeat_params->comp_nodes_len = len;

      ga_out_repeat(repeat_params,

		    /* parameters for this specific repeat family */
		    good /* 1=merged into one "anchor", 0=not merged */
		    );
    }

    if (!good) continue;

    count_used += len;

    /* if all entries are duplicated, just use first one */
    if (good == 2) { good = len = 1; }

    if (len == 1)
      count_single++;
    else
      count_merged++;

    /* output anchor */
#if 0
    fprintf(anc_out,"0"); /* ID # not used */
#endif
    /* ID # from orig anchor */
    fprintf(anc_out,"%d",merged_anch[FIELD_id]);

    /* (merged) anchor coordinates */
    print_block_nocr(gs_params, anc_out, anctable, merged_anch);

    /* comment if merged */
    if (len>1) {
      fprintf(anc_out, " # merged; repeat %d", repeat_params->repeat_id - 1);
    }
    fprintf(anc_out,"\n");
  }

  fclose(repeat_out);
  fclose(anc_out);

  fprintf(report_file,"%d anchors output:\n", count_single + count_merged);
  fprintf(report_file,"   %d as-is anchors\n",count_single);
  fprintf(report_file,"   %d merged anchors\n",count_merged);
  fprintf(report_file,"   Used %d input alignments for those\n",
	  same_species ? count_used/nspecies : count_used);
  fprintf(report_file,"%d repeats output:\n", repeat_params->repeat_id - 1);
  fprintf(report_file,"   %d merged anchors\n",count_merged);
  fprintf(report_file,"   %d non-merged repeat families\n",
	  repeat_params->repeat_id - 1 - count_merged);
  fprintf(report_file,"%d components:\n", count_comp);
  fprintf(report_file,"   %d as-is anchors\n",count_single);
  fprintf(report_file,"   %d merged anchors\n",count_merged);
  fprintf(report_file,"   %d non-merged repeat families\n",
	  repeat_params->repeat_id - 1 - count_merged);

  free((void *) comp_list);
  free((void *) repeat_params->scratch_anc);
  array2d_destroy(repeat_params->intervals_alloc,
		  (void **) repeat_params->intervals);

  free((void *) merged_anch);
}

void ga_out_repeat(REPEAT_PARAMS_t *repeat_params,
		   int merged
		   )
{
  FILE *report_file = repeat_params->report_file;
  FILE *repeat_out = repeat_params->repeat_out;
  ANCTABLE *anctable = repeat_params->anctable;
  GSparams *gs_params = repeat_params->gs_params;
  GANCTREE *ganctree = repeat_params->ganctree;
  int *scratch_anc = repeat_params->scratch_anc;

  int *comp_nodes = repeat_params->comp_nodes;
  int len = repeat_params->comp_nodes_len;

  int **intervals = repeat_params->intervals;
  int intervals_alloc = repeat_params->intervals_alloc;
  int intervals_len = 0;
  int m;


  int nspecies = anctable->nspecies;
  int nanc = anctable->nanc;
  int **anclist = anctable->anclist;
  int fields_per_anchor = anctable->fields_per_anchor;
  int bytes_per_anchor = fields_per_anchor * sizeof(int);
  int fields_per_species = anctable->fields_per_species;

  int *nodelist = ganctree->nodelist;
  int **comp = ganctree->comp;

  int i,j;
  int *anc_i;
  int base_field;
  int chr_ij, start_ij, len_ij, end_ij, sign_ij, comp_ij;
  int w;
  int *cur_int;
  int *c_r = ganctree->comp[0];  /* merged components of all species */

  int anc_sign, ambig;
  int count;

  int same_species = gs_params->same_species;
  int nspecies_s = same_species ? 1 : nspecies;


  /****************************************************************
   * Pre-allocate space for intervals array if chance it's too
   * small.
   * Need to count # components for each species.
   ****************************************************************/

  if (repeat_params->intervals_alloc + REPEAT_INC <= len) {
    /* There's a chance we'll need more space than one increment
     * will create.
     */

    /* Determine exact number of intervals */
    count = 0;

    /* Count # distinct components for species j.
     * Do this by re-sorting the node list.
     * Then put back into normal order (j = -1).
     * (In same species mode, only do count for j=0, then put in order
     * for j=-1.)
     */
    for (j = nspecies_s - 1;
	 j>=-1;
	 j--) {
      COMPLIST_comp = ganctree->comp[j+1];
      qsort((void *) comp_nodes, len, sizeof(int), COMPLIST_cmp);

      /* Final pass is just to put nodes back in order */
      if (j==-1) break;

      count++;  /* component for i=0 node */

      for (i=1; i<len; i++) {
	if (COMPLIST_comp[comp_nodes[i-1]] != COMPLIST_comp[comp_nodes[i]]) {
	  /* bump if new component for node i in species j */
	  count++;
	}
      }
    }

    if (intervals_alloc < count) {
      intervals =
	int_array2d_resize(intervals,
			   intervals_alloc,
			   RP_size,
			   count);
      repeat_params->intervals = intervals;
      intervals_alloc = 
	repeat_params->intervals_alloc = count;
    } 
  }

  /* Create intervals array
   *  intervals[m] = { species, component,  // set up once and for all
   *                   start, end,          // initial values
   *                   sign, count_plus, count_minus // 0
   *                 }
   */

  /****************************************************************
   * Determine coordinates of intervals in each species
   ****************************************************************/

  repeat_params->intervals_len = 0;
  for (i=0; i<len; i++) {
    w = comp_nodes[i];
    anc_i = anclist[nodelist[w]];
    for (j=0; j<nspecies_s; j++) {
      /* foo_ij: anchor #i in the component, species #j */
      comp_ij = comp[j+1][w];

      base_field = j*fields_per_species + FIELD_base;
      start_ij = anc_i[base_field + FIELDR_start];
      len_ij = anc_i[base_field + FIELDR_len];
      end_ij = start_ij + len_ij;
#if 0
      sign_ij = anc_i[base_field + FIELDR_sign];
#endif

      /* If it's root of a component, set up a new interval.
       * Else, expand the start/end for existing interval.
       */
      if (comp_ij == w) {
	/* This node is the root of its component in this species */
	if (repeat_params->intervals_len >= intervals_alloc) {
	  /* need to allocate more space to intervals table */
	  intervals =
	    int_array2d_resize(intervals,
			       intervals_alloc,
			       RP_size,
			       intervals_alloc + REPEAT_INC);
	  repeat_params->intervals = intervals;
	  intervals_alloc = 
	    repeat_params->intervals_alloc =
	    intervals_alloc + REPEAT_INC;
	}

	cur_int = intervals[repeat_params->intervals_len++];

	/* component and species identify this interval */
	cur_int[RP_component] = comp_ij;
	cur_int[RP_species] = j;

	chr_ij = anc_i[base_field + FIELDR_chr];
	cur_int[RP_chr] = chr_ij;

	/* initial values for start, end; will be computed this loop  */
	cur_int[RP_start] = start_ij;
	cur_int[RP_end] = end_ij;

	/* signs will be computed later */
	cur_int[RP_sign] = 0;
	cur_int[RP_ct_p] = 0;
	cur_int[RP_ct_m] = 0;

      } else {
	/* Not the root of the component, just grow [start,end] */
	cur_int = get_interval(repeat_params, comp_ij, j);

	if (cur_int[RP_start] > start_ij) {
	  cur_int[RP_start] = start_ij;
	}
	if (cur_int[RP_end] < end_ij) {
	  cur_int[RP_end] = end_ij;
	}
      }
    }
  }

  /****************************************************************
   * Determine signs of intervals
   ****************************************************************/

  ga_repeat_signs(repeat_params);

  /****************************************************************
   * Output signed intervals in each species
   ****************************************************************/

  /* print header:
   *     family 30: 10 x 4 (25)
   *     "family 30"=ID number
   *     10:    # repeats in organism 1
   *     4:     # repeats in organism 1
   *     (25):  # supporting anchors
   */

  fprintf(repeat_out,"repeat %d:", repeat_params->repeat_id++);
  /* print dimensions in each species */
  for (j=0; j<nspecies_s; j++) {
    count = 0;
    for (m=0 ; m < repeat_params->intervals_len ; m++) {
      cur_int = intervals[m];
      if (cur_int[RP_species] == j) count++;
    }
    if (j==0) {
      fprintf(repeat_out," %d", count);

    } else {
      fprintf(repeat_out," x %d", count);
    }
  }
  fprintf(repeat_out," (%d)", same_species ? (len/nspecies) : len);
  if (merged) {
    fprintf(repeat_out," merged");
  }
  fprintf(repeat_out,"\n");

  /* print intervals per-species */
  for (j=0; j<nspecies_s; j++) {
    fprintf(repeat_out,"species %d:\n",j+1);
    for (m=0 ; m < repeat_params->intervals_len ; m++) {
      cur_int = intervals[m];
      if (cur_int[RP_species] != j) continue;

      scratch_anc[0] = cur_int[RP_chr];                     /* chromosome */
      scratch_anc[1] = cur_int[RP_start];                   /* start */
      scratch_anc[2] = cur_int[RP_end] - cur_int[RP_start]; /* length */
      scratch_anc[3] = cur_int[RP_sign];                    /* sign */

      print_species_coords(gs_params,
			   repeat_out,
			   anctable,
			   scratch_anc,
			   j);

      if (cur_int[RP_ct_p]>0 && cur_int[RP_ct_m]>0) {
	fprintf(repeat_out," # ambiguous sign");
      } else if (cur_int[RP_ct_p]==0 && cur_int[RP_ct_m]==0) {
	fprintf(repeat_out," # ambiguous sign ");
      }

      fprintf(repeat_out,"\n");
    }
  }


#if 0
  /* DEBUGGING */
  fprintf(repeat_out,"DEBUG INTERVALS\n"); fflush(repeat_out);
  for (m=0 ; m < repeat_params->intervals_len ; m++) {
    cur_int = intervals[m];
    fprintf(repeat_out,"INTERVAL %d:", m);
    for (j=0; j<RP_size; j++) {
      fprintf(repeat_out," %d",cur_int[j]);
    }
    fprintf(repeat_out,"\n");
  }
#endif

  /****************************************************************
   * Print out table of alignment coordinates in repeat family
   ****************************************************************/

  fprintf(repeat_out,"support:\n",j);
  for (i=0; i<len; i++) {
    w = comp_nodes[i];

    /* in same_species mode,
     * only print out original coordinates, not shifts
     */
    if (same_species && ((nodelist[w] % nspecies) > 0)) continue;

    anc_i = anclist[nodelist[w]];
    anc_sign = interval_sign(repeat_params, anc_i, w, 1);
    if (anc_sign > 1) {
      ambig = 1;
      anc_sign -= 3;
    } else {
      ambig = 0;
    }

    /* copy original anchor */
    memcpy((void *) scratch_anc,
	   (const void *) anc_i,
	   bytes_per_anchor);

    /* modify the signs */
    if (anc_sign == -1) {
      for (j=0; j<nspecies; j++) {
	base_field = j*fields_per_species + FIELD_base;
	sign_ij = scratch_anc[base_field + FIELDR_sign];
	scratch_anc[base_field + FIELDR_sign] = -sign_ij;
      }
    }

#if 0
    fprintf(repeat_out,"[w=%d,in=%d] ", w, nodelist[w]);
    for (j=-1; j<nspecies; j++) {
      /* foo_ij: anchor #i in the component, species #j */
      comp_ij = comp[j+1][w];
#if 0
      base_field = j*fields_per_species + FIELD_base;
      start_ij = anc_i[base_field + FIELDR_start];
      len_ij = anc_i[base_field + FIELDR_len];
      end_ij = start_ij + len_ij;
      sign_ij = anc_i[base_field + FIELDR_sign];
#endif
      fprintf(repeat_out,"(%d) ", comp_ij);
    }
#endif

    /* ID # from orig anchor */
    fprintf(repeat_out,"%d",scratch_anc[FIELD_id]);

    print_block_nocr(gs_params, repeat_out, anctable, scratch_anc);
    if (ambig) {
      fprintf(repeat_out, " # sign error");
    }
    fprintf(repeat_out,"\n");
  }
  fprintf(repeat_out,"\n");
}

/*
 * Assign all intervals signs +1,0,-1 as appropriate.
 * If there are no inconsistencies, it will get the unique solution.
 * If there are inconsistencies, it will make up to REPEAT_SIGN_TRIES
 * attempts, taking the one with the most assigned +/-1 signs,
 * where ambiguous signs are determined with the sign_threshold parameter.
 */

void ga_repeat_signs(REPEAT_PARAMS_t *repeat_params)
{
  int **intervals = repeat_params->intervals;
  int intervals_len = repeat_params->intervals_len;

  int m;
  int nsigns_best=0, nsigns_best_m=0, last_m=-1;
  int nsigns_cur;

  int ntries = REPEAT_SIGN_TRIES;

  /* Assign some interval a sign, then use alignments its in to
   * determine the signs of other intervals, and iterate till get
   * as many signs as possible.
   * If there are no sign inconsistencies, it will not matter which
   * one we start with.
   * If there are sign inconsistencies, it will depend on which one
   * we start with.
   */

  for (m=0; m<intervals_len && ntries>0; m++) {
    /* Can seed with m=0, or with next m with undetermined sign on previous try */
    if (m>0 && intervals[m][RP_sign] != 0) continue;

    last_m = m;
    nsigns_cur = ga_repeat_signs_interval(repeat_params,m);

    if (2*nsigns_cur >= intervals_len) {
      /* For perfectly consistent signs, we actually have nsigns_cur == intervals_len
       * For inconsistent signs, stop if we get at least half of them
       * (or if we do too many tries)
       */
      return;
    }
    if (nsigns_cur > nsigns_best) {
      /* got more signs than previously */
      nsigns_best = nsigns_cur;
      nsigns_best_m = m;
    }

    ntries--;
  }

  /* if got here, then never got a perfect signage, so use the best one */

  if (nsigns_best_m != last_m)
    ga_repeat_signs_interval(repeat_params, nsigns_best_m);
}


/*
 * Determine all interval signs.
 * Start with interval # seed_m with sign +1.
 * Propagate it to other intervals via intersecting anchors.
 *
 * Return:
 *    Signs +1,0,-1 are set in the intervals.
 *    return value = # intervals assigned +/-1 signs
 */
int ga_repeat_signs_interval(REPEAT_PARAMS_t *repeat_params,
			     int seed_m)
{
  ANCTABLE *anctable = repeat_params->anctable;
  GANCTREE *ganctree = repeat_params->ganctree;
  int **intervals = repeat_params->intervals;
  int intervals_len = repeat_params->intervals_len;
  int *comp_nodes = repeat_params->comp_nodes;
  int len = repeat_params->comp_nodes_len;
  int same_species = repeat_params->gs_params->same_species;

  int nspecies = anctable->nspecies;
  int **anclist = anctable->anclist;
  int fields_per_species = anctable->fields_per_species;

  int *nodelist = ganctree->nodelist;
  int **comp = ganctree->comp;

  int i,j,m;
  int *anc_i;
  int base_field;
  int sign_ij, comp_ij;
  int w;
  int *cur_int;

  int count_signs, count_signs_last;
  int anc_sign, new_sign;
  int count_p, count_m;

  double sign_threshold = repeat_params->gs_params->sign_threshold;
  double sign_frac;

  /****************************************************************
   * Set all signs to unknown and all counts to 0
   ****************************************************************/

  for (m=0; m<intervals_len; m++) {
    cur_int = intervals[m];
    cur_int[RP_sign] = 0;
    cur_int[RP_ct_p] = 0;
    cur_int[RP_ct_m] = 0;
  }

  /* seed sign set to +1 */
  intervals[seed_m][RP_sign] = 1;

  /****************************************************************
   * Successive passes will use signs from previous passes
   ****************************************************************/
  count_signs_last = -1;
  count_signs = 0;

  while (count_signs < intervals_len
	 && count_signs > count_signs_last) {

    count_signs_last = count_signs;

    /****************************************************************
     * Loop over alignments.
     * For each alignment, determine its overall sign if possible,
     * based on consistency with majority of putative signs.
     * If that is possible, then for intervals with unknown signs, 
     * vote for their signs.
     ****************************************************************/

    for (i=0; i<len; i++) {
      w = comp_nodes[i];

      /* in same species mode, skip the mirrored nodes */
      if (same_species && (nodelist[w] % nspecies))
	continue;

      anc_i = anclist[nodelist[w]];
      anc_sign = interval_sign(repeat_params, anc_i, w, 0);
      if (anc_sign == 0) {
	/* Alignment did not intersect with previous intervals with
	 * putative signs, or else they tied
	 */
	continue;
      }

      for (j=0; j<nspecies; j++) {
	/* foo_ij: anchor #i in the component, species #j */
	comp_ij = comp[j+1][w];

	base_field = j*fields_per_species + FIELD_base;
	sign_ij = anc_i[base_field + FIELDR_sign] * anc_sign;

	/* In case anchor is unsigned */
	if (sign_ij == 0) continue;

	cur_int = get_interval(repeat_params, comp_ij, j);

	/* vote for sign */
	cur_int[(sign_ij == 1) ? RP_ct_p : RP_ct_m]++;
      }
    }

    /****************************************************************
     * Compute signs from votes.
     ****************************************************************/

    count_signs = 0;
    for (m=0; m<intervals_len; m++) {
      cur_int = intervals[m];
      count_p = cur_int[RP_ct_p];
      count_m = cur_int[RP_ct_m];

      /* (a) If both are 0, then no decision on sign.
       * (b) If one is 0 and one is nonzero, the nonzero one is the winner.
       * (c) Else, there are sign inconsistencies, and should use a threshold
       * to determine if a vote is cast or if it is left ambiguous.
       *
       * Test below includes (a) & (b),
       * and does (c) with threshold > 2/3 of total votes
       */

#if 0
      if (count_p > 2*count_m) {
	new_sign = 1;
	count_signs++;
      } else if (count_m > 2*count_p) {
	new_sign = -1;
	count_signs++;
      } else {
	new_sign = 0;
      }
#endif
      /* do (a),(b),(c) explicitly instead of combined */
      if (count_p == 0) {
	if (count_m > 0) {
	  new_sign = -1; count_signs++;    /* (b) with winner -1 */
	} else {
	  new_sign = 0;                    /* (a) */
	}
      } else {
	if (count_m == 0) {
	  new_sign = 1; count_signs++;     /* (b) with winner +1 */
	} else {
	  /* both count_p, count_m are nonzero, sign is ambiguous */
	  sign_frac = (double) count_p / (double) (count_m+count_p);

	  if (sign_frac > sign_threshold) {
	    new_sign = 1; count_signs++;   /* (c) with winner +1 */
	  } else if (1-sign_frac > sign_threshold) {
	    new_sign = -1; count_signs++;  /* (c) with winner -1 */
	  } else {
	    new_sign = 0;
	  }
	}
      }

      cur_int[RP_sign] = new_sign;
    }
  }

  /****************************************************************
   * # intervals for which signs were determined
   ****************************************************************/

  return count_signs;
}



/*
 * return values if flag_ambiguous=0:
 *      1 or -1: majority vote for that sign
 *      0: tied vote
 *
 * return values if flag_ambiguous=1:
 *      1 or -1: >0 votes for that sign, none for the opposite sign
 *      0: no votes either direction
 *   add 3 for ambiguities:
 *      2: sign ambiguously -1
 *      3: tied with >0 votes both ways
 *      4: sign ambiguously 1
 */

int interval_sign(REPEAT_PARAMS_t *repeat_params,
		  int *anc_i,
		  int w,
		  int flag_ambiguous)
{
  ANCTABLE *anctable = repeat_params->anctable;
  int nspecies = anctable->nspecies;
  int fields_per_species = anctable->fields_per_species;
  int **intervals = repeat_params->intervals;
  int intervals_len = repeat_params->intervals_len;
  int *cur_int;
  int same_species = repeat_params->gs_params->same_species;

  GANCTREE *ganctree = repeat_params->ganctree;
  int **comp = ganctree->comp;

  int base_field;
  int comp_ij;
  int sign_compare;
  int j;
  int s;



  /* see how many putative interval signs agree (_p) or disagree (_m)
   * with the signs in this alignment */
  int count_p = 0;
  int count_m = 0;

  for (j=0; j<nspecies; j++) {
    base_field = j*fields_per_species + FIELD_base;

    comp_ij = comp[j+1][w];
    cur_int = get_interval(repeat_params, comp_ij, j);

    sign_compare = anc_i[base_field + FIELDR_sign] * cur_int[RP_sign];

    if (sign_compare > 0) {
      count_p++;
    } else if (sign_compare < 0) {
      count_m++;
    }
  }

  s = (count_p > count_m) ? 1 : (count_p < count_m) ? -1 : 0;
  if (flag_ambiguous && count_p > 0 && count_m > 0) { s += 3; }
  return s;
}

/****************************************************************************/

void debug_interval_table(REPEAT_PARAMS_t *repeat_params)
{
  FILE *repeat_out = repeat_params->repeat_out;
  int **intervals = repeat_params->intervals;
  int intervals_len = repeat_params->intervals_len;
  int m, i;
  int *cur_int;

  fprintf(repeat_out, "INTERVALS table (%d)\n", intervals_len);
  for (m=0; m<intervals_len; m++) {
    cur_int = intervals[m];
    fprintf(repeat_out, "[%d]:", m);
    for (i=0; i<RP_size; i++) {
      fprintf(repeat_out, " %d", cur_int[i]);
    }
    fprintf(repeat_out, "\n");
  }
}



/* (int **) a,b;   are cast to (const void *) */
int interval_cmp(const void *a, const void *b)
{
  int *interval_a = *(int **)a;
  int *interval_b = *(int **)b;
  int d = interval_a[RP_component] - interval_b[RP_component];

  if (d == 0) {
    d = interval_a[RP_species] - interval_b[RP_species];
  }
  return d;
}




int *get_interval(REPEAT_PARAMS_t *repeat_params,
		  int component,
		  int species)
{
  int key[RP_size], *key_ptr=key;
  int **intervals_ptr;

  if (species > 0 && repeat_params->gs_params->same_species) {
    /* get equivalent component in species 0 */
    key[RP_component] =
      same_species_mate_comp(repeat_params, component, species);
    key[RP_species] = 0;
  } else {
    key[RP_component] = component;
    key[RP_species] = species;
  }

  intervals_ptr = 
    bsearch((const void *) &key_ptr,
	    (const void *) repeat_params->intervals,
	    repeat_params->intervals_len,
	    sizeof(&key),
	    interval_cmp);

  if (intervals_ptr == (int **) NULL) {
    fprintf(stderr, "get_interval failed! comp=%d, species=%d\n", component, species);
    exit(-1);
  }
  return *intervals_ptr;
}


/* For same-species mode, each input coordinate line is duplicated/shifted:
 *     ID c0 c1 c2
 *     ID c2 c0 c1
 *     ID c1 c2 c0
 * in groups of nspecies lines.
 * Given anchor # nodelist[v], species s,
 *   same_species_mate(...,v,s) = w
 *      where anchor # nodelist[w], species 0, corresponds to it
 *   same_species_mate_comp(...,v,s) = w2
 *      where w2 is the species 0 component index of w
 */
int same_species_mate(REPEAT_PARAMS_t *repeat_params,
		      int v, int s)
{
  ANCTABLE *anctable = repeat_params->anctable;
  int **anclist = anctable->anclist;
  GANCTREE *ganctree = repeat_params->ganctree;
  int *nodelist = ganctree->nodelist;
  int nspecies = anctable->nspecies;

  int anc_no = nodelist[v];   /* actual input anchor number */
  int anc_mate;               /* input anchor number of mate */
  int v_mate;                 /* node # of mate */

  /* for nspecies == 2:
   *   anc_mate = anc_no ^ s
   * In general: write anc_no in base nspecies
   *             subtract s from 1's place w/o borrowing.
   */
  anc_mate = (anc_no % nspecies) - s;
  anc_mate = (anc_mate < 0) ? anc_no-s + nspecies : anc_no-s;

  v_mate = anclist[anc_mate][FIELD_mate];

  return v_mate;
}

int same_species_mate_comp(REPEAT_PARAMS_t *repeat_params,
			   int v, int s)
{
  ANCTABLE *anctable = repeat_params->anctable;
  int **anclist = anctable->anclist;
  GANCTREE *ganctree = repeat_params->ganctree;
  int *nodelist = ganctree->nodelist;
  int nspecies = anctable->nspecies;

  int anc_no = nodelist[v];   /* actual input anchor number */
  int anc_mate;               /* input anchor number of mate */
  int v_mate;                 /* node # of mate */

  /* for nspecies == 2:
   *   anc_mate = anc_no ^ s
   * In general: write anc_no in base nspecies
   *             subtract s from 1's place w/o borrowing.
   */
  anc_mate = (anc_no % nspecies) - s;
  anc_mate = (anc_mate < 0) ? anc_no-s + nspecies : anc_no-s;

  v_mate = anclist[anc_mate][FIELD_mate];

  return repeat_params->ganctree->comp[0+1][v_mate];
}


