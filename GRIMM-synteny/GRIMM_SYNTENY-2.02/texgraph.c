/* texgraph.c
 *    Output a graph drawing for use with LaTeX package bpgraph.sty.
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

/* Last modified on Wed Aug 2, 2006, by Glenn Tesler
 */

/* Glenn Tesler */
/* Output a graph drawing for use with LaTeX package bpgraph.sty */
/* formats:
 * lr_or_cc:
 *        TEXGR_LR: list vertices & edges from left to right
 *        TEXGR_CC: organize by connected components, and then by
 *                  paths/cycles within connected components
 *
 * stage:
 *        1: breakpoint graph G(pi^,gamma^)
 *        2: no tails         G(Pi^,Gamma^)
 *        3: delete edges     G(Pi,Gamma)
 *        4: add edges        G(Pi*,Gamma*)
 *        5: optimal concat   G(Pi**,Gamma**)
 *        6: add tails        G(pi*,gamma*)
 *
 * showlabels:
 *        0: no labels
 *        1: signed entries
 *        2: unsigned vertex #s
 *        3: both
 */

#include <stdio.h>
#include <stdlib.h>
#include "uniinvdist.h"
#include "mcstructs.h"
#include "mcrdist.h"
#include "graph_edit.h"
#include "graph_components.h"
#include "texgraph.h"

int seek_cc(int cur_cc, int gsize, distmem_t *distmem);

int seek_cyc(int cur_cyc, int compno, int gsize, distmem_t *distmem);

void texgraph_drawcyc(int *pilabels, int *galabels,
		      int num_genes, int num_chromosomes,
		      int showtails, int showlabels,
		      distmem_t *distmem,
		      int cyc_index,
		      int lr_or_cc);

void texgraph_cc(int *pilabels, int *galabels,
		 int num_genes, int num_chromosomes,
		 int showtails,
		 int showlabels,
		 distmem_t *distmem,
		 int lr_or_cc);

void texgraph_lr(int *pilabels, int *galabels,
		 int num_genes, int num_chromosomes,
		 int showtails,
		 int showlabels,
		 distmem_t *distmem);

void texgraph_brackets(int *pilabels, int num_genes, int num_chromosomes);

void texgraph_pis(int *pilabels, int num_genes, int num_chromosomes);

int texgraph_cyclen(distmem_t *distmem,
		    int start,
		    int num_chromosomes,

		    int *n_ori,
		    int *n_chrom,
		    int *n_ichrom);

void print_hurdle_spec(int h);


/* options: string with a combination of these letters
 * l: organize by left to right
 * c: organize by cycle/component
 * (use l or c, not both)
 *
 * b: show chromosome brackets
 * u: show unsigned labels
 * s: show signed labels
 *
 * 1-6: stage
 *
 * default: stage=1; organize=l; do not show brackets or (un)signed labels
*/

void draw_texgraph(struct genome_struct *g1,
		   struct genome_struct *g2,
		   int num_genes, int num_chromosomes,
		   char *options)
{
  int stage = 1;
  int lr_or_cc = TEXGR_LR;
  int showlabels = 0;

  char *s = options;
  int c;

  while ((c=*s) != '\0') {
    switch (c) {
    case 'l': case 'L':
      lr_or_cc = TEXGR_LR; break;

    case 'c': case 'C':
      lr_or_cc = TEXGR_CC; break;

    case 'r': case 'R':
      /* "Raw" cycle/path information -- not really TeX */
      /* Honors "a" and stage, ignores b,u,s */
      lr_or_cc = TEXGR_RAW; break;


    case 'b': case 'B':
      showlabels |= TEXGR_B; break;   /* chromosome brackets */

    case 'u': case 'U':
      showlabels |= TEXGR_U; break;   /* "unsigned" doubling of genes */

    case 'a': case 'A':
      /* alternate to unsigned -- show AB format of doubled genes instead */
      showlabels |= TEXGR_U | TEXGR_AB; break;

    case 's': case 'S':
      showlabels |= TEXGR_S; break;   /* signed (orig) entries of pi */

    case '1': case '2': case '3': case '4': case '5': case '6':
      stage = c - '0'; break;



    default:
      break;
      /* TODO: ERROR MESSAGE */
    }

    s++;
  }

  texgraph_nomem(g1,g2,num_genes,num_chromosomes,
		 lr_or_cc, showlabels, stage);
}


void texgraph_nomem(struct genome_struct *g1,
		   struct genome_struct *g2,
		   int num_genes, int num_chromosomes,
		   int lr_or_cc,
		   int showlabels,
		   int stage)
{
  distmem_t distmem;

  /* allocate memory */
  mcdist_allocmem(num_genes, num_chromosomes,
		  &distmem);

  texgraph(g1,g2,num_genes,num_chromosomes,
	   lr_or_cc,showlabels,stage,
	   &distmem);
}

void texgraph(struct genome_struct *g1,
	     struct genome_struct *g2,
	     int num_genes, int num_chromosomes,
	     int lr_or_cc,
	     int showlabels,
	     int stage,
	     distmem_t *distmem)
{
  graph_t G;
  pathcounts_t pcounts_G;

  /* draw tail cycles or not? */
  int showtails = 0;

  /* pi, gamma as signed permutations */
  int *pilabels = g2->genes;
  int *galabels = g1->genes;

  struct genome_struct new_g1, new_g2;


  if (num_chromosomes == 0) {
    stage = 1;
  }

  if (stage<1 || stage>6)
    stage = 1;

  switch (stage) {
  case 1:    /* 1: breakpoint graph G(pi^,gamma^) */
    showtails = 1;
    invdist_noncircular(g1,g2,0,num_genes,distmem);
    break;

  case 2:    /* 2: no tails         G(Pi^,Gamma^) */
    showtails = 0;
    invdist_noncircular(g1,g2,0,num_genes,distmem);
    break;

  case 3:    /* 3: delete edges     G(Pi,Gamma) */
    showtails = 0;
    mcdist_noncircular_graph_d(g1,g2,num_genes,num_chromosomes,distmem);
    break;

  case 4:    /* 4: add edges        G(Pi*,Gamma*) */
    showtails = 0;
    mcdist_capgraph_connections(g1,g2,num_genes,num_chromosomes,
				distmem,
				&G, &pcounts_G,
				0);
    break;

  case 5:    /* 5: optimal concat   G(Pi**,Gamma**) */
  case 6:    /* 6: add tails        G(pi*,gamma*) */
    showtails = (stage == 5) ? 0 : 1;

    mcdist_capgraph(g1,g2,num_genes,num_chromosomes,distmem,0);
    pilabels = distmem->cappedp2;
    galabels = distmem->cappedp1;

    if (stage == 6) {
      new_g1.genes = galabels;
      new_g2.genes = pilabels;

      /* add in the tails by rebuilding graph */
      /* TODO: just add the necessary edges to the current graph */
      showtails = 1;
      invdist_noncircular(&new_g1,&new_g2,0,num_genes,distmem);
      break;

    }

    break;
  }

  fprintf(outfile,"%%\n%%\n\\begin{bpgraph}\n");

  if (lr_or_cc == TEXGR_RAW) {
    texgraph_cc(
		pilabels, galabels,
		num_genes, num_chromosomes,
		showtails,
		showlabels,
		distmem,
		lr_or_cc);
    return;
  }

  if (lr_or_cc == TEXGR_LR) {
    texgraph_lr(pilabels, galabels,
		num_genes, num_chromosomes,
		showtails,
		showlabels,
		distmem);
  } else {
    texgraph_cc(
		pilabels, galabels,
		num_genes, num_chromosomes,
		showtails,
		showlabels,
		distmem,
		lr_or_cc);
  }

  /* draw brackets underneath chromosomes */
  if (num_chromosomes > 0 &&
      (showlabels & TEXGR_B))
    texgraph_brackets(pilabels, num_genes, num_chromosomes);

  fprintf(outfile,"%%\n%%\n\\end{bpgraph}\n");

}


int vlabel_u(int *pilabels, int gsize, int vnum)
{
  int pilab;

  if (vnum == 0 || vnum == gsize-1)
    pilab = vnum;
  else {
    pilab = pilabels[(vnum-1)/2];

    if (vnum % 2 == 0)
      pilab = -pilab;

    if (pilab < 0)
      pilab = -2*pilab;
    else
      pilab = 2*pilab - 1;
  }

  return pilab;
}


char *vlabel_ab(int *pilabels, int gsize, int vnum)
{
  int pilab;
  pilab = vlabel_u(pilabels, gsize, vnum);
  return (pilab % 2 == 0) ? "b" : "a";
}

char *vlabel(int *pilabels, int gsize, int vnum)
{
  static char vmem[10];
  int pilab = vlabel_u(pilabels,gsize,vnum);
  int n = (pilab+1)/2;
  int ab = (pilab+1) % 2;

  sprintf(vmem,"%d%c", n, "ab"[ab]);

  return vmem;
}

int showv(int *pilabels, int gsize, int vnum, int lowcap, int showtails)
{
  int lowcap_u, pi_u;

  if (showtails)
    return 1;

  lowcap_u = 2*lowcap - 1;

  pi_u = vlabel_u(pilabels, gsize, vnum);

  if (pi_u == 0)
    return 0;

  if (pi_u < lowcap_u)
    return 1;

  pi_u -= lowcap_u;
  pi_u &= 3;

  if (pi_u == 1 || pi_u == 2)
    return 1;

  return 0;
}


/* draw graph in left-to-right fashion
 * pilabels,galabels = pi, gamma as signed perms
 * showtails=0: don't show tails; 1: do
 * distmem: has a representation of the graph
 */


/* #define SHOWV(i) (showtails || 
 *                  (i>0 && i<gsize-1 && 
 *		   pilabels[(i-1)/2]<lowcap && pilabels[(i-1)/2]>-lowcap))
 */

#define SHOWV(i) showv(pilabels, gsize, i, lowcap, showtails)

void texgraph_lr(int *pilabels, int *galabels,
		 int num_genes, int num_chromosomes,
		 int showtails,
		 int showlabels,
		 distmem_t *distmem)
{
  int *greyEdges  = distmem->greyEdges;


  int lowcap = num_genes - 2*num_chromosomes + 1;
  int gsize = 2*num_genes + 2;
  int countout;

  int i=0,j=0;

  int pilab, s;



  /* draw vertices */
  fprintf(outfile,"%% vertices\n");

  for (i=0; i<gsize; i++) {
    if (i % 4 != 0)
      fprintf(outfile," ");

    if (SHOWV(i)) {
      fprintf(outfile,"\\vert{%d}{%s}", i, vlabel(pilabels,gsize,i));
    } else {
      fprintf(outfile,"\\vert[\\pvert]{%d}{%s}", i, vlabel(pilabels,gsize,i));
    }

    if (i % 4 == 0) {
      fprintf(outfile,"%%\n");
    }
  }

  fprintf(outfile,"%%\n%%\n");


  /* draw black edges */
  fprintf(outfile,"%% black edges\n");

  countout = 0;
  for (i=0; i<gsize; i += 2) {
    if (SHOWV(i)) {
      if (countout % 4 != 0)
	fprintf(outfile," ");

      fprintf(outfile,"\\bedge{%s}", vlabel(pilabels,gsize,i));
      fprintf(outfile,"{%s}", vlabel(pilabels,gsize,i+1));

      countout++;
    } else {
    }

    if (countout % 4 == 0) {
      fprintf(outfile,"%%\n");
    }
  }

  fprintf(outfile,"%%\n%%\n");

  /* draw gray edges
   * draw from the left vertex to the right vertex
   */
  fprintf(outfile,"%% gray edges\n");

  countout = 0;
  for (i=0; i<gsize; i++) {
    if (!SHOWV(i)) continue;

    j = greyEdges[i];

    if (j<0) {
      /* no edge to show if at end of a path */
      if (j == V_PICAP || j == V_GTAIL)
	continue;

      /* adjacency or tail, follow edge */
      if (i % 2 == 0)
	j = i+1;
      else
	j = i-1;
    }

    /* only show edges left to right */
    if (i>j) continue;

    if (countout % 4 > 0)
      fprintf(outfile," ");

    fprintf(outfile,"\\gedge{%s}", vlabel(pilabels, gsize, i));
    fprintf(outfile,"{%s}", vlabel(pilabels, gsize, j));

    countout++;

    if (countout % 4 == 0)
      fprintf(outfile,"%%\n");
  }

  fprintf(outfile,"%%\n%%\n");


  /* mark tails, pi-caps, gamma-tails */
  if (num_chromosomes > 0) {
    fprintf(outfile,"%% mark tails, pi-caps, gamma-tails\n");
    for (i=0; i<gsize; i++) {
      j = greyEdges[i];
      if (j>=0) continue;

      switch (j) {
      case V_ADJ:
	/* tails may still be adjacencies in some stages */
	if (i>0 && i<gsize-1) {
	  s = 1;
	  pilab = pilabels[(i-1)/2];
	  if ((i & 1) == 0)
	    s = -s;
	  if (pilab < 0) {
	    pilab = -pilab;
	    s = -s;
	  }

	  if (pilab < lowcap)
	    break;

	  pilab -= lowcap;
	  if ((pilab & 1) == 1)
	    s = -s;

	  if (s == -1)
	    break;
	}
      case V_TAIL:
	if (!showtails)
	  fprintf(outfile,"\\marktail{%s}%%\n", vlabel(pilabels,gsize,i));
	break;
      case V_PICAP:
	fprintf(outfile,"\\markpi{%s}%%\n", vlabel(pilabels,gsize,i));
	break;
      case V_GTAIL:
	fprintf(outfile,"\\markga{%s}%%\n", vlabel(pilabels,gsize,i));
	break;
      }
    }
  }

  fprintf(outfile,"%%\n%%\n");
  

  /* signed entries of pi */
  if (showlabels & TEXGR_S)
    texgraph_pis(pilabels, num_genes, num_chromosomes);


  /* unsigned entries of pi */
  if (showlabels & TEXGR_U) {
    fprintf(outfile,"%% pi, doubled entries\n");

    countout = 0;
    for (i=0; i<gsize; i++) {

      if (countout % 3 != 0)
	fprintf(outfile, " ");

      if (showlabels & TEXGR_AB) {
	fprintf(outfile,
		"\\showpiu{%s}{%s}",
		vlabel(pilabels, gsize, i),
		vlabel_ab(pilabels, gsize, i));
      } else {
	fprintf(outfile,
		"\\showpiu{%s}{%d}",
		vlabel(pilabels, gsize, i),
		vlabel_u(pilabels, gsize, i));
      }

      countout++;

      if (countout % 3 == 0)
	fprintf(outfile,"%%\n");
    }
  }
  fprintf(outfile,"%%\n%%\n");
}


/* draw the graph, ordered by connected components.
 * Within connected component, do each path/cycle separately.
 * Do each path/cycle by following from one end to another.
 */

#define DRAWCYC(cyc) \
    texgraph_drawcyc(pilabels, galabels, num_genes, num_chromosomes, \
                     showtails, showlabels, distmem, cyc, lr_or_cc)

void texgraph_cc(int *pilabels, int *galabels,
		 int num_genes, int num_chromosomes,
		 int showtails,
		 int showlabels,
		 distmem_t *distmem,
		 int lr_or_cc)
{
  int *cc    = distmem->cc;                /* con. comp. # of each vertex */

  /* flags on each component */
  component_t *components = distmem->components;
  component_t *comp;

  int gsize = 2*num_genes + 2;

  int cur_cc = -2;
  int cur_cyc;

  int has_sgrk;

  graph_t G0,   *G = &G0;
  pathcounts_t  pcounts_G0,   *pcounts_G = &pcounts_G0;

  hurdlecounts_t   hurdlesIU_G,  hurdlesRU_G;

  int j;
  int *greyEdges = distmem->greyEdges;


  /* build the paths and cycles */
  G->distmem = distmem;
  G->size = 2*num_genes+2;
  G->components_good = FALSE;

  repairG_chrom(num_genes, num_chromosomes, distmem, G);
  classify_cycle_path(G, pcounts_G);
  form_connected_components(G);
  classify_connected_components(G);

  if (lr_or_cc == TEXGR_RAW) {
    if (num_chromosomes == 0) {
      /* unichromosomal */
      calc_num_hurdles_and_fortress(G, C_U, &hurdlesIU_G);
    } else {
      /* multichromosomal */
      /* calculate # semiknots and other hurdle info */

      /* IU no longer needed, but do it anyway for annotation*/
      calc_num_hurdles_and_fortress(G, C_IU, &hurdlesIU_G);

      calc_num_hurdles_and_fortress(G, C_RU, &hurdlesRU_G);

      /* mark the components that are semi-knots or simple */
#if 0
      calc_num_semiknots(G);
#endif
      calc_num_semirealknots(G,&has_sgrk);

      /* only done for annotation */
      calc_num_simple(G);
    }
  }



  /* go in increasing order of left vertex of components */
  while ((cur_cc = seek_cc(cur_cc+2, gsize, distmem)) >= 0) {
    if (lr_or_cc == TEXGR_CC) {
      fprintf(outfile,"%% Component %d\n", cur_cc);
    } else {
      /* Raw format */
      fprintf(outfile,"%% Component @ %d:", cur_cc);

      /* flags about its type */
      if (cc[cur_cc] >= 0) {
	comp = &components[cc[cur_cc]];
	fprintf(outfile, "%s",
		comp->oriented ? " o" : " u");

	if (num_chromosomes > 0) {
	  fprintf(outfile, "%s",
		  comp->intrachrom ? " intra" : " inter");
	  fprintf(outfile, "%s",
		  /*		  comp->real ? " r" : " nr" */
		  comp->realf ? " nr" : " r"
		  );
	}


#if 0
	/* this is reset by each call to calc_num_hurdles_and_fortress,
	 * its value at this point is useless
	 */
	fprintf(outfile, " b:%d", comp->blocks);
#endif

	
	/* hurdlemania */

	/*	fprintf(outfile, " h:%o", comp->hurdle);
	 */
	fprintf(outfile, " h:");
	print_hurdle_spec(comp->hurdle);
      } else {
	/* a trivial component */
	j = greyEdges[cur_cc];
	if (j == V_TAIL) {
	  /* tail */
	  fprintf(outfile, " t");
	} else if (j == V_ADJ) {
	  /* adjacency */
	  fprintf(outfile, " a");
	} else {
	  /* bare edge */
	  fprintf(outfile, " b");
	}
      }
      fprintf(outfile, "\n");
    }

    /* within component, do cycles/paths in increasing order of
     * leftmost vertex */
    cur_cyc = cur_cc;
    DRAWCYC(cur_cyc);

    /* if adjacency, tail, bare edge: component is over.
     * else, find next path/cycle in it. */
    if (cc[cur_cc] >= 0)
      while ((cur_cyc = seek_cyc(cur_cyc+2, cur_cc, gsize, distmem)) >= 0) {
	DRAWCYC(cur_cyc);
      }
  }


  /* signed pi labels */
  if ((lr_or_cc == TEXGR_CC) && (showlabels & TEXGR_S))
    texgraph_pis(pilabels, num_genes, num_chromosomes);


}

void print_hurdle_spec(int h)
{
  /* codes for the hurdle flags in the bits as given in mcstructs.h */ 
  const char *specs[] = { "h","gh","sh",
			  "k","gk","sk",
			  "rk","grk","srk",
			  "sem","sim" };
  int nprops=0;
  int i;
  for (i=0;   
       i<11  &&  h != 0;
       i++) {
    if (h & 1) {
      if (nprops > 0) { fprintf(outfile,","); }
      fprintf(outfile, "%s", specs[i]);
      nprops++;
    }
    h >>= 1;
  }
  if (nprops == 0) { fprintf(outfile,"-"); }
}


/* Find first component starting at cur_cc or later */
int seek_cc(int cur_cc, int gsize, distmem_t *distmem)
{
  int *cc         = distmem->cc;        /* connected component w/vertex */
  component_t *components = distmem->components; /* list of components */
  int cc_ind;
  int i;

  for (i=cur_cc; i<gsize; i += 2) {
    cc_ind = cc[i];

    if (cc_ind < 0 ||       /* adjacency, tail, bare edge */
	i == components[cc_ind].index
                            /* leftmost vertex of normal component */
	)
      return i;
  }

  /* no more components found */
  return -1;
}

/* Find path or cycle starting at cur_cc or later */
int seek_cyc(int cur_cyc, int compno, int gsize, distmem_t *distmem)
{
  int *cycle      = distmem->labeled;     /* cycle index */
  int *cc         = distmem->cc;          /* connected component w/vertex */
  component_t *components = distmem->components; /* list of components */
  int cyc_ind, cc_ind;
  int i;

  /* paths/cycles of length 2 are in their own component, and
   * would be filtered out by the cc_ind >=0 test below,
   * so catch them here
   */
  if (cur_cyc == cycle[cur_cyc] && cur_cyc == compno)
    return cur_cyc;

  for (i=cur_cyc; i<gsize; i += 2) {
    cyc_ind = cycle[i];
    cc_ind = cc[i];

    if (i == cyc_ind         /* leftmost vertex of cycle */
	&& cc_ind >= 0
	&& components[cc_ind].index == compno   /* in correct component */
	)
      return i;
  }

  /* no more cycles found */
  return -1;
}

/* Traverse path or cycle whose leftmost vertex is cyc_ind.
 * lr_or_cc = TEXGR_CC:   produce TeX commands to draw it.
 *            TEXGR_RAW:  list path/cycle in raw text format
 */
void texgraph_drawcyc(int *pilabels, int *galabels,
		      int num_genes, int num_chromosomes,
		      int showtails, int showlabels,
		      distmem_t *distmem,
		      int cyc_index,
		      int lr_or_cc)
{

  int gsize = 2*num_genes + 2;
  int lowcap = num_genes - 2*num_chromosomes + 1;

  int *greyEdges  = distmem->greyEdges;

  /*  int *cycle      = distmem->labeled; */

  int *Bcap       = distmem->Bcap;
  int *Gcap       = distmem->Gcap;
  int *ptype      = distmem->ptype;

  int start;
  int i, lastv, nextv, j;

  int color;
  int visible;

  /* tail cycle? */
  int pilab;
  int s;
  int istail = FALSE;

  int len;               /* length of path or cycle */
  int n_ori=0;           /* number of oriented edges */
  int n_chrom_used;      /* number of chromomes used in path/cycle */
  int n_ichrom;          /* number of interchromosomal gray edges */


  int pt = ptype[cyc_index];  /* path or cycle type */

  /* In some stages of the graph, tails may still be adjacencies.
   * In the left-to-right drawing we had to test each vertex separately.
   * Here, the whole cycle is that way, so do the test once.
   */


  j = greyEdges[cyc_index];
  istail = (j == V_TAIL);
  if (j == V_ADJ && num_chromosomes > 0) {
    /* adjacency at start or end is a tail */
    if (cyc_index == 0 || cyc_index == gsize-1)
      istail = TRUE;

    else {
      /* test if adjacency joins rcap to lcap */

      s = 1;
      pilab = pilabels[(cyc_index-1)/2];

      /* normally cyc_index is even, so this step should always occur */
      if ((cyc_index & 1) == 0)
	s = -s;

      if (pilab < 0) {
	pilab = -pilab;
	s = -s;
      }

      if (pilab >= lowcap) {
	pilab -= lowcap;
	if ((pilab & 1) == 1)
	  s = -s;

	if (s == 1)
	  istail = TRUE;
      }
    }
  }



  /* if it's a path, start from one end */
  /* pi-gamma or gamma-pi path: start from pi; else arbitrary */
  start = cyc_index;
  if (pt == P_GP) {
    start = Gcap[cyc_index];
  } else if (pt == P_PG) {
    start = Bcap[cyc_index];
  } else if (Bcap[cyc_index] >= 0)
    start = Bcap[cyc_index];



  if (lr_or_cc == TEXGR_CC) {
    if (pt == P_NONE || pt == P_CYCLE)
      fprintf(outfile,"%% cycle %d\n", cyc_index);
    else
      fprintf(outfile,"%% path %d\n", cyc_index);
  } else {
    len = texgraph_cyclen(distmem, start, num_chromosomes,
			  &n_ori, &n_chrom_used, &n_ichrom);
    switch(pt) {
    case P_PP:
      fprintf(outfile, "%d pp", len); break;   /* pi-pi path */
    case P_PG: case P_GP:
      fprintf(outfile, "%d pg", len); break;   /* pi-gamma path */
    case P_GG:
      fprintf(outfile, "%d gg", len); break;   /* gamma-gamma path */
    case P_NONE: case P_CYCLE:
    default:
      if (istail) {
	fprintf(outfile, "%d t", len); break;    /* simple tail */
      } else {
	fprintf(outfile, "%d c", len); break;    /* cycle */
      }
    }
    fprintf(outfile," ");
  }





  /* lastv = NULL
   * COLOR = gray
   * while (...) {
   *    Draw vertex i
   *    Draw labels attached to vertex i
   *    if (lastv != NULL)
   *       draw COLOR edge from lastv to i
   *    figure out next vertex, update lastv and i
   *    COLOR = !COLOR
   */

  i = start;
  lastv = -1;
  color = 1;

  visible = SHOWV(i);

  while (1) {
    /* draw the vertex, unless we are closing the cycle */
    if (i == start && lastv >= 0) {
    } else {
      /* draw the vertex */
      if (lr_or_cc == TEXGR_CC) {
	if (visible) {
	  fprintf(outfile,"\\vert{%d}{%s}",
		  i, vlabel(pilabels,gsize,i));
	} else {
	  fprintf(outfile,"\\vert[\\pvert]{%d}{%s}",
		  i, vlabel(pilabels,gsize,i));
	}


	/* doubled label */
	if (showlabels & TEXGR_U) {
	  if (showlabels & TEXGR_AB) {
	    fprintf(outfile,
		    "\\showpiu{%s}{%s}",
		    vlabel(pilabels, gsize, i),
		    vlabel_ab(pilabels, gsize, i));
	  } else {
	    fprintf(outfile,
		    "\\showpiu{%s}{%d}",
		    vlabel(pilabels, gsize, i),
		    vlabel_u(pilabels, gsize, i));
	  }
	}
      } else { /* lr_or_cc == TEXGR_RAW */
	  /* doubled label */
	  if (showlabels & TEXGR_AB) {
	    fprintf(outfile, "%s ", vlabel(pilabels, gsize, i));
	  } else {
	    fprintf(outfile, "%d ", vlabel_u(pilabels, gsize, i));
	  }
      }


      /* mark if tail, pi-cap, or gamma-tail */
      if (lr_or_cc == TEXGR_CC) {
	j = greyEdges[i];
	if (j<0) {
	  switch (j) {
	  case V_PICAP:
	    fprintf(outfile,"\\markpi{%s}%%\n", vlabel(pilabels,gsize,i));
	    break;
	  case V_GTAIL:
	    fprintf(outfile,"\\markga{%s}%%\n", vlabel(pilabels,gsize,i));
	    break;

	  default:
	    if (!istail)
	      break;
	  case V_TAIL:
	    if (!showtails)
	      fprintf(outfile,"\\marktail{%s}%%\n", vlabel(pilabels,gsize,i));
	    break;
	  }
	}
      }
    }

    /* draw edge preceding vertex */
    if ((lr_or_cc == TEXGR_CC) && visible && lastv >= 0) {
      fprintf(outfile,
	      color ? "\\gedge" : "\\bedge");
      fprintf(outfile,"{%s}", vlabel(pilabels, gsize, lastv));
      fprintf(outfile,"{%s}", vlabel(pilabels, gsize, i));
    }

    /* follow edge to next vertex */
    if (color == 0) {
      /* last edge was black, this one is gray */

      color = 1;    /* next edge is gray */

      nextv = greyEdges[i];
      if (nextv < 0) {
	if (nextv == V_ADJ || nextv == V_TAIL) {
	  if (i % 2 == 0)
	    nextv = i+1;
	  else
	    nextv = i-1;
	} else
	  /* end of a path, stop */
	  break;
      }
    } else {
      /* last edge was gray, this one is black */

      /* just closed cycle, stop */
      if (i == start && lastv >= 0)
	break;

      /* next edge is black */
      color = 0;

      if (i % 2 == 0)
	nextv = i+1;
      else
	nextv = i-1;
    }

    lastv = i;
    i = nextv;
  }

  if (lr_or_cc == TEXGR_CC) {
    fprintf(outfile,"%%\n");
  } else { /* TEXGR_RAW */
    /* print statistics about path/cycle */
    fprintf(outfile,"\t%% ");

    if (n_ori == 0) {
      fprintf(outfile, "u");
    } else {
      fprintf(outfile, "o%d", n_ori);
    }

    if (num_chromosomes > 0 && n_chrom_used > 1) {
      fprintf(outfile, " i%d ch%d", n_ichrom, n_chrom_used);
    }

    
    fprintf(outfile,"\n");
  }
}

int texgraph_cyclen(
		    distmem_t *distmem,
		    int start,           /* vertex to start traversal from */
		    int num_chromosomes,

		    /* return values: */
		    int *n_ori,             /* # oriented edges */
		    int *n_chrom,           /* # chromosomes in span */
		    int *n_ichrom           /* # interchrom grey edges */
		    )
{
  int len = 0;
  int i = start;
  int nextv;
  int *greyEdges  = distmem->greyEdges;
  int *chromNum   = distmem->chromNum;
  int *chromUsed  = distmem->capMark;
  int c;

  /* mark the chromosomes */
  /* note chromNum[i] = chromosome number 1,2,...,num_chromosomes (1-offset)
   * but chromUsed[c] uses c=0,1,...,num_chromosomes-1 (0-offset)
   */
  if (num_chromosomes > 0) {
    for (c=0; c<num_chromosomes; c++) { chromUsed[c] = 0; }
  }

  *n_ori = 0;            /* count of number of oriented edges */
  *n_ichrom = 0;         /* count of number of interchrom grey edges */


  while (1) {
    /* mark chromosome */
    if (num_chromosomes > 0) {
      chromUsed[chromNum[i]-1] = 1;
    }

    /* next edge: black edge */

    /* end condition for cycle */
    if (i == start && len > 0)
      break;

    len +=2;  /* add both verts of the black edge */

    /* follow black edge */
    if (i % 2 == 0)
      i++;
    else
      i--;


    /* mark chromosome */
    if (num_chromosomes > 0) {
      chromUsed[chromNum[i]-1] = 1;
    }



    /* follow gray edge */
    nextv = greyEdges[i];
    if (nextv < 0) {
      if (nextv == V_ADJ || nextv == V_TAIL) {
	if (i % 2 == 0)
	  nextv = i+1;
	else
	  nextv = i-1;
      } else
	/* end condition for path */
	break;
    }

    if ((i-nextv) % 2 == 0) { ++*n_ori; }
    if (chromNum[i] != chromNum[nextv]) { ++*n_ichrom; }

    i = nextv;
  }


  /* count # chromosomes used */
  *n_chrom = 0;
  if (num_chromosomes > 0) {
    for (c=0; c<num_chromosomes; c++) {
      if (chromUsed[c]) { ++*n_chrom; }
    }
  }
  
  return len;
}



void texgraph_brackets(int *pilabels, int num_genes, int num_chromosomes)
{
  int lowcap = num_genes - 2*num_chromosomes + 1;

  int lc_ptr = 0;
  int rc_ptr = -1;

  int chr;
  int gsize = 2*num_genes + 2;


  fprintf(outfile,"%% chromosome brackets\n");
  for (chr = 1; chr <= num_chromosomes; chr++) {
    /* left cap is one past previous right cap */
    lc_ptr = rc_ptr + 1;

    /* find right cap */
    for (rc_ptr = lc_ptr+1;
	 pilabels[rc_ptr]<lowcap && pilabels[rc_ptr] > -lowcap;
	 rc_ptr++) {};

    fprintf(outfile, "\\underbracket{%s}",
	    vlabel(pilabels, gsize, 2*lc_ptr+1));
    fprintf(outfile, "{%s}", vlabel(pilabels, gsize, 2*rc_ptr+2));
    fprintf(outfile,"{chromosome %d}\n", chr);
  }
}


/* draw signed entries of pi */
void texgraph_pis(int *pilabels, int num_genes, int num_chromosomes)
{
  int lowcap = num_genes - 2*num_chromosomes + 1;
  int gsize = 2*num_genes + 2;
  int countout = 0;
  int i;
  int pilab;
  
  fprintf(outfile,"%% pi, signed entries\n");

  countout = 0;
  for (i=1; i<gsize-2; i += 2) {
    pilab = pilabels[(i-1)/2];

    if (countout % 3 != 0)
      fprintf(outfile," ");

    fprintf(outfile,"\\showpis{%s}", vlabel(pilabels, gsize, i));
    fprintf(outfile,"{%s}", vlabel(pilabels,gsize,i+1));
    if ((pilab <= -lowcap) || (pilab >= lowcap))
      fprintf(outfile,"{\\bpcap{%d}}", pilab);
    else
      fprintf(outfile,"{%d}", pilab);
    
    countout++;

    if (countout % 3 == 0) {
      fprintf(outfile,"%%\n");
    }
  }

  fprintf(outfile,"%%\n%%\n");
}
