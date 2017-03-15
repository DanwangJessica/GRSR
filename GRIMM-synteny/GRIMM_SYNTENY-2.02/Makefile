# Makefile
#    Makefile for GRIMM-Synteny.
#
# Copyright (C) 2001-2008 The Regents of the University of California
# by Glenn Tesler
#
# See file COPYRIGHT for details.
#############################################################################
# This file is part of GRIMM-Synteny.
#
# GRIMM-Synteny is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License, Version 2,
# dated June 1991, as published by the Free Software Foundation.
#
# GRIMM-Synteny is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#
# Last modified on Wed Mar 19, 2008, by Glenn Tesler


VERS_MAJ = 2
VERS_MIN = 02



CC = gcc
CCOPT = -O2
#CCOPT = -g -Wall
#CCOPT = -g 
#INCDIRS = -I${GRIMM}
INCDIRS = 
CFLAGS = ${CCOPT} ${INCDIRS} -DVERS=$(VERS_MAJ).$(VERS_MIN)
LIBS = -lm
PROGRAMS = grimm_synt

SOURCES_GRIMM_SYNT = mcmain uniinvdist mcrdist graph_edit graph_components mcread_input scenario testrev write_data unsigned e_malloc circ_align texgraph opt_scenario
SOURCES_GRIMM = uniinvdist mcrdist graph_edit graph_components mcread_input scenario testrev write_data unsigned e_malloc circ_align texgraph opt_scenario

GRIMM_O= uniinvdist.o mcrdist.o graph_edit.o graph_components.o mcread_input.o scenario.o testrev.o write_data.o unsigned.o e_malloc.o circ_align.o texgraph.o opt_scenario.o

#GRIMM_H =

all: ${PROGRAMS}

clean:
	/bin/rm -f *.o

.c.o:
	$(CC) $(CFLAGS) -c $<

OBJS_grimm_synt = grimm_synt.o ckalloc.o hash.o parsestr.o time_stamp.o \
	gscomp.o gs_grimm.o gsfile.o grimm_anch.o RS.o \
	$(GRIMM_O)
grimm_synt: $(OBJS_grimm_synt)
	$(CC) $(CCOPT) -o grimm_synt $(OBJS_grimm_synt) $(LIBS)



# GRIMM-Synteny:
RS.o: RS.c ckalloc.h RS.h
ckalloc.o: ckalloc.c ckalloc.h
grimm_anch.o: grimm_anch.c ckalloc.h time_stamp.h hash.h anctable.h gsy.h grimm_anch.h gsfile.h gscomp.h
grimm_synt.o: grimm_synt.c time_stamp.h hash.h ckalloc.h parsestr.h anctable.h gsy.h gscomp.h gsfile.h grimm_anch.h
gs_grimm.o: gs_grimm.c ckalloc.h mcstructs.h mcrdist.h gs_grimm.h
gscomp.o: gscomp.c ckalloc.h hash.h anctable.h gsy.h gscomp.h mcstructs.h mcrdist.h gs_grimm.h gsfile.h RS.h
gsfile.o: gsfile.c gsfile.h ckalloc.h
hash.o: hash.c hash.h ckalloc.h
parsestr.o: parsestr.c parsestr.h
time_stamp.o: time_stamp.c time_stamp.h

# GRIMM:
circ_align.o: circ_align.c mcstructs.h scenario.h circ_align.h
e_malloc.o: e_malloc.c e_malloc.h mcstructs.h
graph_components.o: graph_components.c mcstructs.h graph_components.h
graph_edit.o: graph_edit.c mcstructs.h graph_edit.h
mcrdist.o: mcrdist.c uniinvdist.h mcstructs.h mcrdist.h graph_edit.h graph_components.h e_malloc.h
mcread_input.o: mcread_input.c mcread_input.h mcstructs.h e_malloc.h
opt_scenario.o: opt_scenario.c uniinvdist.h mcstructs.h mcrdist.h scenario.h testrev.h opt_scenario.h write_data.h mcread_input.h graph_edit.h
scenario.o: scenario.c uniinvdist.h mcstructs.h scenario.h mcrdist.h e_malloc.h
testrev.o: testrev.c uniinvdist.h mcstructs.h scenario.h mcrdist.h mcread_input.h e_malloc.h write_data.h testrev.h
texgraph.o: texgraph.c uniinvdist.h mcstructs.h mcrdist.h graph_edit.h graph_components.h texgraph.h
uniinvdist.o: uniinvdist.c uniinvdist.h mcstructs.h graph_edit.h graph_components.h mcrdist.h
unsigned.o: unsigned.c uniinvdist.h mcstructs.h scenario.h mcrdist.h unsigned.h graph_edit.h graph_components.h e_malloc.h circ_align.h
write_data.o: write_data.c mcstructs.h
mcread_input.h: mcstructs.h
