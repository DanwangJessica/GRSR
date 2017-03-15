# Makefile
#    Makefile for GRIMM.
#
# Copyright (C) 2001-2006 The Regents of the University of California
# by Glenn Tesler
#
# Contains code from Makefile in GRAPPA 1.02
# Copyright (C) 2000-2001  The University of New Mexico and
#                          The University of Texas at Austin
# by David A. Bader, Bernard M.E. Moret, Tandy Warnow, Stacia K Wyman, Mi Yan
#
# See file COPYRIGHT for details.
#############################################################################
# This file is part of GRIMM.
#
# GRIMM is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License, Version 2,
# dated June 1991, as published by the Free Software Foundation.
#
# GRIMM is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# Last modified on Tue Aug 1, 2006, by Glenn Tesler


VERS_MAJ = 2
VERS_MIN = 01

EXEC	= grimm

OS    = $(shell uname -s | sed -e "s/\-/_/g" -e "s/\.//g")
ZIP   = gzip
MKDIR = mkdir
TAR   = tar
CP    = cp
RM    = rm
LS    = ls

# GNU C
CC = gcc

# Microsoft C
# CC = cl

# Intel C
# CC = icl

# Portland Group C
# CC = pgcc

# MPI C
# CC = mpicc -DMPBPA

# Default CFLAGS
CFLAGS := -D${OS}
CFLAGS := $(CFLAGS) -DVERS=$(VERS_MAJ).$(VERS_MIN)

# DEBUG options for CFLAGS
#CFLAGS := -DDEBUG $(CFLAGS)
#CFLAGS := -DDEBUGu $(CFLAGS)      # debugging unsigned.c
#CFLAGS := -DDEBUGf $(CFLAGS)      # debugging proper flipping/cocapping
#CFLAGS := -DDEBUGc $(CFLAGS)	   # debugging optimal concatenates
#CFLAGS := -DDEBUGb $(CFLAGS)      # debugging bonding
#CFLAGS := -g $(CFLAGS)


# TESTING option for CFLAGS
# CFLAGS := -DTESTING $(CFLAGS)

ifeq ($(strip $(CC)), gcc)
CFLAGS := -O3 -Wall -DGCC $(CFLAGS)
endif
ifeq ($(strip $(CC)), pgcc)
CFLAGS := -fast $(CFLAGS)
endif
ifeq ($(strip $(CC)), cl)
CFLAGS := /Ot /Ox $(CFLAGS)
endif
ifeq ($(strip $(CC)), icl)
CFLAGS := -O2 $(CFLAGS)
endif

ifeq ($(strip $(CC)), gcc)
ifeq ($(OS), Linux)
#CFLAGS := -march=pentiumpro $(CFLAGS)
CFLAGS := $(CFLAGS)
endif
ifeq ($(OS), SunOS)
CFLAGS := -mv8 $(CFLAGS)
endif
endif

OBJ_SUFFIX := o
EXE_SUFFIX := 

ifeq ($(OS), CYGWIN_NT_40)
OBJ_SUFFIX := obj
EXE_SUFFIX := exe
endif


SOURCES = mcmain uniinvdist mcrdist graph_edit graph_components mcread_input scenario testrev write_data unsigned e_malloc circ_align texgraph opt_scenario ext_function countperms unsignedhc matrixmisc

HEADERS = structs.h $(filter-out mcmain.h,$(addsuffix .h,$(SOURCES)))

OBJS = $(addsuffix .${OBJ_SUFFIX},$(SOURCES))

TSTAMP   := $(shell date +"%y%m%d-%H%M")
DISTBASE  = MPBPA
DISTNAME  = $(DISTBASE)-$(TSTAMP)

RELBASE   = GRIMM
RELNAME   = $(RELBASE)-$(VERS_MAJ).$(VERS_MIN)

LIBS 		:= -lm
ifeq ($(OS), SunOS)
LIBS		:= -lsocket -lnsl ${LIBS}
endif

.PHONY: dist release release_code clean

all: ${EXEC} ${SYM}

${EXEC} :  ${OBJS}
	${CC}  ${CFLAGS}  -o $@  ${OBJS} ${LIBS}

%.${OBJ_SUFFIX}:%.c
	$(CC) $(CFLAGS) -c $< -o $@

${SYM}: ${EXEC}
	ln -s $< $@

dist:
	$(RM) -fr $(DISTNAME) $(DISTNAME).tar $(DISTNAME).tar.gz
	$(MKDIR) $(DISTNAME)
	-$(CP) Makefile COPYRIGHT README CHANGES $(addsuffix .c,$(SOURCES)) $(HEADERS) \
	   $(DISTNAME)/.
	$(TAR) cvf $(DISTNAME).tar $(DISTNAME)
	$(ZIP) $(DISTNAME).tar
	$(RM) -fr $(DISTNAME)
	$(LS) -l $(DISTNAME).tar.gz

release:
	$(RM) -fr $(RELNAME) $(RELNAME).tar $(RELNAME).tar.gz
	$(MKDIR) $(RELNAME)
	-$(CP) Makefile COPYRIGHT README CHANGES $(RELNAME)/.
	$(MAKE) $(addprefix $(RELNAME)/,$(addsuffix .c,$(SOURCES)) $(HEADERS))
	$(TAR) cvf $(RELNAME).tar $(RELNAME)
	$(ZIP) $(RELNAME).tar
	$(RM) -fr $(RELNAME)
	$(LS) -l $(RELNAME).tar.gz

$(RELNAME)/%.c : %.c
	-$(CP) $< $(RELNAME)/.

$(RELNAME)/%.h : %.h
	-$(CP) $< $(RELNAME)/.

clean:
	-rm -f *.${OBJ_SUFFIX} *~ ${EXEC} ${EXEC}.${EXE_SUFFIX} ${SYM} \
	       core gmon.out

# DO NOT DELETE

$(OBJS): Makefile
circ_align.o: circ_align.c mcstructs.h scenario.h circ_align.h
countperms.o: countperms.c uniinvdist.h mcstructs.h scenario.h mcrdist.h unsigned.h graph_edit.h graph_components.h e_malloc.h countperms.h
e_malloc.o: e_malloc.c e_malloc.h mcstructs.h
ext_function.o: ext_function.c mcstructs.h uniinvdist.h mcrdist.h unsigned.h e_malloc.h texgraph.h ext_function.h
graph_components.o: graph_components.c mcstructs.h graph_components.h
graph_edit.o: graph_edit.c mcstructs.h graph_edit.h
matrixmisc.o: matrixmisc.c matrixmisc.h e_malloc.h mcread_input.h
mcmain.o: mcmain.c mcstructs.h uniinvdist.h mcrdist.h mcread_input.h write_data.h scenario.h testrev.h opt_scenario.h unsigned.h e_malloc.h circ_align.h texgraph.h ext_function.h countperms.h unsignedhc.h matrixmisc.h
mcrdist.o: mcrdist.c uniinvdist.h mcstructs.h mcrdist.h graph_edit.h graph_components.h e_malloc.h
mcread_input.o: mcread_input.c mcread_input.h mcstructs.h e_malloc.h
opt_scenario.o: opt_scenario.c uniinvdist.h mcstructs.h mcrdist.h scenario.h testrev.h opt_scenario.h write_data.h mcread_input.h graph_edit.h
scenario.o: scenario.c uniinvdist.h mcstructs.h scenario.h mcrdist.h e_malloc.h
testrev.o: testrev.c uniinvdist.h mcstructs.h scenario.h mcrdist.h mcread_input.h e_malloc.h write_data.h testrev.h
texgraph.o: texgraph.c uniinvdist.h mcstructs.h mcrdist.h graph_edit.h graph_components.h texgraph.h
uniinvdist.o: uniinvdist.c uniinvdist.h mcstructs.h graph_edit.h graph_components.h mcrdist.h
unsigned.o: unsigned.c uniinvdist.h mcstructs.h scenario.h mcrdist.h unsigned.h graph_edit.h graph_components.h e_malloc.h circ_align.h
unsignedhc.o: unsignedhc.c mcstructs.h uniinvdist.h unsignedhc.h scenario.h mcrdist.h e_malloc.h write_data.h matrixmisc.h
write_data.o: write_data.c mcstructs.h
circ_align.h: mcstructs.h
mcread_input.h: mcstructs.h
uniinvdist.h: mcstructs.h
