# Useful directories

#this is the path where the samtools source package was built
#so all the samtools headers (*.h) and (*.a) library files are in there
SAM := ./samtools-0.1.18

GDIR := ../gclib

#SEARCHDIRS := -I${SAM} -I${GDIR}
SEARCHDIRS := -I${GDIR}
SYSTYPE :=     $(shell uname)

MACHTYPE :=     $(shell uname -m)

CC      := g++
BASEFLAGS  = -Wall ${SEARCHDIRS} $(MARCH) -D_FILE_OFFSET_BITS=64 \
-D_LARGEFILE_SOURCE -fno-exceptions -fno-rtti -fno-strict-aliasing \
-D_REENTRANT

ifneq (,$(findstring release,$(MAKECMDGOALS)))
  CFLAGS = -O2 -msse2 -DNDEBUG $(BASEFLAGS)
  #LDFLAGS = -L${SAM}
else
  CFLAGS = -g -DGDEBUG -DDEBUG $(BASEFLAGS)
  LDFLAGS = -g ${LDFLAGG}
endif

%.o : %.cpp
	${CC} ${CFLAGS} -c $< -o $@

# C/C++ linker

LINKER := g++
LIBS := 

OBJS := ${GDIR}/GBase.o ${GDIR}/GStr.o ${GDIR}/GArgs.o ${GDIR}/gdna.o ./GapAssem.o

#ifndef NOTHREADS
# OBJS += ${GDIR}/GThreads.o 
#endif

#ifdef GDEBUG
# OBJS += ${GDIR}/proc_mem.o
#endif

.PHONY : all debug release
all: paf2msa
#all:    pwasm pwaor
debug : all
release : all

paf2msa :  ./paf2msa.o ./GapAssem.o ${OBJS}
	${LINKER} -o $@ ${filter-out %.a %.so, $^} $(LDFLAGS) ${LIBS}

bamcons :  ./bamcons.o ${GDIR}/GFastaIndex.o ${GDIR}/GFaSeqGet.o ${GDIR}/GBam.o ${OBJS} 
	${LINKER} $(LDFLAGS) -o $@ ${filter-out %.a %.so, $^} -L${SAM} ${LIBS} -lbam

pwasm :  ./pwasm.o ${GDIR}/GCdbYank.o ${GDIR}/gcdb.o ${OBJS}
	${LINKER} $(LDFLAGS) -o $@ ${filter-out %.a %.so, $^} ${LIBS}

nrcl:  ./nrcl.o ${GDIR}/GBase.o ${GDIR}/GStr.o ${GDIR}/GArgs.o
	${LINKER} $(LDFLAGS) -o $@ ${filter-out %.a %.so, $^}

tclust:  ./tclust.o ${GDIR}/GBase.o ${GDIR}/GStr.o ${GDIR}/GArgs.o
	${LINKER} $(LDFLAGS) -o $@ ${filter-out %.a %.so, $^}

sclust:  ./sclust.o ${GDIR}/GBase.o ${GDIR}/GStr.o ${GDIR}/GArgs.o
	${LINKER} $(LDFLAGS) -o $@ ${filter-out %.a %.so, $^}

./GapAssem.o: GapAssem.h

pwaor :  ./pwaor.o ./GapAssem.o ${GDIR}/GCdbYank.o ${GDIR}/gcdb.o ${OBJS}
	${LINKER} -o $@ ${filter-out %.a %.so, $^} $(LDFLAGS) ${LIBS}

# target for removing all object files

.PHONY : tidy
tidy::
	${RM} core* pwasm pwasm.exe pwaor pwaor.exe nrcl nrcl.exe *clust *clust.exe *.o ${OBJS}

# target for removing all object files

.PHONY : clean
clean:: tidy
	
