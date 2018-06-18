# Useful directories

GDIR := ../gclib

INCDIRS := -I${GDIR}
SYSTYPE :=     $(shell uname)

MACHTYPE :=     $(shell uname -m)

CC      := g++

BASEFLAGS  = -Wall -Wextra ${INCDIRS} -fno-exceptions -fno-rtti \
 -D_REENTRANT -fno-strict-aliasing

# C/C++ linker and core libs

LINKER := g++
LDFLAGS := 
LIBS := 

ifneq (,$(filter %release %static, $(MAKECMDGOALS)))
  # -- release build
  CFLAGS := -O2 -msse2 -DNDEBUG $(BASEFLAGS)
  LDFLAGS := $(LDFLAGS)
  LIBS := $(LIBS)
  ifneq (,$(findstring static,$(MAKECMDGOALS)))
    LDFLAGS += -static-libstdc++ -static-libgcc
  endif
else # debug build
  ifneq (,$(filter %memcheck %memdebug, $(MAKECMDGOALS)))
     #make memcheck : use the statically linked address sanitizer in gcc 4.9.x
     GCCVER49 := $(shell expr `g++ -dumpversion | cut -f1,2 -d.` \>= 4.9)
     ifeq "$(GCCVER49)" "0"
       $(error gcc version 4.9 or greater is required for this build target)
     endif
     CFLAGS := -fno-omit-frame-pointer -fsanitize=undefined -fsanitize=address
     GCCVER5 := $(shell expr `g++ -dumpversion | cut -f1 -d.` \>= 5)
     ifeq "$(GCCVER5)" "1"
       CFLAGS += -fsanitize=bounds -fsanitize=float-divide-by-zero -fsanitize=vptr
       CFLAGS += -fsanitize=float-cast-overflow -fsanitize=object-size
       #CFLAGS += -fcheck-pointer-bounds -mmpx
     endif
     CFLAGS += $(BASEFLAGS)
     CFLAGS := -g -DDEBUG -D_DEBUG -DGDEBUG -fno-common -fstack-protector $(CFLAGS)
     LDFLAGS := -g $(LDFLAGS)
     #LIBS := -Wl,-Bstatic -lasan -lubsan -Wl,-Bdynamic -ldl $(LIBS)
     LIBS := -lasan -lubsan -ldl $(LIBS)
  else
     # regular debug build
     CFLAGS := -g -DDEBUG -D_DEBUG -DGDEBUG $(BASEFLAGS)
     LDFLAGS := -g $(LDFLAGS)
     LIBS := $(LIBS)
  endif
endif



%.o : %.cpp
	${CC} ${CFLAGS} -c $< -o $@

OBJS := ${GDIR}/GBase.o ${GDIR}/GStr.o ${GDIR}/GArgs.o ${GDIR}/gdna.o \
 ${GDIR}/codons.o ${GDIR}/GFastaIndex.o ${GDIR}/GFaSeqGet.o ${GDIR}/GapAssem.o

#ifndef NOTHREADS
# OBJS += ${GDIR}/GThreads.o 
#endif

#ifdef GDEBUG
# OBJS += ${GDIR}/proc_mem.o
#endif

.PHONY : all debug release
all: pafreport
#all:    pwasm pwaor
debug : all
release : all
memcheck : all
memdebug : all 
static : all

pafreport :  ./pafreport.o ${OBJS}
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}

${GDIR}/GapAssem.o: ${GDIR}/GapAssem.h

# target for removing all object files

.PHONY : tidy
tidy::
	${RM} pafreport pafreport.exe pwasm pwasm.exe pwaor pwaor.exe nrcl nrcl.exe *clust *clust.exe *.o ${OBJS} ${GDIR}/codons.o

# target for removing all object files

.PHONY : clean
clean:: tidy
	
