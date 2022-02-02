SHELL = /bin/sh

FAST_LOCAL_DIRECTORY = /usr/tmp/minrms_compile_temp

#default parameters (can be overriden by caller)
# for SGI:
#ANSI_C = cc
#ANSI_CPP = CC
#L_COMP = CC -ar -o
#CFLAGS = -c -LANG:ansi-for-init-scope=on
#LFLAGS = 
# for Alpha:
#ANSI_C = cc
#ANSI_CPP = cxx -ptr $(FAST_LOCAL_DIRECTORY) -std strict_ansi -DCOMPILER_BUG1
#L_COMP = ar -q
#CFLAGS = -c
#LFLAGS = 
INSTALL_PATH = .

#GLOBAL_UTILS = ../utils/utils.h ../utils/debug.h ../utils/dbg_mem_leak.h

# clear the list of suffixes
.SUFFIXES:

.SUFFIXES: .c .cc .o .a
.cc.o:; $(ANSI_CPP) $(CFLAGS) $(INCLUDES) $(DEFINES) $< 


INCLUDES = \
-I../dyn \
-I../pair_alignment \
-I../biopolymer \
-I../superimpose \
-I../fast_rot_metric \
-I../vect_3d \
-I../global_utils \
-I../pdb++

#Default target: create a library
libor_gen.a: or_gen.o orientation_filter.o
	$(L_COMP) \
	libor_gen.a \
	or_gen.o \
    orientation_filter.o

install:
	$(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" INSTALL_PATH="$(INSTALL_PATH)" libor_gen.a
	-cp -f libor_gen.a $(INSTALL_PATH)

checkout:
	co -q $(RCS_FILES)

GENERATED_FILES = *.o *.a core ii_files *.tmp

clean:
	rm -rf $(GENERATED_FILES)

distclean:	clean
	rm -f libor_gen.a


OBJECT_SRC = or_gen.cc orientation_filter.cc

RCS_FILES = $(OBJECT_SRC) \
or_gen.h \
orientation_filter.h


depend:
	mv Makefile Makefile.tmp
# The next line erases everything after the special "DO NOT MOVE...###"
# line below and copies it to Makefile.tmp
	sed -n '1,/DE\PEND/p' < Makefile.tmp > Makefile
# Generate the dependencies (using the compiler's -M option) and append
# them to the makefile
	$(ANSI_CPP) -M $(CFLAGS) $(INCLUDES) $(OBJECT_SRC) >> Makefile
	rm -f Makefile.tmp


# "make depend" requires that you...
# DO NOT MOVE OR DELETE (or place your own text after) THE FOLLOWING LINE:
### DEPEND
