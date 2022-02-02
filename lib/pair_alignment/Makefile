SHELL = /bin/sh

#FAST_LOCAL_DIRECTORY = /usr/tmp/minrms_compile_temp

#default parameters (can be overriden by caller)
# for SGI:
#ANSI_C = cc
#ANSI_CPP = CC -g -DDEBUG -LANG:std -DSUPPORTS_STD_NAMESPACE
#L_COMP = CC -ar -o
#CFLAGS = -c
#LFLAGS = 
# for Alpha:
#ANSI_C = cc
#ANSI_CPP = cxx -ptr $(FAST_LOCAL_DIRECTORY) -std strict_ansi -DCOMPILER_BUG1 -DDEBUG
#L_COMP = ar -q
#CFLAGS = -c
#LFLAGS = 

INSTALL_PATH = .

# clear the list of suffixes
.SUFFIXES:

.SUFFIXES: .c .cc .o .a
.cc.o:; $(ANSI_CPP) $(CFLAGS) $(INCLUDES) $(DEFINES) $< 


INCLUDES = \
-I../dyn \
-I../biopolymer \
-I../superimpose \
-I../fast_rot_metric \
-I../eigen \
-I../vect_3d \
-I../parse_utils \
-I../global_utils \
-I../pdb


#Default target: create a library
libpair_alignment.a: pair_alignment.o pair_align_settings.o
	$(L_COMP)        \
    libpair_alignment.a \
    pair_alignment.o \
    pair_align_settings.o

install:
	$(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" INSTALL_PATH="$(INSTALL_PATH)" libpair_alignment.a
	-cp -f libpair_alignment.a $(INSTALL_PATH)

rCS_ILES = pair_alignment.cc pair_alignment.h \
 pair_align_settings.cc pair_align_settings.h

checkout:
	co -q $(RCS_FILES)

GENERATED_FILES = *.o *.a core lib*.a ii_files *.tmp

clean:
	rm -rf $(GENERATED_FILES)

distclean:	clean
	rm -f libpair_alignment.a

OBJECT_SRC = pair_alignment.cc pair_align_settings.cc

depend:
	mv Makefile Makefile.tmp
# The next line erases everything after the special "DO NOT MOVE...###"
# line below and copies it to Makefile.tmp
	sed -n '1,/DE\PEND/p' < Makefile.tmp > Makefile
# Generate the dependencies (using the compiler's -M option) and append
# them to the makefile
	$(ANSI_CPP) -M $(CFLAGS) $(INCLUDES) $(DEFINES) $(OBJECT_SRC) >> Makefile
	rm -f Makefile.tmp


# "make depend" requires that you...
# DO NOT MOVE OR DELETE (or place your own text after) THE FOLLOWING LINE:
### DEPEND
