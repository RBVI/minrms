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
#ANSI_CPP = cxx -ptr $(FAST_LOCAL_DIRECTORY) -std strict_ansi -DSUPPORTS_STD_NAMESPACE -DCOMPILER_BUG1
#L_COMP = ar -q
#CFLAGS = -c
#LFLAGS = 
INSTALL_PATH = .

# clear the list of suffixes
.SUFFIXES:

.SUFFIXES: .c .cc .o .a
.cc.o:; $(ANSI_CPP) $(CFLAGS) $(INCLUDES) $(DEFINES) $< 


INCLUDES = \
-I../biopolymer \
-I../pair_alignment \
-I../superimpose \
-I../eigen \
-I../vect_3d \
-I../global_utils \
-I../pdb


#Default target: create a library
libfast_rot_metric.a: \
fast_rot_metric.o \
principle_axis.o
	$(L_COMP)        \
	libfast_rot_metric.a   \
	fast_rot_metric.o \
	principle_axis.o

install:
	$(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" INSTALL_PATH="$(INSTALL_PATH)" libfast_rot_metric.a
	-cp -f libfast_rot_metric.a $(INSTALL_PATH)

RCS_FILES = \
fast_rot_metric.cc \
fast_rot_metric.h \
principle_axis.cc \
principle_axis.h

checkout:
	co -q $(RCS_FILES)

GENERATED_FILES = *.o *.a core lib*.a demo demo.pure ii_files *.tmp

clean:
	rm -rf $(GENERATED_FILES)

distclean:	clean
	rm -f libfast_rot_metric.a

OBJECT_SRC = fast_rot_metric.cc principle_axis.cc

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
