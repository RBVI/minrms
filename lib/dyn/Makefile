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
#ANSI_CPP = cxx -ptr $(FAST_LOCAL_DIRECTORY) -std strict_ansi
#L_COMP = ar -q
#CFLAGS = -c
#LFLAGS = 
INSTALL_PATH = .

# clear the list of suffixes
.SUFFIXES:

.SUFFIXES: .c .cc .o .a
.cc.o:; $(ANSI_CPP) $(CFLAGS) $(INCLUDES) $(DEFINES) $< 

LIBS = \
-L../pair_alignment -lpair_alignment\
-L../superimpose/libsuperimpose.a \
-L../biopolymer -lbiopolymer \
-L../eigen -leigen \
-L../global_utils -lglobal_utils \
-L../vect_3d -lvect_3d \
-L../pdb++ -lpdb++

LIB_FILES = \
../pair_alignment/libpair_alignment.a \
../superimpose/libsuperimpose.a \
../biopolymer/libbiopolymer.a \
../eigen/libeigen.a \
../global_utils/libglobal_utils.a \
../vect_3d/libvect_3d.a \
../pdb/libpdb++.a


INCLUDES = \
-I../pair_alignment \
-I../superimpose \
-I../biopolymer \
-I../eigen \
-I../global_utils \
-I../vect_3d \
-I../pdb++

#Default target: create a library
libdyn.a: \
distance_metric.o \
dyn3d.o \
dynNW.o
	$(L_COMP)        \
    libdyn.a   \
	distance_metric.o \
	dyn3d.o \
	dynNW.o

install:
	$(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" INSTALL_PATH="$(INSTALL_PATH)" libdyn.a
	-cp -f libdyn.a $(INSTALL_PATH)


#create the demonstration program:
demoNW: demoNW.o $(LIB_FILES)
	$(ANSI_CPP) $(LFLAGS) -o demoNW \
    demoNW.o \
    $(LIBS)


GENERATED_FILES = *.o core lib*.a demo*.o ii_files *.tmp
clean:
	rm -rf $(GENERATED_FILES)

distclean:	clean
	rm -f libdyn.a

OBJECT_SRC = \
demoNW.cc \
distance_metric.cc \
dyn3d.cc \
dynNW.cc


RCS_FILES = $(OBJECT_SRC)\
distance_metric.h \
dyn.h \
dyn3d.h \
dynNW.h


checkout:
	co -q $(RCS_FILES)


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
