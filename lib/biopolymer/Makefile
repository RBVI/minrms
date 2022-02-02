SHELL = /bin/sh

#FAST_LOCAL_DIRECTORY = /tmp/minrms_compile_temp

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
#L_COMP = ar -r
#CFLAGS = -c
#LFLAGS = 
INSTALL_PATH = .

# clear the list of suffixes
.SUFFIXES:

.SUFFIXES: .c .cc .o .a
.cc.o:; $(ANSI_CPP) $(CFLAGS) $(INCLUDES) $(DEFINES) $< 

LIBS = \
-L../pdb++ -lpdb++ \
-L../vect_3d -lvect_3d \
-L../parse_utils -lparse_utils \
-L../global_utils -lglobal_utils

LIB_FILES = \
../pdb++/libpdb++.a \
../vect_3d/libvect_3d.a \
../parse_utils/libparse_utils.a \
../global_utils/libglobal_utils.a


INCLUDES = \
-I../pdb++ \
-I../vect_3d \
-I../parse_utils \
-I../global_utils


#Default target: create a library
libbiopolymer.a: \
apply_trans.o \
atom.o \
intervals.o \
biopolymer.o \
load_pdb.o \
mol2sequence.o \
parse_set.o \
pdb_id.o \
residue_names.o
	$(L_COMP) \
        libbiopolymer.a \
	apply_trans.o \
	atom.o \
	intervals.o \
	biopolymer.o \
	load_pdb.o \
	mol2sequence.o \
	parse_set.o \
	pdb_id.o \
	residue_names.o


install:
	$(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" INSTALL_PATH="$(INSTALL_PATH)" libbiopolymer.a
	-cp -f libbiopolymer.a $(INSTALL_PATH)

RCS_FILES = \
apply_trans.cc \
apply_trans.h \
atom.cc \
atom.h \
intervals.cc \
intervals.h \
biopolymer.cc \
biopolymer.h \
load_pdb.cc \
load_pdb.h \
mol2sequence.cc \
mol2sequence.h \
parse_set.cc \
parse_set.h \
pdb_id.cc \
pdb_id.h \
residue_names.cc \
residue_names.h

checkout:
	co -q $(RCS_FILES)

GENERATED_FILES = *.o *.a core lib*.a ii_files/* *.tmp

clean:
	rm -f $(GENERATED_FILES)

distclean:
	$(MAKE) clean
	rm -f $(INSTALL_PATH)/libbiopolymer.a

OBJECT_SRC = \
apply_trans.cc \
atom.cc \
intervals.cc \
biopolymer.cc \
load_pdb.cc \
mol2sequence.cc \
parse_set.cc \
pdb_id.cc \
residue_names.cc

depend:
	mv Makefile Makefile.tmp
# The next line erases everything after the special "DO NOT MOVE...###"
# line below and copies it to Makefile.tmp
	sed -n '1,/DE\PEND/p' < Makefile.tmp > Makefile
# Generate the dependencies (using the compiler's -M option) and append
# them to the makefile
	$(ANSI_CPP) -M $(CFLAGS) $(INCLUDES) $(OBJECT_SRC) >> Makefile
	rm -f Makefile.tmp
# Here's the version from gregc's makefile
#	${ANSI_CPP} -M ${CFLAGS} $(INCLUDES) ${OBJECT_SRC} | sed s\|/usr/local/otf\|\$$\{OTFROOT\}\| >> Makefile
# more generally, use:
#	${ANSI_CPP} -M ${CFLAGS} $(INCLUDES) ${OBJECT_SRC} | sed s\|REPLACE_THIS\|\$WITH_THIS\| >> Makefile


# "make depend" requires that you...
# DO NOT MOVE OR DELETE (or place your own text after) THE FOLLOWING LINE:
### DEPEND
