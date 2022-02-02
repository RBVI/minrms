SHELL = /bin/sh

FAST_LOCAL_DIRECTORY = /usr/tmp/minrms_compile_temp

#default parameters (can be overriden by caller)
# for SGI:
#ANSI_C = cc
#ANSI_CPP = CC -g -LANG:std
#L_COMP = CC -ar -o
#CFLAGS = -c -LANG:ansi-for-init-scope=on
#LFLAGS = 
# for Alpha:
#ANSI_C = cc
#ANSI_CPP = cxx -ptr $(FAST_LOCAL_DIRECTORY) -std strict_ansi -DSUPPORTS_STD_NAMESPACE -DCOMPILER_BUG1
#L_COMP = ar -r
#CFLAGS = -c
#LFLAGS = 

INTERNAL_LIB_PATH = ../../lib

INSTALL_PATH = .

# clear the list of suffixes
.SUFFIXES:

.SUFFIXES: .c .cc .o .a
.cc.o:; $(ANSI_CPP) $(CFLAGS) $(INCLUDES) $(DEFINES) $< 

LIBS = \
-L$(INTERNAL_LIB_PATH)/dyn -ldyn \
-L$(INTERNAL_LIB_PATH)/or_gen -lor_gen \
-L$(INTERNAL_LIB_PATH)/pair_alignment -lpair_alignment\
-L$(INTERNAL_LIB_PATH)/biopolymer -lbiopolymer \
-L$(INTERNAL_LIB_PATH)/superimpose -lsuperimpose \
-L$(INTERNAL_LIB_PATH)/fast_rot_metric -lfast_rot_metric \
-L$(INTERNAL_LIB_PATH)/eigen -leigen \
-L$(INTERNAL_LIB_PATH)/vect_3d -lvect_3d \
-L$(INTERNAL_LIB_PATH)/parse_utils -lparse_utils \
-L$(INTERNAL_LIB_PATH)/global_utils -lglobal_utils \
-L$(INTERNAL_LIB_PATH)/pdb++ -lpdb++

LIB_FILES = \
$(INTERNAL_LIB_PATH)/dyn/libdyn.a \
$(INTERNAL_LIB_PATH)/or_gen/libor_gen.a \
$(INTERNAL_LIB_PATH)/pair_alignment/libpair_alignment.a \
$(INTERNAL_LIB_PATH)/biopolymer/libbiopolymer.a \
$(INTERNAL_LIB_PATH)/superimpose/libsuperimpose.a \
$(INTERNAL_LIB_PATH)/fast_rot_metric/libfast_rot_metric.a \
$(INTERNAL_LIB_PATH)/eigen/libeigen.a \
$(INTERNAL_LIB_PATH)/vect_3d/libvect_3d.a \
$(INTERNAL_LIB_PATH)/parse_utils/libparse_utils.a \
$(INTERNAL_LIB_PATH)/global_utils/libglobal_utils.a \
$(INTERNAL_LIB_PATH)/pdb++/libpdb++.a

INCLUDES = \
-I$(INTERNAL_LIB_PATH)/dyn \
-I$(INTERNAL_LIB_PATH)/or_gen \
-I$(INTERNAL_LIB_PATH)/pair_alignment \
-I$(INTERNAL_LIB_PATH)/biopolymer \
-I$(INTERNAL_LIB_PATH)/superimpose \
-I$(INTERNAL_LIB_PATH)/fast_rot_metric \
-I$(INTERNAL_LIB_PATH)/eigen \
-I$(INTERNAL_LIB_PATH)/vect_3d \
-I$(INTERNAL_LIB_PATH)/parse_utils \
-I$(INTERNAL_LIB_PATH)/global_utils \
-I$(INTERNAL_LIB_PATH)/pdb++

OBJECT_FILES = \
main.o \
minrms_parser.o \
search_or.o \
search_or_dyn3d.o \
search_or_nw.o

OBJECT_SRC = \
main.cc \
minrms_parser.cc \
search_or.cc \
search_or_dyn3d.cc \
search_or_nw.cc

rCS_FILES = \
$(OBJECT_FILES) \
minrms_settings.h \
minrms_parser.h \
search_or.h \
search_or_dyn3d.h \
search_or_nw.h

#Default target: the minrms binary
minrms: $(OBJECT_FILES) $(LIB_FILES)
	$(ANSI_CPP) $(LFLAGS) -o minrms \
	$(OBJECT_FILES) \
    $(LIBS)

install:
	$(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" INSTALL_PATH="$(INSTALL_PATH)" minrms
	cp -f minrms $(INSTALL_PATH)

checkout:
	co -q $(RCS_FILES)

GENERATED_FILES = $(OBJECT_FILES) *.a core minrms minrms.pure ii_files *.tmp

clean:
	rm -rf $(GENERATED_FILES)

distclean:	clean
	rm -f minrms

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
