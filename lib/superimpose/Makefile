SHELL = /bin/sh

FAST_LOCAL_DIRECTORY = /usr/tmp/compile_temp

INSTALL_PATH = .

# clear the list of suffixes
.SUFFIXES:

.SUFFIXES: .c .cc .o .a
.cc.o:; $(ANSI_CPP) $(CFLAGS) $(INCLUDES) $(DEFINES) $< 

LIBS = \
-L../eigen -leigen \
-L../vect_3d -lvect_3d \
-L../global_utils -lglobal_utils

LIB_FILES = \
../eigen/libeigen.a \
../vect_3d/libvect_3d.a \
../global_utils/libglobal_utils.a


INCLUDES = \
-I../eigen \
-I../vect_3d \
-I../global_utils



#Default target: create a library
libsuperimpose.a: \
superimpose.o     
	$(L_COMP)        \
    libsuperimpose.a   \
    superimpose.o

install:
	$(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" INSTALL_PATH="$(INSTALL_PATH)" libsuperimpose.a
	-cp -f libsuperimpose.a $(INSTALL_PATH)

RCS_FILES = superimpose.h superimpose.cc superimpose.h

checkout:
	co -q $(RCS_FILES)

GENERATED_FILES = *.o *.a core lib*.a demo demo.pure ii_files/* *.tmp

clean:
	rm -f $(GENERATED_FILES)

distclean:
	$(MAKE) clean
	rm -f $(INSTALL_PATH)/libsuperimpose.a

OBJECT_SRC = superimpose.cc

depend:
	mv Makefile Makefile.tmp
# The next line erases everything after the special "DO NOT MOVE...###"
# line below and copies it to Makefile.tmp
	sed -n '1,/DE\PEND/p' < Makefile.tmp > Makefile
# Generate the dependencies (using the compiler's -M option) and append
# them to the makefile
	$(ANSI_CPP) -M $(CFLAGS) $(INCLUDES) $(DEFINES) $(OBJECT_SRC) >> Makefile
	rm -f Makefile.tmp
# Here's the version from gregc's makefile
#	${ANSI_CPP} -M ${CFLAGS} $(INCLUDES) $(DEFINES) ${OBJECT_SRC} | sed s\|/usr/local/otf\|\$$\{OTFROOT\}\| >> Makefile
# more generally, use:
#	${ANSI_CPP} -M ${CFLAGS} $(INCLUDES) $(DEFINES) ${OBJECT_SRC} | sed s\|REPLACE_THIS\|\$WITH_THIS\| >> Makefile

# "make depend" requires that you...
# DO NOT MOVE OR DELETE (or place your own text after) THE FOLLOWING LINE:
### DEPEND
