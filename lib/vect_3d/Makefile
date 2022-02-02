SHELL = /bin/sh

FAST_LOCAL_DIRECTORY = /usr/tmp/compile_temp

INSTALL_PATH = .


# clear the list of suffixes
.SUFFIXES:

.SUFFIXES: .c .cc .o .a
.cc.o:; $(ANSI_CPP) $(CFLAGS) $(INCLUDES) $(DEFINES) $< 


LIBS = \
-L../global_utils -lglobal_utils

LIB_FILES = \
-L../global_utils/libglobal_utils.a

INCLUDES = \
-I../global_utils

#Default target: create a library
libvect_3d.a: \
vect_3d.o
	$(L_COMP) \
	libvect_3d.a \
	vect_3d.o

install:
	$(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" INSTALL_PATH="$(INSTALL_PATH)" libvect_3d.a
	-cp -f libvect_3d.a $(INSTALL_PATH)

RCS_FILES = \
vect_3d.cc \
vect_3d.h

checkout:
	co -q $(RCS_FILES)

GENERATED_FILES = *.o *.a core ii_files/* *.tmp

clean:
	rm -f $(GENERATED_FILES)

distclean:
	$(MAKE) clean
	rm -f $(INSTALL_PATH)/libvect_3d.a


OBJECT_SRC = vect_3d.cc

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
#	${ANSI_CPP} -M ${CFLAGS} ${OBJECT_SRC} | sed s\|/usr/local/otf\|\$$\{OTFROOT\}\| >> Makefile
# more generally, use:
#	${ANSI_CPP} -M ${CFLAGS} $(INCLUDES) ${OBJECT_SRC} | sed s\|REPLACE_THIS\|\$WITH_THIS\| >> Makefile


# "make depend" requires that you...
# DO NOT MOVE OR DELETE (or place your own text after) THE FOLLOWING LINE:
### DEPEND
