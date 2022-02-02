SHELL = /bin/sh
FAST_LOCAL_DIRECTORY = /usr/tmp/compile_temp

#INSTALL_PATH = .

# clear the list of suffixes
.SUFFIXES:

.SUFFIXES: .c .cc .o .a
.cc.o:; $(ANSI_CPP) $(CFLAGS) $(INCLUDES) $(DEFINES) $< 
.c.o:; $(ANSI_C) $(CFLAGS) $(INCLUDES) $(DEFINES) $< 

INCLUDES =

#Default target: create a library
libeigen.a: jacobi.o nrutil.o
	$(L_COMP) \
	libeigen.a \
	jacobi.o \
    nrutil.o

install:
	$(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" INSTALL_PATH="$(INSTALL_PATH)" libeigen.a
	-cp -f libeigen.a $(INSTALL_PATH)

GENERATED_FILES = *.o core lib*.a *.tmp

clean:
	rm -f $(GENERATED_FILES)

distclean:
	$(MAKE) clean
	rm -f $(INSTALL_PATH)/libeigen.a

OBJECT_SRC = jacobi.c nrutil.c

depend:
	mv Makefile Makefile.tmp
# The next line erases everything after the special "DO NOT MOVE...###"
# line below and copies it to Makefile.tmp
	sed -n '1,/DE\PEND/p' < Makefile.tmp > Makefile
# Generate the dependencies (using the compiler's -M option) and append
# them to the makefile
	$(ANSI_C) -M $(CFLAGS) $(INCLUDES) $(DEFINES) $(OBJECT_SRC) >> Makefile
	rm -f Makefile.tmp
# Here's the version from gregc's makefile
#	${ANSI_C} -M ${CFLAGS} $(INCLUDES) $(DEFINES) ${OBJECT_SRC} | sed s\|/usr/local/otf\|\$$\{OTFROOT\}\| >> Makefile
# more generally, use:
#	${ANSI_C} -M ${CFLAGS} $(INCLUDES) $(DEFINES) ${OBJECT_SRC} | sed s\|REPLACE_THIS\|\$WITH_THIS\| >> Makefile


# "make depend" requires that you...
# DO NOT MOVE OR DELETE (or place your own text after) THE FOLLOWING LINE:
### DEPEND
