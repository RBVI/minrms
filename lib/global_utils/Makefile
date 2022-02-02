SHELL = /bin/sh

FAST_LOCAL_DIRECTORY = /usr/tmp/compile_temp

INSTALL_PATH = .

# clear the list of suffixes
.SUFFIXES:

.SUFFIXES: .c .cc .o .a
.cc.o:; $(ANSI_CPP) $(CFLAGS) $(INCLUDES) $(DEFINES) $< 


LIBS = 


INCLUDES = 

#Default target: create a library
libglobal_utils.a: \
stream_bugfixes.o
	$(CPP_PRELINKER_COMMAND)
	$(L_COMP) \
	libglobal_utils.a \
	stream_bugfixes.o

#	$(LINKER_TEMP_FILES)


install:
	$(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" INSTALL_PATH="$(INSTALL_PATH)" libglobal_utils.a
	-cp -f libglobal_utils.a $(INSTALL_PATH)

RCS_FILES = \
DEBUG_MSG.h \
ERR.h \
global_utils.h \
file_constants.h \
binary_file.h \
short_string.h \
simple_numeric_types.h \
simple_numeric_utils.h \
stl_utils.h \
stream_bug_fixes.h \
stream_bug_fixes.cc

checkout:
	co -q $(RCS_FILES)

GENERATED_FILES = *.o *.a core $(COMPILER_TEMP_FILES) $(LINKER_TEMP_FILES)

clean:
	rm -f $(GENERATED_FILES)

distclean:
	$(MAKE) clean
	rm -f $(INSTALL_PATH)/libglobal_utils.a


OBJECT_SRC = stream_bug_fixes.cc

depend:
	mv Makefile Makefile.tmp
# The next line erases everything after the special "DO NOT MOVE...###"
# line below and copies it to Makefile.tmp
	sed -n '1,/DE\PEND/p' < Makefile.tmp > Makefile
# Generate the dependencies (using the compiler's -M option) and append
# them to the makefile
	$(ANSI_CPP) -M $(CFLAGS)  $(INCLUDE) $(OBJECT_SRC) >> Makefile
	rm -f Makefile.tmp
# Here's the version from gregc's makefile
#	${ANSI_CPP} -M ${CFLAGS} $(INCLUDE) ${OBJECT_SRC} | sed s\|/usr/local/otf\|\$$\{OTFROOT\}\| >> Makefile
# more generally, use:
#	${ANSI_CPP} -M ${CFLAGS} $(INCLUDE) ${OBJECT_SRC} | sed s\|REPLACE_THIS\|\$WITH_THIS\| >> Makefile


# "make depend" requires that you...
# DO NOT MOVE OR DELETE (or place your own text after) THE FOLLOWING LINE:
### DEPEND
