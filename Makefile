SHELL = /bin/sh

INSTALL_PATH = /usr/local/minrms

#The "SRC_DIRS" directories contain source code for different binaries
#including the "minrms" binary.
#This makefile constructs each binary from the code in these directories
#and the libraries in the "LIB_DIRS" directories.
SRC_DIRS = \
	lib \
	bin \
	share

install:
	#create all libraries and all binaries as well as the directories they
	#need to reside in, and copy them there.
	-mkdir $(INSTALL_PATH)
	for i in $(SRC_DIRS) ; do \
		(cd $$i; \
		$(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" \
			L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" \
			LFLAGS="$(LFLAGS)" INSTALL_PATH=$(INSTALL_PATH) \
			install \
		);\
	done
	echo "\n  ********  MINRMS compilation finished.  ********\n" \
	     "  Check the ${INSTALL_PATH}/bin directory to\n" \
	     "  see whether it contains a copy of the \"minrms\" binary.\n" \
             "  (Other optional binaries might be included there too.)\n" \
	     "  If not, then compilation was not successful.)"

depend:
	for i in $(SRC_DIRS) ; do \
		(cd $$i; \
		$(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" \
			L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" \
			depend \
		); \
	done

install_public_directories:
	-mkdir $(INSTALL_PATH)

clean:
	for i in $(SRC_DIRS) ; do \
		(cd $$i; $(MAKE) clean) ;\
	done

distclean:
	for i in $(SRC_DIRS) ; do \
		(cd $$i; \
		 $(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" \
			L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" \
			LFLAGS="$(LFLAGS)" INSTALL_PATH=$(INSTALL_PATH) \
			distclean \
		) ;\
	done

distribution:	distclean
	gnutar --gzip --create --exclude=RCS --exclude=not_implemented_yet \
		--directory=.. --file=/usr/tmp/conrad/minrms.tar.gz minrms
