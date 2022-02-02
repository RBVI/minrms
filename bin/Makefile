SHELL = /bin/sh

#default parameters
# for SGI:
#ANSI_C = cc
#ANSI_CPP = CC
#L_COMP = CC -ar -o
#CFLAGS = -c -LANG:ansi-for-init-scope=on
#LFLAGS = 
# for Alpha:
#FAST_LOCAL_DIRECTORY = /usr/tmp/minrms_compile_temp
#ANSI_C = cc
#ANSI_CPP = cxx -ptr $(FAST_LOCAL_DIRECTORY) -std strict_ansi -DCOMPILER_BUG1
#L_COMP = ar -r
#CFLAGS = -c
#LFLAGS = 


#The "SRC_DIRS" directories contain source code for different binaries
#including the "minrms" binary.
#This makefile constructs each binary from the code in these directories
#and the libraries in the "LIB_DIRS" directories.
SRC_DIRS = \
align2msf \
fssp2msf \
minrms \
msf2sequence \
msf2stat3d \
msf_compare \
pdb2sequence \
rotate_pdb \
pdb_select


default:
#create all binaries and leave them in thier own directories.
	for i in $(SRC_DIRS) ; do \
		(cd $$i; \
         $(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" INSTALL_PATH=$(INSTALL_PATH) install \
		) ;\
	done

#When the "install" target is used, some of the libraries
#get copied to a directory designated by the INSTALL_PATH variable.
install:
#1)Create the directory where the binaries will be stored.
	$(MAKE) INSTALL_PATH="$(INSTALL_PATH)" install_public_directories
#2)install binaries in that directory
	for i in $(SRC_DIRS) ; do \
		(cd $$i; \
         $(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" L_COMP="$(L_COMP)" LFLAGS="$(LFLAGS)" INSTALL_PATH="$(INSTALL_PATH)/bin" install ); \
	done

depend:
	for i in $(SRC_DIRS) ; do \
	     (cd $$i; \
	     $(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" depend ); \
	done

install_public_directories:
	-mkdir $(INSTALL_PATH)/bin

clean:
	for i in $(SRC_DIRS) ; do \
		(cd $$i; $(MAKE) clean) ;\
	done

distclean:	clean
