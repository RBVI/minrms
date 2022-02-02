INSTALL_PATH = .
FILES_TO_COPY = \
README.txt \
atoms_used.txt \
res_code_dict.txt

default:
	echo "WARNING: No default behavior for copying shared files.\n         Must specify destination path."

install:
	$(MAKE) INSTALL_PATH="$(INSTALL_PATH)" install_public_libraries
	-cp -f $(FILES_TO_COPY) $(INSTALL_PATH)/share

install_public_libraries:
	-mkdir $(INSTALL_PATH)/share

clean:
#	echo "WARNING: No default behavior for deleting local shared files."
#	echo "Must specify destination path."

distclean:	clean

depend:
	echo "No dependencies for shared files."
