#!/usr/bin/env bash

test_simple() {
    cd tests/
    wget http://www.rcsb.org/pdb/files/4hhb.pdb
    wget http://www.rcsb.org/pdb/files/1a6m.pdb
    mv 4hhb.pdb hemoglobin.pdb
    mv 1a6m.pdb myoglobin.pdb
    rm -f index*.html
    ../bin/minrms/minrms -fm 4 -of 8.0,0.33 -HS \
                         -minN 0.33 -max-rmsd 3.5 -ir -r \
                         hemoglobin.pdb,"*.A" myoglobin.pdb

    assertTrue "align141.msf file not created" "[ -s align141.msf ]"
    rm -f align*.msf align_chimera.*
    #rm -f *.pdb  (commenting out. we might need these files later.)
  cd ../
}

. shunit2/shunit2
