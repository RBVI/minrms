#!/usr/bin/env bash

test_simple() {
    cd tests/
    # commenting out:
    # (we already downloaded these files in an earlier test)
    #wget http://www.rcsb.org/pdb/files/4hhb.pdb
    #wget http://www.rcsb.org/pdb/files/1a6m.pdb
    #mv 4hhb.pdb hemoglobin.pdb
    #mv 1a6m.pdb myoglobin.pdb
    # ...but let's rename them back to their original names:
    mv hemoglobin.pdb 4hhb.pdb 
    mv myoglobin.pdb 1a6m.pdb 
    python3 ../bin/align2msf/align2msf.py align_server_output/1a6m+4hhbA.txt 1a6m.pdb 4hhb.pdb > align2msf_output.msf
    FILE_SIZE=`wc -c align2msf_output.msf | awk '{print $1}'`
    assertTrue "align2msf.py failed" "[ -eq $FILE_SIZE 1578 ]"
    #rm -f align2msf_output.msf
  cd ../
}

. shunit2/shunit2
