#!/usr/bin/env python
#
# Copyright (c) 2002 The Regents of the University of California.
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#   1. Redistributions of source code must retain the above copyright
#      notice, this list of conditions, and the following disclaimer.
#   2. Redistributions in binary form must reproduce the above
#      copyright notice, this list of conditions, and the following
#      disclaimer in the documentation and/or other materials provided
#      with the distribution.
#   3. Redistributions must acknowledge that this software was
#      originally developed by the UCSF Computer Graphics Laboratory
#      under support by the NIH National Center for Research Resources,
#      grant P41-RR01081.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
# OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
# OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import string
import pprint

_ThreeOneMap = {
    'ALA':'A',
    'ARG':'R',
    'ASN':'N',
    'ASP':'D',
    'CYS':'C',
    'GLU':'E',
    'GLN':'Q',
    'GLY':'G',
    'HIS':'H',
    'HYP':'O',
    'ILE':'I',
    'LEU':'L',
    'LYS':'K',
    'MET':'M',
    'PHE':'F',
    'PRO':'P',
    'SER':'S',
    'THR':'T',
    'TRP':'W',
    'TYR':'Y',
    'VAL':'V',
    'ASX':'B',
    'GLX':'Z',
}

def readAlign(fname):
    aList = []
    f = open(fname)
    while True:
        line = f.readline()
        if not line:
            break
        if line[:5] == 'EQUIV':
            break
    while True:
        line = f.readline()
        if not line:
            break
        v = line.split()
        if len(v) != 10:
            continue
        aList.append(((int(v[1]), v[2], int(v[3])),
                      (int(v[7]), v[8], int(v[9]))))
    f.close()
    return aList

def readPDB(fname):
    aList = []
    f = open(fname)
    while True:
        line = f.readline()
        if not line:
            break
        if line[:4] != 'ATOM' or line[12:16] != ' CA ':
            continue
        aList.append((line[17:20], int(line[22:26])))
    f.close()
    return aList

def verify(align, pdb1, pdb2):
    for a in align:
        (ai, seqai, i), (aj, seqaj, j) = a
        r1, s1 = pdb1[i - 1]
        if r1 != seqai:
            print('PDB1 mismatch:', i, seqai, r1)
        r2, s2 = pdb2[j - 1]
        if r2 != seqaj:
            print('PDB2 mismatch:', j, seqai, r2)

def makeLine(part):
    fragments = []
    for i in range(0, len(part), 10):
        fragments.append(''.join(part[i:i + 10]))
    return ' '.join(fragments)

def generate(align, pdb1, pdb2, n1, n2):
    print('..')
    print('//')
    index1 = 1
    index2 = 1
    list1 = []
    list2 = []
    for a in align:
        (ai, seqai, i), (aj, seqaj, j) = a
        while index1 < i:
            list1.append(_ThreeOneMap[pdb1[index1 - 1][0]])
            list2.append('.')
            index1 = index1 + 1
        while index2 < j:
            list1.append('.')
            list2.append(_ThreeOneMap[pdb2[index2 - 1][0]])
            index2 = index2 + 1
        list1.append(_ThreeOneMap[seqai])
        index1 = index1 + 1
        list2.append(_ThreeOneMap[seqaj])
        index2 = index2 + 1
    while index1 <= len(pdb1):
        list1.append(_ThreeOneMap[pdb1[index1 - 1][0]])
        list2.append('.')
        index1 = index1 + 1
    while index2 <= len(pdb2):
        list1.append('.')
        list2.append(_ThreeOneMap[pdb2[index2 - 1][0]])
        index2 = index2 + 1
    format = '%%%ds  %%s' % max(len(n1), len(n2))
    for i in range(0, len(list1), 50):
        part1 = list1[i:i + 50]
        part2 = list2[i:i + 50]
        print('\n')
        #print format % (n2, makeLine(part2))
        #print format % (n2, makeLine(part2))
        print(format % (n1, makeLine(part1)))
        print(format % (n2, makeLine(part2)))

def main():
    import sys
    if len(sys.argv) != 4:
        print('Usage:', sys.argv[0], 'align_file PDB_file1 PDB_file2')
        raise(SystemExit, 0)
    align = readAlign(sys.argv[1])
    pdb1 = readPDB(sys.argv[2])
    pdb2 = readPDB(sys.argv[3])
    verify(align, pdb1, pdb2)
    generate(align, pdb1, pdb2, sys.argv[2], sys.argv[3])

if __name__ == '__main__':
    main()
