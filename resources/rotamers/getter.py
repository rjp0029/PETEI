#!/usr/bin/env python

import os
import sys

headerLine = "REMARK 220 SIDE CHAIN COORDINATES GENERATED FROM 2010 DUNBRACK ROTAMER LIBRARY\n"
headerLine += """REMARK 220 SHAPOVALOV, M.S. AND DUNBRACK, R.L., JR. (2011). A SMOOTHED 
REMARK 220 BACKBONE-DEPENDENT ROTAMER LIBRARY FOR PROTEINS DERIVED FROM ADAPTIVE
REMARK 220 KERNEL DENSITY ESTIMATES AND REGRESSIONS. STRUCTURE, 19, 844-858
REMARK 220 DOI: 10.1016/J.STR.2011.03.019\n"""

aminos = ['ala', 'cys', 'asp', 'glu', 'phe', 'gly', 'hsd', 'ile', 'lys', 'leu',
          'met', 'asn', 'gln', 'arg', 'ser', 'thr', 'val', 'trp', 'tyr']

SOURCE = "/home/rjp0029/PETEI/resources/rotlib/independ/"

for AA in aminos:
    for i in range(1, 201):
        fileName = AA + str(i) + ".pdb"
        try:
            f = open(SOURCE + fileName, "r")
        except IOError:
            break
        lines = f.readlines()
        f.close()
        atomNum = 0
        output = headerLine
        for line in lines:
            if line.startswith("REMARK"):
                continue
            elif line.startswith("END"):
                output += line.strip() + "\n"
            elif line.startswith("ATOM"):
                atomNum += 1
                text = "ATOM" + str(atomNum).rjust(7) + line[11:]
                output += text.strip()[:-3] + "ROT\n"
        f = open(fileName, "w")
        f.write(output)
        f.close()
