#!/usr/bin/env python

import argparse
import re
import os.path
import sys
import subprocess

parser = argparse.ArgumentParser(description="Hello, this tool maps values from the CHAIN program output to a PDB file viewable in Pymol. "
                                             "You will need to convert the RTF to TXT using 'unoconv'. "
                                             "You can do this via the commands: 'sudo apt-get install unoconv' "
                                             "and 'unoconv -f .txt mychainoutput.rtf'. "
                                             "PLEASE DOUBLE CHECK YOUR RESULTS! rtf conversion can be an unreliable process. "
                                             "Example: main.py -g f -m c /home/me/Documents/chain2pymol/chainOut.txt /home/me/Documents/chain2pymol/myPymol.pdb HUMAN /home/chan/Documents/chain2pymol/")
parser.add_argument("chainIn", help="specify absolute path and name for input file .txt")
parser.add_argument("pdbIn", help="specify absolute path and name for input .pdb")
parser.add_argument("target", help="a keyword contained in the sequence in the rtf that matches the pdb, your query (case sensitive)")
parser.add_argument("output", help="specify absolute path but not the filename for output pdb")
parser.add_argument("-g", "--group", choices=['f','b'], help="select from foreground or background group")
parser.add_argument("-m", "--metric", choices=['c','i','d'], help="select conservation, insertion, or deletion")
args = parser.parse_args()

chainTxt = open(args.chainIn, "r")
lines = chainTxt.readlines()

WeightCounter = 0
InsertCounter = 0
DeleteCounter = 0
TargetCounter = 0

BackWeight = ""
BackInsert = ""
BackDelete = ""
ForeWeight = ""
ForeInsert = ""
ForeDelete = ""

TargetName = args.target
TargetSequence = ""

for index, line in enumerate(lines):
    if "wt_res_freqs" in line:
        WeightCounter += 1
        if WeightCounter % 2 == 0:
            ForeWeight += line[21:]
        else: BackWeight += line[21:124] + "\n"
    if "insertions" in line:
        InsertCounter += 1
        if InsertCounter % 2 == 0:
            ForeInsert += line[12:]
        else: BackInsert += line[12:]
    if "deletions" in line:
        DeleteCounter += 1
        if DeleteCounter % 2 == 0:
            ForeDelete += line[11:114]
        else: BackDelete += line[11:114]
    if TargetName in line:
        TargetCounter += 1
        if TargetCounter % 3 == 0:
            TargetSequence += line[17:124] + "\n"

lengthTotal = len(BackWeight) - 1

BackWeight = BackWeight.replace(" ", "0").replace("\n", "").replace("\t", "")
ForeWeight = ForeWeight.replace(" ", "0").replace("\n", "").replace("\t", "")
BackInsert = BackInsert.replace(" ", "0").replace("\n", "").replace("\t", "")
ForeInsert = ForeInsert.replace(" ", "0").replace("\n", "").replace("\t", "")
BackDelete = BackDelete[:lengthTotal].replace(" ", "0").replace("\n", "").replace("\t", "")
ForeDelete = ForeDelete[:lengthTotal].replace(" ", "0").replace("\n", "").replace("\t", "")
TargetSequence = re.sub(r'\d', "", TargetSequence).replace("\t", "").replace("\n", "")[:-1]

if args.group == "f":
    if args.metric == "c":
        DesiredMetric = ForeWeight
    elif args.metric == "i":
        DesiredMetric = ForeInsert
    elif args.metric == "d":
        DesiredMetric = ForeDelete
    else:
        print "incorrect metric argument, see --help"
elif args.group == "b":
    if args.metric == "c":
        DesiredMetric = BackWeight
    elif args.metric == "i":
        DesiredMetric = BackInsert
    elif args.metric == "d":
        DesiredMetric = BackDelete
    else:
        print "incorrect metric argument, see --help"
else:
    print "incorrect group argument, see --help"

combined = zip(TargetSequence, DesiredMetric)
filtered = filter(lambda (a, b): a != '-' or '*' or '.', combined)
TargetSequence, DesiredMetric = zip(*filtered)

pdbFile = open(args.pdbIn, "r")
pdbLines = pdbFile.read().splitlines()
StartIndex = 0
StopIndex = 0
TotalBlockLength = 0

for pdbIndex, pdbline in enumerate(pdbLines):
    if pdbline.startswith("ATOM"):
        StartIndex = pdbIndex
        break
for pdbIndex, pdbline in enumerate(pdbLines[StartIndex:]):
    if pdbline.startswith("ATOM"):
        continue
    else: StopIndex = pdbIndex + StartIndex
    break

for pdbIndex, pdbline in enumerate(pdbLines[StartIndex:]):
    if pdbline.startswith("ATOM") or pdbline.startswith("TER"):
        continue
    else: TotalBlockLength = pdbIndex
    break

AtomBlock = pdbLines[StartIndex:TotalBlockLength + StartIndex]
Interval = StopIndex - StartIndex
currentResidue = 0
newAtomBlock = ""

for atomIndex, atomline in enumerate(AtomBlock):
    if "TER" in atomline:
        newAtomBlock += atomline + "\n"
    else:
        currentResidue = int(atomline[23:26].replace(" ", ""))
        newAtomBlock += atomline[:61] + DesiredMetric[currentResidue - 1] + "0.00" + atomline[66:78] + "\n"


toFile = "\n".join(pdbLines[:StartIndex]) + "\n" + newAtomBlock + "\n".join(pdbLines[TotalBlockLength + StartIndex:])

f = open(os.path.join(args.output, os.path.basename('toolOut.pdb')), 'w')
f.write(toFile)
f.close()

print "complete!"