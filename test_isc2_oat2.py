#!/usr/bin/env python

import oat2
import os, glob


phase_file = "raw/isc_phase/1960s.txt"

files = glob.glob("raw/isc_phase/*.txt")
files.sort()

# split isc file
for phase_file in files:
    # if phase_file == "raw/isc_phase/1960s.txt":
    #     continue

    oat2.split_isc_phase(iscfile=phase_file, outdir="data/isc_phase")
# oat2.isc2oat(iscfile=phase_file, oatfile="", cal_dist=True, cal_traveltime=True)
os._exit(-1)




files = glob.glob("raw/isc_phase/*")

# count = 0
files.sort()
for cedcfile in files:

    if cedcfile[:3] == "CB":
        continue

    print cedcfile

    oat2.cedc2oat2(cedcfile=cedcfile, oat2file=oat2file, station2file=station2file, cal_dist=True, cal_traveltime=True)
    # count += 1

    # if count>10:
    #     os._exit(-1)