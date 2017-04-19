#!/usr/bin/env python

import oat2
import os, glob

cedcfile = "raw/cedc_phase/ZJ.201409270830.0004.C.001"
cedcfile = "raw/cedc_phase/SC.201008180939.0002.C.001"
oat2file = ""
station2file = "china.station2"

files = glob.glob("raw/cedc_phase/*")

# count = 0
files.sort()
for cedcfile in files:

    if cedcfile[:3] == "CB":
        continue

    print cedcfile

    oat2.cedc2oat2(cedcfile=cedcfile, oat2file=oat2file, station2file=station2file, cal_dist=True, cal_traveltime=True, utctime=True)

    os._exit(-1)
    # count += 1

    # if count>10:
    #     os._exit(-1)