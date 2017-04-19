#!/usr/bin/env python

import oat2
import os, glob
from multiprocessing import Pool


oat2file = ""
station2file = "china.station2"

def p_cedc2oat2(cedcfile):
    oat2.cedc2oat2(cedcfile=cedcfile, oat2file=oat2file, station2file=station2file, cal_dist=True, cal_traveltime=True, utctime=True)


files2 = glob.glob("raw/cedc_phase/*")
files  = []
for line in files2:
    bname = os.path.basename(line)
    if bname[:2] != "CB":
        files.append(line)

# print len(files2), len(files)

# count = 0
files.sort()

# Make the Pool of workers
pool = Pool()
pool.map(p_cedc2oat2, files)

# close the pool and wait for the work to finish
pool.close()
pool.join()






# for cedcfile in files:
#
#     if cedcfile[:3] == "CB":
#         continue
#
#     print cedcfile
#
#     oat2.cedc2oat2(cedcfile=cedcfile, oat2file=oat2file, station2file=station2file, cal_dist=True, cal_traveltime=True)
    # count += 1

    # if count>10:
    #     os._exit(-1)