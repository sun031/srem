#!/usr/bin/env python

import oat2
import os, glob
from multiprocessing import Pool


files = glob.glob("data/isc_phase/*.txt")
files.sort()

# a signal file test
# count = 0
# files.sort()
# for file in files:
#
#
#     oat2.isc2oat(iscfile=file, oatfile="", cal_dist=True, cal_traveltime=True)
#
#     count += 1
#
#     if count>0:
#         os._exit(-1)


def p_isc2oat(file):
    oat2.isc2oat(iscfile=file, oatfile="", cal_dist=True, cal_traveltime=True)

# Make the Pool of workers
pool = Pool()
pool.map(p_isc2oat, files)

# close the pool and wait for the work to finish
pool.close()
pool.join()
