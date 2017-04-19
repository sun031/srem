#!/usr/bin/env python

"""
This contains class oat2 and related functions

"""

import os
import glob
from obspy import UTCDateTime
from obspy.taup import TauPyModel
from obspy.geodetics.base import *

class oat2:

    def __init__(self):
        """
        Initialise the class.

        """
        ## event origin
        self.origin = UTCDateTime("1970-01-01")

        ## event latitude
        self.evla = -1.0
        self.evlo = -1.0
        self.evdp = -1.0
        self.mag = -1.0
        self.magt = "???"
        self.knetwk = "???"
        self.kstnm = "???"
        self.kcmpnm = "???"
        self.stla = -1.0
        self.stlo = -1.0
        self.stel = -1.0
        self.dist = -1.0
        self.gcarc = -1.0
        self.az = -1.0
        self.baz = -1.0
        self.phase = "???"
        self.tak135 = -1.0
        self.tobs = UTCDateTime("1970-01-01")

    def read_line(self, line):
        """
        Read a line in oat2 format

        :param self:
        :param line:
        :return:
        """
        row = line.split()
        self.origin = UTCDateTime(row[0])
        self.evla = float(row[1])
        self.evlo = float(row[2])
        self.evdp = float(row[3])
        self.mag = float(row[4])    ## \param mag magnitude
        self.magt = row[5].strip()
        self.knetwk = row[6].strip()
        self.kstnm = row[7].strip()
        self.kcmpnm = row[8].strip()
        self.stla = float(row[9])
        self.stlo = float(row[10])
        self.stel = float(row[11])
        self.dist = float(row[12])
        self.gcarc = float(row[13])
        self.az = float(row[14])
        self.baz = float(row[15])
        self.phase = row[16].strip()
        self.tak135 = float(row[17])
        self.tobs = UTCDateTime(row[18])

    def write(self, fp):
        """
        Write oat into file

        :param fp: file pointer
        :return:
        """
        fp.write('%s %12.4f %12.4f %8.3f %8.3f %-8s %-8s %-8s %-8s %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %-8s %12.4f %s\n'%(self.origin, self.evla, self.evlo, self.evdp, self.mag, self.magt, self.knetwk, self.kstnm, self.kcmpnm, self.stla, self.stlo, self.stel, self.dist, self.gcarc, self.az, self.baz, self.phase, self.tak135, self.tobs))



def read_station2(station_file):
    """
    Read station information from station2 file

    :param station_file:
    :return: stations, a dictionary
    """

    # read
    with open(station_file, "r") as fp:
        lst = fp.readlines()

    stations = {}
    for line in lst:
        row = line.split()
        key = row[0]
        lon = row[1]
        lat = row[2]
        el  = row[3]
        # key, lon, lat, el = line.split()[0:5]
        value = {
            "stlo" : float(lon),
            "stla" : float(lat),
            "stel" : float(el),
        }
        stations[key] = value

    return stations

def station2_info(key, stations):
    """
    read station information from stations dictionary

    :param key: name of station, e.g., "AH.ANQ"
    :param stations: a dictionary
    :return: stlo, stla, stel
    """
    dict = stations[key]
    stla = dict["stla"]
    stlo = dict["stlo"]
    stel = dict["stel"]

    return stlo, stla, stel


def traveltime_first_arrival(source_depth_in_km, distance_in_degree, phase):

    # calculate the theoretical traveltime
    model = TauPyModel(model="ak135")

    if phase[0] in ["P", "p"]:
        arrivals = model.get_travel_times(source_depth_in_km=source_depth_in_km, distance_in_degree=distance_in_degree,
                                          phase_list=["Pg", "Pn", "P", "p"])
    elif phase[0] in ["S", "s"]:
        arrivals = model.get_travel_times(source_depth_in_km=source_depth_in_km, distance_in_degree=distance_in_degree,
                                          phase_list=["Sg", "Sn", "S", "s"])
    else:
        return -1.0

    if len(arrivals)<1:
        print "no %sf arrivals available." % phase[0]
        t_pf = -1.0
        # name_pf = "-" % phase[0]
    else:
        tlst = []
        for arr in arrivals:
            tlst.append(arr.time)
        maxx = min(tlst)
        imax = tlst.index(maxx)

        arr = arrivals[imax]
        t_pf = arr.time
        # name_pf = arr.name

    return t_pf

def cedc2oat2(cedcfile, station2file, oat2file="", cal_dist=False, cal_traveltime=False, utctime=True):
    """
    Convert CEDC phase into oat2 format.

    :param cedcfile: Note each cedcfile contains only one event
    :param oat2file: if oat2file == "", then save the oatfile under data/oat/cedcfile.oat2
    :param station2file:
    :param cal_dist: if False, only have event and station information no dist az information
                 if True, calculate dist, az baz, gcarc and write into oat2 file
    :param cal_traveltime:  if False, it would not calculate the first arrival traveltime
    :param utctime: the original CEDC phase gives Beijing Time, if True, then convert Beijing Time to UTCTime

    :return: None
    """

    # read CEDC phase file
    with open(cedcfile, "r") as fp:
        phlst = fp.readlines()

    # read stations
    stations = read_station2(station_file=station2file)


    oat = oat2()

    # read event infomation
    for line in phlst:
        row = line.split()
        if line[:3] != "DBO":
            continue


        date = row[2]
        time = row[3]
        oat.origin = UTCDateTime(date + "T" + time)

        oat.evla = float(row[4])
        oat.evlo = float(row[5])
        oat.evdp = float(row[6])

        oat.magt = row[7]
        oat.mag  = float(row[8])

    if utctime:
        oat.origin = oat.origin - 8.0 * 3600.0

    if oat2file == "":
        bname = oat.origin.datetime.strftime("%Y%m%d%H%M%S")
        try:
            os.makedirs("data/oat_cedc")
        except:
            pass

        oat2file = "data/oat_cedc/" + bname + ".oat2"

    # open oatfile
    f = open(oat2file, "w")

    # read phase information
    for line in phlst:
        row = line.split()
        if line[:3] != "DPB":
            continue


        oat.knetwk = line[3:7].strip()
        oat.kstnm = line[6:13].strip()
        oat.kcmpnm= line[12:17].strip()
        oat.phase = line[20:29].strip()

        if "M" in oat.phase:
            continue

        date = line[33:45].strip()
        time = line[44:57].strip()
        # print date, time

        try:
            oat.tobs = UTCDateTime(date + "T" + time)
        except:
            continue


        try:
            key = oat.knetwk +"." + oat.kstnm
            stlo, stla, stel = station2_info(key, stations)
        except:
            print "Cannot find station", key, oat.origin
            continue

        oat.stlo = stlo
        oat.stla = stla
        oat.stel = stel

        if cal_traveltime:
            cal_dist = True

        if cal_dist:
            # calculate distance for traveltime
            dist, az, baz = gps2dist_azimuth(oat.evla, oat.evlo, oat.stla, oat.stlo)
            oat.dist = dist/1000.0
            oat.gcarc = kilometer2degrees(oat.dist)
            oat.az   = az
            oat.baz  = baz
            # print oat.knetwk, oat.kstnm, oat.kcmpnm, oat.phase, oat.tobs, oat.az, oat.baz

        if cal_traveltime:
            oat.tak135 = traveltime_first_arrival(source_depth_in_km=oat.evdp, distance_in_degree=oat.gcarc, phase=oat.phase[0])

        if utctime:
            oat.tobs = oat.tobs - 8.0 * 3600.0

        oat.write(f)


    f.close()

def split_isc_phase(iscfile, outdir="data/isc_phase"):
    """
    Split ISC phase containing numerous events into ISC phase containing one event for one file.

    :param iscfile:
    :param outdir: default is "data/isc_phase", the output filename should be origin.txt
    :return:
    """

    # open ISC phase file
    with open(iscfile, "r") as fp:
        isclst = fp.readlines()

    try:
        os.makedirs(outdir)
    except:
        pass

    namelst = []

    for line in isclst:
        # print line
        row = line.split(',')

        try:
            evid = int(row[0].strip())
        except:
            continue

        date = row[18].strip()
        time = row[19].strip()
        origin = UTCDateTime(date+"T"+time)

        filename = outdir + "/" + origin.datetime.strftime("%Y%m%d%H%M%S") + ".txt"
        fp = open(filename, "a")
        fp.write(line)
        fp.close()

        if filename not in namelst:
            namelst.append(filename)
            print filename

    for name in namelst:
        os.system("sort -u %s > temp" % (name))
        os.system("mv temp %s" % name)


def isc2oat(iscfile, oatfile="", cal_dist=False, cal_traveltime=False):
    """
    Convert ISC phase into oat2 format, each iscfile contain one event.

    :param iscfile:
    :param oatfile:
    :param cal_dist:
    :param cal_traveltime:
    :return:
    """

    # open ISC phase file
    with open(iscfile, "r") as fp:
        isclst = fp.readlines()

    # geting oat2 filename
    if oatfile == "":
        try:
            os.makedirs("data/oat_isc")
        except:
            pass
        bname = os.path.basename(iscfile).replace(".txt", ".oat2")
        oname = "data/oat_isc/" + bname
    else:
        oname = oatfile

    print oname

    fp = open(oname, "w")
    oat = oat2()

    for line in isclst:
        # print line
        row = line.split(',')

        try:
            evid = int(row[0].strip())
        except:
            continue

        oat.__init__()
        oat.kstnm = row[2].strip()
        oat.stla = float(row[3].strip())
        oat.stlo = float(row[4].strip())
        oat.stel = float(row[5].strip())
        oat.kcmpnm = row[6].strip()

        oat.gcarc = float(row[7].strip())
        oat.dist = oat.gcarc*111.195
        oat.baz = float(row[8].strip())
        oat.phase = row[10].strip()
        if oat.phase == "":
            continue

        if oat.knetwk == "???":
            oat.knetwk = "ISC"

        date = row[11].strip()
        time = row[12].strip()
        oat.tobs = UTCDateTime(date+"T"+time)

        # res = float(row[13].strip())

        date = row[18].strip()
        time = row[19].strip()
        oat.origin = UTCDateTime(date+"T"+time)

        # oat.tobs = tobs - orgi
        oat.tak135 = -1.0


        try:
            oat.evdp = float(row[22].strip())
            oat.evlo = float(row[21].strip())
            oat.evla = float(row[20].strip())
        except:
            continue

        try:
            oat.mag = float(row[25].strip())
        except:
            oat.mag = -1.0

        if row[24].strip() != "":
            oat.magt = row[24].strip()

        if cal_traveltime:
            cal_dist = True

        if cal_dist:
            # calculate distance for traveltime
            dist, az, baz = gps2dist_azimuth(oat.evla, oat.evlo, oat.stla, oat.stlo)
            oat.dist = dist/1000.0
            oat.gcarc = kilometer2degrees(oat.dist)
            oat.az   = az
            oat.baz  = baz
            # print oat.knetwk, oat.kstnm, oat.kcmpnm, oat.phase, oat.tobs, oat.az, oat.baz

        if cal_traveltime:
            oat.tak135 = traveltime_first_arrival(source_depth_in_km=oat.evdp, distance_in_degree=oat.gcarc, phase=oat.phase[0])

        oname = "data/oat_isc/" + oat.origin.datetime.strftime("%Y%m%d%H%M%S") + ".oat2"
        # oname = oat.origin.datetime.strftime("%Y%m%d%H%M%S")
        # print oname

        oat.write(fp)

    fp.close()

def oat2_station2_event2(oat2_file, station2_file="", event2_file=""):

    with open(oat2_file, "r") as fp:
        lst = fp.readlines()

    evtlst = []
    stalst = []

    oat = oat2()
    for line in lst:
        oat.read_line(line)

        evt = [oat.origin, oat.evlo, oat.evla, oat.evdp, oat.mag , oat.magt]
        if evt not in evtlst:
            evtlst.append(evt)

        sta = [oat.knetwk+"."+oat.kstnm, oat.stlo, oat.stla, oat.stel, -12345.0]
        if sta not in stalst:
            stalst.append(sta)

    evtlst.sort()
    stalst.sort()

    if station2_file == "":
        station2_file = oat2_file.replace(".oat2", ".station2")

    if event2_file == "":
        event2_file = event2_file.replace(".oat", ".event2")


    fp = open(station2_file, "w")
    for line in stalst:
        fp.write("%s\t%12.4f\t%12.4f\t%12.4f\t%12.4f\n" % (line[0], line[1], line[2], line[3], line[4]))
    fp.close()

    fp = open(event2_file, "w")
    for line in event2_file:
        fp.write("%s\t%12.4f\t%124.f\t%12.4f\t%12.4f\t%s\n" % (line[0], line[1], line[2], line[3], line[4], line[5]))
    fp.close()

def merge_isc_cedc(cedc_file, isc_dir, origin_error=60, location_error=0.5):

    # get cedc origin time
    filename = os.path.basename(cedc_file).split(".")[0]
    print filename
    cedc_origion = UTCDateTime(filename)

    # time windows
    t1 = cedc_origion - origin_error
    t2 = cedc_origion + origin_error

    a1 = cedc_origion - 3600.0
    a2 = cedc_origion + 3600.0

    n1 = a1.datetime.strftime("%Y%m%d%H")
    n2  = cedc_origion.datetime.strftime("%Y%m%d%H")
    n3 = a2.datetime.strftime("%Y%m%d%H")
    # print n1, n2, n3

    file1 = glob.glob(isc_dir+"/"+n1+"*.oat2")
    file2 = glob.glob(isc_dir+"/"+n2+"*.oat2")
    file3 = glob.glob(isc_dir+"/"+n3+"*.oat2")

    # print isc_dir+"/"+n1+"*.oat"
    # ISC events in the time window
    file = file1 + file2 + file3
    files = list(set(file))
    files.sort()

    if len(files)<1:
        print "Same event not found", cedc_file


    print len(files)

    # find events according to the orgin error
    tlst = []
    for file in files:
        filename = os.path.basename(file).split(".")[0]
        tfile = UTCDateTime(filename)

        if tfile>=t1 and tfile<=t2:
            tlst.append(file)

    print len(tlst)

    if len(tlst)<1:
        return

    with open(cedc_file, "r") as fp:
        line = fp.readline()

    print line

    oat = oat2()
    oat.read_line(line)
    evla_cedc = oat.evla
    evlo_cedc = oat.evlo

    # print evla_cedc, evlo_cedc

    # select event using location error
    samefile = []
    for iscfile in tlst:

        with open(iscfile, "r") as fp:
            line = fp.readline()

        oat.read_line(line)
        evla_isc = oat.evla
        evlo_isc = oat.evlo

        dist, az, baz = gps2dist_azimuth(evla_cedc, evlo_cedc, evla_isc, evlo_isc)
        dist = dist/1000.0
        gcarc = kilometer2degrees(dist)
        if gcarc > location_error:
            continue

        samefile.append(iscfile)

    if len(samefile)==1:
        print samefile[0]
    else:
        print "No same file found.", cedc_file
        return

    # merge event
    isc_file = samefile[0]
    with open(isc_file, "r") as fp:
        isc_lst = fp.readlines()

    with open(cedc_file, "r") as fp:
        cedc_lst = fp.readlines()
















    pass