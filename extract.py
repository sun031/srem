#!/usr/bin/env python


"""
Extract wavefroms from MINISEED and save as SAC format for events

This now works for China Reference Model.

"""
import os
import glob

from obspy import UTCDateTime
from obspy.taup import TauPyModel
from obspy.geodetics.base import *
from obspy import read
from obspy.io.sac import SACTrace

def read_event2(file):
    """
    :param file: file name of catalogue file
    :return: event list, [origin, evlo, evla, evdp, mag]
    """

    with open(file, "r") as fp:
        lst = fp.readlines()

    events = []
    for line in lst:
        origin, long, lat, dep, mag, magt = line.split()
        origin  = UTCDateTime(origin)
        evlo = float(long)
        evla = float(lat)
        evdp = float(dep)
        mag = float(mag)

        e = [origin, evlo, evla, evdp, mag]
        if e not in events:
            events.append(e)

    return events

def read_events(event_file):
    with open(event_file, "r") as fp:
        lst = fp.readlines()

    events = {}
    for line in lst:
        key, long, lat, dep, mag, magt = line.split()
        value = {
            "evlo" : float(long),
            "evla" : float(lat),
            "evdp" : float(dep),
            "mag" : float(mag),
            "magt" : magt
        }
        events[key] = value

    return events

def event_info(key, events):
    """
    read event information from a dictionary type variable
    :param key:
    :param events:
    :return:
    """

    dict = events[key]
    evla = dict["evla"]
    evlo = dict["evlo"]
    evdp = dict["evdp"]
    mag  = dict["mag"]

    return evlo, evla, evdp, mag

def read_event_line(e):
    """
    a test function
    """
    orig = e[0]
    evlo = e[1]
    evla = e[2]
    evdp = e[3]
    mag  = e[4]
    return orig, evlo, evla, evdp, mag

def read_stations(station_file):

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

def station_info(key, stations):
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

def cal_travel_time(evla, evlo, evdp, stla, stlo):

    # calculate distance for traveltime
    dist, az, baz = gps2dist_azimuth(evla, evlo, stla, stlo)

    # convert meter to km
    dist /= 1000.0
    gcarc = kilometer2degrees(dist)
    # print dist, az, baz, gcarc


    # calculate the theoretical traveltime
    model = TauPyModel(model="ak135")
    arrivals_pf = model.get_travel_times(source_depth_in_km=evdp, distance_in_degree=gcarc, phase_list=["Pg", "Pn", "P", "p"])
    arrivals_sf = model.get_travel_times(source_depth_in_km=evdp, distance_in_degree=gcarc, phase_list=["Sg", "Sn", "S", "s"])

    if len(arrivals_pf)<1:
        print "no Pf arrivals available."
        t_pf = -12345.0
        name_pf = "-P"
    else:
        tlst = []
        for arr in arrivals_pf:
            tlst.append(arr.time)
        maxx = max(tlst)
        imax = tlst.index(maxx)

        arr = arrivals_pf[imax]
        t_pf = arr.time
        name_pf = arr.name
    # print arr.time, arr.name, arr.source_depth, arr.distance, gcarc

    if len(arrivals_sf)<1:
        print "no Sf arrivals available."
        t_sf = -12345.0
        name_sf = "-S"
    else:
        tlst = []
        for arr in arrivals_sf:
            tlst.append(arr.time)
        maxx = max(tlst)
        imax = tlst.index(maxx)

        arr = arrivals_sf[imax]
        t_sf = arr.time
        name_sf = arr.name

    return dist, gcarc, az, baz, t_pf, name_pf, t_sf, name_sf


# def get_files(path, suffix=".mseed"):
#
#     files = glob.glob(path)
#     pass

def read_data(file, file_merge=""):
    """
    read raw data
    :param filename:
    :return:
    """
    st = read(file)
    # print st

    if file_merge=="":
        pass
    else:
        st_merge = read(file_merge)
        st += st_merge

    st2 = st.copy()
    try:
        st2.merge()
    except:
        pass
    try:
        st2 = st2.split()
    except:
        pass
    tr = st2[0]
    tr.detrend(type="linear")
    return tr

def get_event_waveform(key, events, stations, datapath, outpath, suffix=".mseed", origin=False, tmin=-200, tmax=1200):
    """
    get wavefrom for an event, the output format is SAC.
    :param key: origin time in string of UTCDateTime, 2016-09-26T15:48:50.000000Z
    :param events: a tuple containing event information
    :param stations: a tuple containing station information
    :param datapath: datapath of miniseed
    :param outpath: path saving SAC file
    :param suffix: the suffix of raw data in the datapath, here is ".mseed"
    :param origin: True or False, if True cut data from origin, else using the P first arrival as reference time
    :param tmin: relative time
    :param tmax: relative time
    :return: None
    """


    # event information
    evlo, evla, evdp, mag = event_info(key, events)
    print "Event", key, evlo, evla, evdp, mag

    # station information

    # check mseed existed or not
    # files = glob.glob(datapath+"/*.mseed")

    # make directory for saving waveforms
    o = UTCDateTime(key)
    o_bjt = o + 8.0 * 3600
    temp = o.datetime
    temp = temp.strftime("%Y%m%d%H%M%S")
    prefix = outpath+"/"+temp
    # print prefix
    try:
        os.makedirs(prefix)
    except:
        pass

    count = 0
    for key_st in stations.keys():
        # files = glob.glob(datapath+"/*%s*%s" % (key_st, suffix))
        # print files
        #
        # # in this case, file not exist.
        # if len(files)<1:
        #     continue

        count += 1

        # station information
        stlo, stla, stel = station_info(key=key_st, stations=stations)
        print "\tStation", key_st, stlo, stla, stel, count


        ######################################
        # test breakpoint resume
        tempfiles = glob.glob(prefix+"/*"+key_st+"*")
        if len(tempfiles)>=1:
            print "\t\t", key_st, "already extracted. Continue."
            continue

        # calculate theoretical traveltime
        try:
            dist, gcarc, az, baz, t_pf, name_pf, t_sf, name_sf = cal_travel_time(evla, evlo, evdp, stla, stlo)
            print "\t\t", dist, gcarc, az, baz, t_pf, name_pf, t_sf, name_sf
        except:
            print key_st, "failed calculation of traveltime. Continue."
            continue

        if t_pf<0:
            print "\t\tPf phase does not exist:", key, key_st
            continue

        # if t_sf<0:
        #     print "Sf phase does not exist:", key, key_st

        # if gcarc>8.0:
        #     continue

        # print dist, gcarc, az, baz, t_pf, name_pf, t_sf, name_sf


        # determine time window for extracting waveforms
        if not origin:  # using the theoretical traveltime of first arrival
            starttime = o + t_pf + tmin
            endtime = o + t_pf + tmax
        else:  # using origin
            starttime = o + tmin
            endtime = o + tmax

        starttime_bjt = starttime + 8.0 * 3600
        endtime_bjt = endtime + 8.0 * 3600

        t1 = o_bjt
        t2 = o_bjt + 24 * 3600

        t1 = UTCDateTime(str(t1)[:10])
        t2 = UTCDateTime(str(t2)[:10])

        yesterday = False
        tomorrow = False

        if starttime_bjt < t1:
            print "\t\tneed to merge data of yesterday. Today is", t1
            yesterday = True

        if endtime_bjt > t2:
            print "\t\tneed to merge data of tomorrow. Today is", t1
            tomorrow = True

        # print t1, t2

        # Beijing time, since the folder named after Beijing Time
        # folder_name = str(o_bjt)[:10].replace("-","")
        folder_name = o_bjt.strftime("%Y%m%d")
        # print folder_name

        if yesterday:
            o_bjt_yes = o - 24*3600
            folder_name_merge = o_bjt_yes.strftime("%Y%m%d")
            # print folder_name_merge

        if tomorrow:
            o_bjt_tmr = o + 24*3600
            folder_name_merge = o_bjt_tmr.strftime("%Y%m%d")
            # print folder_name_merge


        files = glob.glob(datapath+"/%s/*%s*%s" % (folder_name, key_st, suffix))
        # files_merge = glob.glob(datapath+"/%s/*%s*%s" % (folder_name_merge, key_st, suffix))
        # print files
        # print files_merge

        # os._exit(-1)

        # in this case, file not exist.
        if len(files)<1:
            print "\t\t", key_st, "no %s file found." % suffix
            continue


        # # temp = str(key)[:19]
        # temp = o.datetime
        # temp = temp.strftime("%Y%m%d%H%M%S")
        # # temp = temp.replace("-", "")
        # # temp = temp.replace(":", "")
        # # temp = temp.replace("T", "")
        # # print "temp", key, temp
        #
        # prefix = outpath+"/"+temp
        # # print prefix
        # try:
        #     os.makedirs(prefix)
        # except:
        #     pass

        #  read mseed, set header and write sac file
        for file in files:

            # bname = os.path.basename(file)
            # bname = bname.replace(suffix, ".sac")

            # print bname
            print "\t\t", file

            # the to-be-merged file might have different name with current one
            if yesterday or tomorrow:
                file2 = file.replace(folder_name, folder_name_merge)
                row = file2.split(".")
                row[-2]="*"
                file_merge = glob.glob(".".join(row))
                if len(file_merge)<1:
                    file_merge=""
                    print "\t\tno to-be-merged file found"
                else:
                    file_merge = file_merge[0]
            else:
                file_merge = ""
                # print file_merge


            # continue

            # st = read(file)
            # # print st
            #
            # st2 = st.copy()
            # try:
            #     st2.merge()
            # except:
            #     pass
            # try:
            #     st2 = st2.split()
            # except:
            #     pass
            # tr = st2[0]
            # tr.detrend(type="linear")
            # tr.trim(starttime, endtime, pad=True, fill_value=0.0)

            # # it is quite slow to read data directly from Internet, I test to copy data to local disk
            # os.system("cp %s ./" % file)
            # file_local = os.path.basename(file)
            # if file_merge!="":
            #     os.system("cp %s ./" % file_merge)
            #     file_merge_local = os.path.basename(file_merge)

            # tr = read_data(file=file_local, file_merge=file_merge_local)

            # in some case, the data is un-readable.
            try:
                tr = read_data(file=file, file_merge=file_merge)
            except:
                print "\t\t%s is not readable" % file
                continue
            # tr2 = read_data(filename=files_merge)

            tr.trim(starttime, endtime, pad=True, fill_value=0.0)

            temp = o.datetime
            temp = ".".join([temp.strftime("%Y.%j.%H.%M.%S"), "0000", tr.id, "M", "SAC"])


            oname = prefix+"/"+temp
            tr.write(oname, format="SAC")
            print "\t\t", oname

            tr = read(oname)[0]

            # tr = SACTrace.from_obspy_trace(trace=tr)
            #
            # tr.evdp = evdp
            # tr.evla = evla
            # tr.evlo = evlo
            # tr.mag = mag
            # tr.stla = stla
            # tr.stlo = stlo
            # tr.stel = stel
            # tr.az = az
            # tr.baz = baz
            # tr.dist = dist
            # tr.gcarc = gcarc
            #
            # if tr.kcmpnm[-1] == "E":
            #     tr.cmpaz = 90
            #     tr.cmpinc = 90
            # elif tr.kcmpnm[-1] == "N":
            #     tr.cmpaz = 0
            #     tr.cmpinc = 90
            # elif tr.kcmpnm[-1] == "Z":
            #     tr.cmpaz = 0
            #     tr.cmpinc = 0
            # else:
            #     pass
            #     # logger.warn("Not E|N|Z component")
            #
            # # set header
            # tr.t1 = t_pf
            # tr.t2 = t_sf
            # tr.kt1 = name_pf
            # tr.kt2 = name_sf
            #
            # tr.o = 0
            #
            # # change reference time
            # tr.nzyear = o.year
            # tr.nzjday = o.julday
            # tr.nzhour = o.hour
            # tr.nzmin = o.minute
            # tr.nzsec = o.second
            # tr.nzmsec = o.microsecond / 1000
            # tr.iztype = "io"    # Event origin time
            # tr.write(oname)

            tr.stats.sac.evdp = evdp
            tr.stats.sac.evla = evla
            tr.stats.sac.evlo = evlo
            tr.stats.sac.mag = mag
            tr.stats.sac.stla = stla
            tr.stats.sac.stlo = stlo
            tr.stats.sac.stel = stel
            tr.stats.sac.az = az
            tr.stats.sac.baz = baz
            tr.stats.sac.dist = dist
            tr.stats.sac.gcarc = gcarc

            if tr.stats.channel[-1] == "E":
                tr.stats.sac.cmpaz = 90
                tr.stats.sac.cmpinc = 90
            elif tr.stats.channel[-1] == "N":
                tr.stats.sac.cmpaz = 0
                tr.stats.sac.cmpinc = 90
            elif tr.stats.channel[-1] == "Z":
                tr.stats.sac.cmpaz = 0
                tr.stats.sac.cmpinc = 0
            else:
                pass
                # logger.warn("Not E|N|Z component")

            # set header
            tr.stats.sac.t1 = t_pf
            tr.stats.sac.t2 = t_sf
            tr.stats.sac.kt1 = name_pf
            tr.stats.sac.kt2 = name_sf

            tr.stats.sac.o = 0

            # change reference time
            tr.stats.sac.nzyear = o.year
            tr.stats.sac.nzjday = o.julday
            tr.stats.sac.nzhour = o.hour
            tr.stats.sac.nzmin = o.minute
            tr.stats.sac.nzsec = o.second
            tr.stats.sac.nzmsec = o.microsecond / 1000
            tr.stats.sac.iztype = 11    # Event origin time
            tr.write(oname, format="SAC")

        # os._exit(-1)



    pass




if __name__== "__main__":


    # eventfile = "2016.event2"
    # eventfile = "to_be_merged.event2"


    # raw miniseed datapath
    # mseedpath = "raw"

    # output sac datapath

    eventfile = "2016_mag6.event2"
    stationfile = "china.station2"
    mseedpath = "/data6/China_Data_2016"
    sacpath   = "data"




    # evts = read_event2(eventfile)
    # print evts

    # for e in evts:
    #     orig, evlo, evla, evdp, mag = read_event_line(e)
        # print orig, evlo, evla, evdp, mag

    # read event
    # e = evts[0]
    # orig, evlo, evla, evdp, mag = read_event_line(e)
    # print orig, evlo, evla, evdp, mag

    events = read_events(event_file=eventfile)
    stations = read_stations(station_file=stationfile)
    print "Number of events:", len(events)
    print "Number of stations:", len(stations)

    # for key in events.keys():
    #     print key, events[key]
    #     print {key:events[key]}

    # keys = events.keys()
    # get_event_waveform(key=keys[3], events=events, stations=stations, datapath=mseedpath, outpath=sacpath)

    keys = events.keys()
    keys.sort()
    for key in keys:
        get_event_waveform(key=key, events=events, stations=stations, datapath=mseedpath, outpath=sacpath,
                           origin=False, tmin=-100, tmax=1100)




    # for key in keys:
    #     get_event_waveform(key, events, stations=stations, datapath="raw/20160824", outpath="")





    # keys = events.keys()
    # key = keys[0]
    # evlo, evla, evdp, mag = event_info(key, events)
    # print key, evlo, evla, evdp, mag
    # print len(events)
    #
    # # os._exit(-1)
    #
    #
    # key = "AH.ANQ"
    # stlo, stla, stel = station_info(key, stations)
    # print key, stlo, stla, stel
    #
    # dist, gcarc, az, baz, t_pf, name_pf, t_sf, name_sf = cal_travel_time(evla, evlo, evdp, stla, stlo)
    # print dist, gcarc, az, baz, t_pf, name_pf, t_sf, name_sf



    # # test dict
    # es = test_dict(eventfile)
    # # print es
    #
    # for key in es.keys():
    #     # print key
    #     dict = es.get(key)
    #     print key, dict.values()


