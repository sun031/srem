#!/usr/bin/env python


import gmt

st2file = "data/oat_cedc_crude.oat2"
psfile  = "figs/oat_cedc_crude.ps"

st2file = "data/oat_isc_crude.oat2"
psfile  = "figs/oat_isc_crude_station.ps"



# plot the epicenter
gmt5 = gmt.Gmt()

Jh = "L105/36/17/55/10c"
Rh = "70/140/17/55"

gmt5.set("MAP_FRAME_WIDTH", "0.3p")
gmt5.comment("etopo")
gmt5.shell("grdgradient etopo_cn.grd -Gjunk.grd -A0/270 -Nt1.0 -V")
gmt5.cmd("makecpt", "-Cglobe -T-10000/10000/100 -Z > colors.cpt")
gmt5.cmd("grdimage", "etopo_cn.grd -Ijunk.grd -J%s -R%s -K -Ccolors.cpt > %s" % (Jh, Rh, psfile))

gmt5.comment("coast and tectonic lines")
gmt5.cmd("pscoast", "-J%s -R%s -W1/0.5p -N1/0.5p -N2 -K -O -BWSNE -Bxa10 -Bya10 >> %s" % (Jh, Rh, psfile) )
gmt5.cmd("psxy", "China_tectonic.dat -J%s -R%s  -K -O -W0.5p,blue,- >> %s" % (Jh, Rh, psfile))

# plot permanent stations
gmt5.shell("cat %s | awk '{print $11,$10}' > temp.xy" % st2file)
gmt5.shell("sort -u temp.xy > temp2.xy")
gmt5.cmd("psxy", "temp2.xy -J%s -R%s  -K -O -St0.1c -W0.35p,red  -Gred >> %s" % (Jh, Rh, psfile))
# gmt5.cmd("psxy", "events_in_china.xy -J%s -R%s  -K -O -Sc0.15c -Gred -W0.35p,red  >> %s" % (Jh, Rh, psfile))

gmt5.comment("end")
gmt5.cmd("psxy", "-J -R -O -T >> %s" % psfile)
gmt5.cmd("psconvert", "-A -P -Tj %s" % psfile)
gmt5.cmd("psconvert", "-A -P -Tf %s" % psfile)
gmt5.execute()

