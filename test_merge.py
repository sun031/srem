#!/usr/bin/env python


import oat2

cedc_file = "data/oat_cedc/20120614030808.oat2"

oat2.merge_isc_cedc(cedc_file, isc_dir="data/oat_isc", origin_error=60, location_error=0.5)
