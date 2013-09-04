#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Find potential anchor nights."""
import numpy as np
from dr2 import constants

qc = constants.IPHASQC


BLACKLIST = [20040806, 20050920, 20051106, 20061110, 20061127]


def parse_extinction():
    result = []
    for value in qc['ext_r_carlsberg']:
        if value == "":
            result.append(0)
        else:
            result.append(float(value[0:5]))
    return np.array(result)

extinction = parse_extinction()


def consider_night(night):
    if night in BLACKLIST:
        return False

    mask_night = (qc['night'] == night)

    # Night must have low humidity
    # Carlsberg must agree 
    # ext_r_carlsberg < 0.2
    # hours_nonphot_carlsberg < 0.1
    # hour_phot_carlsberg > 5
    # avg_seeing_r < 1.5
    # lost_weather == "00:00"
    # (qc['ext_r_carlsberg'] < 0.2) &
    mask_carlsberg = ((qc['hours_nonphot_carlsberg'] < 0.3) &
                      (qc['hours_phot_carlsberg'] > 4))
    mask_weather = ((qc['lost_weather'] == "00:00") &
                    (qc['hum_avg'] < 70))
    #mask_moon = (qc['moon_phase'] < 100)
    mask_seeing = np.mean(qc['seeing_r'][mask_night]) < 1.7

    apass_r, apass_i = get_avg_apass_shifts(mynight)
    sdss_r, sdss_i = get_avg_sdss_shifts(mynight)
    mask_surveys = ((np.abs(apass_r) < 0.05) &
                    (np.abs(apass_i) < 0.05) &
                    (np.abs(sdss_r) < 0.05) &
                    (np.abs(sdss_i) < 0.05))
    mask_hasanchors = (qc['anchor'] == 1)

    mask_anchor = (mask_carlsberg & mask_weather & 
                   mask_seeing & mask_surveys & mask_hasanchors)
    if np.any(mask_anchor[mask_night]):
        return True
    else:
        return False

def get_weather_comments(night):
    mask_night = (qc['night'] == night)
    if mask_night.sum() > 0:
        return qc['comments_weather'][mask_night][0]
    return ''

def get_avg_apass_shifts(night):
    mask_night = (qc['night'] == night)
    mask_apass = ((qc['rmatch_apassdr7'] > 30) &
                  (qc['imatch_apassdr7'] > 30))
    mask_use = mask_night & mask_apass
    if mask_use.sum() > 0:
        return (np.mean(qc['rshift_apassdr7'][mask_use]),
                np.mean(qc['ishift_apassdr7'][mask_use]))
    return (0.00, 0.00)

def get_avg_sdss_shifts(night):
    mask_night = (qc['night'] == night)
    mask_apass = ((qc['rmatch_sdss'] > 30) &
                  (qc['imatch_sdss'] > 30))
    mask_use = mask_night & mask_apass
    if mask_use.sum() > 0:
        return (np.mean(qc['rshift_sdss'][mask_use]),
                np.mean(qc['ishift_sdss'][mask_use]))
    return (0.00, 0.00)

if __name__ == '__main__':
    anchornights = []

    n_fields, n_anchors = 0, 0
    for mynight in np.unique(qc['night']):
        is_anchor = consider_night(mynight)    
        if is_anchor:
            anchornights.append(mynight)
            n_fields_night = (qc['is_dr2'] & (qc['night'] == mynight)).sum()
            n_anchors_night = (qc['is_dr2'] & (qc['night'] == mynight) & (qc['anchor'] == 1)).sum()
            apass_r, apass_i = get_avg_apass_shifts(mynight)
            sdss_r, sdss_i = get_avg_sdss_shifts(mynight)
            print "{0} ({1:03d}/{2:03d}) [{3:+.2f} {4:+.2f}] [{5:+.2f} {6:+.2f}] {7}".format(
                    mynight,
                    n_anchors_night,
                    n_fields_night,
                    apass_r, apass_i,
                    sdss_r, sdss_i,
                    get_weather_comments(mynight))

            n_fields += n_fields_night
            n_anchors += n_anchors_night

    print "Total: {0} fields ({1} already anchor)".format(n_fields, n_anchors)

    for mynight in anchornights:
        print mynight
