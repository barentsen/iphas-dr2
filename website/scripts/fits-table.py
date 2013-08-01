#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Creates the table of IPHAS FITS catalogues for download."""
import httplib
import numpy as np


def filesize(host, url):
    conn = httplib.HTTPConnection(host)
    conn.request("HEAD", url)
    res = conn.getresponse()
    try:
        size = float(res.getheader('content-length')) / (1024.*1024.)
    except Exception, e:
        print 'Failed {0}'.format(url)
        print e
    if size < 950:
        result = '{0:.0f}&thinsp;MB'.format(size)
    else:
        result = '{0:.01f}&thinsp;GB'.format(size/1024.)
    print result
    return result   


if __name__ == '__main__':

    row = """
                <tr>
                  <td rowspan='2'>
                    {lon1}&deg;&thinsp;&le;&thinsp;l&thinsp;&lt;&thinsp;{lon2}&deg;</td>
                  <td>
                    b&thinsp;&lt;&thinsp;0&deg;</td>
                  <td>
                    <a href='http://{host}{light_url_a}'>iphas-dr2-{lon1:03d}a-light.fits</a> (<small>{light_size_a}</small>)
                  </td>
                  <td>
                    <a href='http://{host}{full_url_a}'>iphas-dr2-{lon1:03d}a.fits.gz</a> (<small>{full_size_a}</small>)
                  </td>
                </tr>
                <tr>
                  <td>
                    b&thinsp;&ge;&thinsp;0&deg;
                  </td>
                  <td>
                    <a href='http://{host}{light_url_b}'>iphas-dr2-{lon1:03d}b-light.fits</a> (<small>{light_size_b}</small>)
                  </td>
                  <td>
                    <a href='http://{host}{full_url_b}'>iphas-dr2-{lon1:03d}b.fits.gz</a> (<small>{full_size_b}</small>)
                  </td>
                </tr>
    """



    html = open('fits-table.html', 'w')
    for lon in np.arange(25, 220, 5):

        host = "stri-cluster.herts.ac.uk"
        root = "/~gb/iphas-dr2-catalogue"
        light_url_a = "{0}/light/iphas-dr2-{1:03d}a-light.fits".format(root, lon)
        light_size_a = filesize(host, light_url_a)
        light_url_b = "{0}/light/iphas-dr2-{1:03d}b-light.fits".format(root, lon)
        light_size_b = filesize(host, light_url_b)
        full_url_a = "{0}/full/iphas-dr2-{1:03d}a.fits.gz".format(root, lon)
        full_size_a = filesize(host, full_url_a)
        full_url_b = "{0}/full/iphas-dr2-{1:03d}b.fits.gz".format(root, lon)
        full_size_b = filesize(host, full_url_b)
        myrow = row.format(**{'lon1': lon, 'lon2': lon+5, 
                              'host':host, 
                              'light_url_a':light_url_a, 
                              'light_size_a':light_size_a, 
                              'light_url_b':light_url_b, 
                              'light_size_b':light_size_b,
                              'full_url_a': full_url_a, 
                              'full_size_a': full_size_a,
                              'full_url_b': full_url_b, 
                              'full_size_b': full_size_b})
        html.write(myrow)
    html.close()

