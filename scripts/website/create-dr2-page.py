"""Produces the DR2 page for the IPHAS website."""
import os
import json
import datetime
import httplib
from astropy import log
from jinja2 import Environment, FileSystemLoader

from dr2 import constants

# Load the column definitions from the JSON file
filename_columns = os.path.join(constants.LIBDIR, 'columns.json')
columns = json.loads(open(filename_columns).read())


def get_filesize(host, url):
    """Returns the filesize of a remote document."""
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
    log.info('Filesize of {0} = {1}'.format(url, result))
    return result


args = {'host': 'www.iphas.org',
        'dir_light': '/data/dr2/light',
        'dir_full': '/data/dr2/full',
        'last_update': datetime.datetime.utcnow().strftime("%Y-%m-%d"),
        'columns': columns}

if __name__ == '__main__':
    env = Environment(loader=FileSystemLoader('.'),
                      trim_blocks=True,
                      lstrip_blocks=True)
    env.globals['get_filesize'] = get_filesize
    template = env.get_template('dr2-template.html')

    with open("dr2.shtml", "w") as f:
        f.write(template.render(args))
