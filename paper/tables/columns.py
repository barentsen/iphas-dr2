"""Creates a TeX table explaining the IPHAS DR2 catalogue columns."""
import os
import json
import re

from dr2 import constants

def html2tex(text):
    """Returns TeX-formatted strings."""
    text = text.replace('&alpha;', r'$\alpha$')
    text = text.replace('&gt;', r'$>$')
    text = text.replace('&lt;', r'$<$')
    text = text.replace('&ge;', r'$\geq$')
    text = text.replace('&le;', r'$\leq$')
    text = text.replace('&amp;', r'\&')
    text = text.replace('&thinsp;', r'\,')
    text = text.replace('#', r'$\#$')
    text = text.replace('_', r'\_')
    text = re.sub('<[^<]+?>', '', text) # Remove HTML
    return text

# Load the column definitions from the JSON file in the DR2 package
filename = os.path.join(constants.LIBDIR, 'columns.json')
data = json.loads(open(filename).read())

# Write TeX table
out = open('columns.tex', 'w')
for i, column in enumerate(data):
    column['n'] = i + 1
    column['desc'] = html2tex(column['desc'])
    out.write("{n} & {name} & {type} & {unit} & {desc} \\\\\n".format(**column))
out.close()

