# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 14:56:51 2016

@author: jdorvinen
"""

import html5lib
from bs4 import BeautifulSoup
import requests
transects = ['F'+str(i+1) for i in range(84)]
data_sources = "http://dune.seagrant.sunysb.edu/nyshore/profile_viewer/index.jsp?&currentProfile={}"
data_source = data_sources.format(transects[0])
# <codecell>
r = requests.get(data_source)

# <codecell>
r.connection.close()