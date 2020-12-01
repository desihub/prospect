# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
======================================
prospect.scripts.prepare_cmx_htmlfiles
======================================

Write html index pages from existing static pages/images produced by other scripts.
"""


import os, glob, stat
import argparse
from desiutil.log import get_logger

from jinja2 import Environment, FileSystemLoader

def _parse():

    parser = argparse.ArgumentParser(description="Write html index pages for CMX exposures")
    parser.add_argument('--webdir', help='Base directory for webpages', type=str)
    parser.add_argument('--template_dir', help='Template directory', type=str, default=None)
    parser.add_argument('--nspecperfile', help='Number of spectra in each prospect html page', type=int, default=50)
    args = parser.parse_args()
    return args

def main():

    args = _parse()
    log = get_logger()

    webdir = args.webdir
    template_dir = args.template_dir
    if template_dir is None :
        template_dir = os.path.join(os.path.dirname(__file__),os.pardir,os.pardir,os.pardir,"templates")

    env = Environment(loader=FileSystemLoader(template_dir))
    template_index = env.get_template('template_index.html')
    template_expolist = env.get_template('template_cmx_expo_list.html')

    exposures = os.listdir( os.path.join(webdir,"exposures") )
    for expo in exposures :
        expo_dir = os.path.join(webdir,"exposures",expo)

        available_subsets = {}
        for specnum in ["0","1","2","3","4","5","6","7","8","9"] :
            subsets = []
            i_sub = 1
            while os.path.exists(expo_dir+"/specviewer_"+expo+"_spectro"+specnum+"_"+str(i_sub)+".html") :
                subsets.append(str(i_sub))
                i_sub += 1
            available_subsets[specnum] = subsets
        pagetext = template_expolist.render(expo=expo, subset_dict=available_subsets, nspecperfile=args.nspecperfile)

        with open( os.path.join(expo_dir,"index_"+expo+".html"), "w") as fh:
            fh.write(pagetext)
            fh.close()
        for x in glob.glob(expo_dir+"/*.html") :
            st = os.stat(x)
            os.chmod(x, st.st_mode | stat.S_IROTH) # "chmod a+r "
        st = os.stat(expo_dir)
        os.chmod(expo_dir, st.st_mode | stat.S_IROTH | stat.S_IXOTH) # "chmod a+rx "

        log.info("Subdirectory done : "+expo)

    pagetext = template_index.render(exposures=exposures)
    indexfile = os.path.join(webdir,"index.html")
    with open(indexfile, "w") as fh:
        fh.write(pagetext)
        fh.close()
        st = os.stat(indexfile)
        os.chmod(indexfile, st.st_mode | stat.S_IROTH) # "chmod a+r"
        log.info("Main index done")
