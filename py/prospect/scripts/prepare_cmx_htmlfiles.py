"""
prospect.scripts.prepare_cmx_htmlfiles
===================================

Derived from prepare_htmlfiles
"""


import os, glob, stat
import argparse
from desiutil.log import get_logger

from jinja2 import Environment, FileSystemLoader

def parse() :

    parser = argparse.ArgumentParser(description="Write html index pages")
    parser.add_argument('--webdir', help='Base directory for webpages', type=str)
    parser.add_argument('--template_dir', help='Template directory', type=str, default=None)
    args = parser.parse_args()
    return args

def main(args) :

    log = get_logger()

    webdir = args.webdir
    template_dir = args.template_dir
    if template_dir is None : template_dir="../templates" # TMP, to define better.

    env = Environment(loader=FileSystemLoader(template_dir))
    template_index = env.get_template('template_index.html')
    template_expolist = env.get_template('template_cmx_expo_list.html')

    exposures = os.listdir( os.path.join(webdir,"exposures") )
    for expo in exposures :
        expo_dir = os.path.join(webdir,"exposures",expo)
        
        available_subsets = []
        for specnum in ["0","1","2","3","4","5","6","7","8","9"] :
            subsets = []
            for i_sub in range(1,11) :
                if os.path.exists(expo_dir+"/specviewer_"+expo+"_spectro"+specnum+"_"+str(i_sub)+".html") :
                    subsets.append(str(i_sub))
            available_subsets.append(subsets)
        pagetext = template_expolist.render(expo=expo, subsets_0=available_subsets[0], subsets_1=available_subsets[1], subsets_2=available_subsets[2], subsets_3=available_subsets[3], subsets_4=available_subsets[4], subsets_5=available_subsets[5], subsets_6=available_subsets[6], subsets_7=available_subsets[7], subsets_8=available_subsets[8], subsets_9=available_subsets[9])

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

