"""
prospect.scripts.prepare_htmlfiles
===================================

Write html index pages from existing pages/images produced by plotframes
"""


import os, glob, stat
import argparse
from desiutil.log import get_logger

from jinja2 import Environment, FileSystemLoader


def parse() :

    parser = argparse.ArgumentParser(description="Write html index pages")
    parser.add_argument('--webdir', help='Base directory for webpages', type=str, default=None)
    parser.add_argument('--template_dir', help='Template directory', type=str, default=None)
    parser.add_argument('--pixels', help='Pixel-based directory structure', action='store_true')
    parser.add_argument('--targets', help='Target-based directory structure', action='store_true')
    parser.add_argument('--exposures', help='Exposure-based directory structure', action='store_true')
    args = parser.parse_args()
    return args


def prepare_subdir(subdir, entry, template_index, template_vignette, target=None, do_expo=False) :

    # TODO Needs restructure !!!
    
    pattern = entry
    if target is not None :
        pattern = target+"_"+entry
    
    spec_pages = glob.glob( subdir+"/specviewer_"+pattern+"_*.html" )
    subsets = [ x[len(subdir+"/specviewer_"+pattern)+1:-5] for x in spec_pages ]
    subsets.sort(key=int)
    img_list = glob.glob( subdir+"/vignettes/*.png" )
    nspec = len(img_list)
    if target is None and do_expo==False :
        pagetext = template_index.render(pixel=entry, subsets=subsets, nspec=nspec)
    elif target is None and do_expo==True :
        pagetext = template_index.render(expo=entry, subsets=subsets, nspec=nspec)
    else :
        pagetext = template_index.render(pixel=entry, subsets=subsets, nspec=nspec, target=target)

    with open( os.path.join(subdir,"index_"+entry+".html"), "w") as fh:
        fh.write(pagetext)
        fh.close()
    for subset in subsets :
        img_sublist = [ os.path.basename(x) for x in img_list if entry+"_"+subset in x ]
        pagetext = template_vignette.render(set=entry, i_subset=subset, n_subsets=len(subsets), imglist=img_sublist)
        with open( os.path.join(subdir,"vignettelist_"+entry+"_"+subset+".html"), "w") as fh:
            fh.write(pagetext)
            fh.close()
    for x in glob.glob(subdir+"/*.html") :
        st = os.stat(x)
        os.chmod(x, st.st_mode | stat.S_IROTH) # "chmod a+r "
    for thedir in [subdir, os.path.join(subdir,"vignettes") ] :
        st = os.stat(thedir)
        os.chmod(thedir, st.st_mode | stat.S_IROTH | stat.S_IXOTH) # "chmod a+rx "
    for x in glob.glob(subdir+"/vignettes/*.png") :
        st = os.stat(x)
        os.chmod(x, st.st_mode | stat.S_IROTH) # "chmod a+r "


def main(args) :

    log = get_logger()

    webdir = args.webdir
    if webdir is None : webdir = os.environ["DESI_WWW"]+"/users/armengau/svdc2019c" # TMP, for test
    template_dir = args.template_dir
    if template_dir is None : template_dir="../templates" # TMP, to define better.

    env = Environment(loader=FileSystemLoader(template_dir))
    template_index = env.get_template('template_index.html')
    template_expolist = env.get_template('template_expo_list.html')
    template_pixellist = env.get_template('template_pixel_list.html')
    template_targetlist = env.get_template('template_target_list.html')
    template_vignettelist = env.get_template('template_vignettelist.html')

    if args.exposures :
        exposures = os.listdir( os.path.join(webdir,"exposures") )
        for expo in exposures :
            expo_dir = os.path.join(webdir,"exposures",expo)
            prepare_subdir(expo_dir, expo, template_expolist, template_vignettelist, do_expo=True)
            log.info("Subdirectory done : "+expo)
    else : exposures = None

    if args.pixels :
        pixels = os.listdir( os.path.join(webdir,"pixels") )
        for pix in pixels :
            pixel_dir = os.path.join(webdir,"pixels",pix)
            prepare_subdir(pixel_dir, pix, template_pixellist, template_vignettelist)
            log.info("Subdirectory done : "+pix)
    else : pixels = None

    if args.targets :
        target_pixels = dict()
        target_dict = { # TODO locate elsewhere - no hardcode
                "BGS_ANY" : "bgs_targets",
                "ELG" : "elg_targets",
                "LRG" : "lrg_targets",
                "QSO" : "qso_targets",
                "MWS_ANY" : "mws_targets"}
        for target_cat, target_dir in target_dict.items() :
            target_pixels[target_cat] = os.listdir( os.path.join(webdir,target_dir) )
            for pix in target_pixels[target_cat] :
                pixel_dir = os.path.join(webdir,target_dir,pix)
                
                prepare_subdir(pixel_dir, pix, template_targetlist, template_vignettelist, target=target_cat )
                log.info("Subdirectory done : "+pix)
    else : target_pixels={'BGS_ANY':None,'ELG':None,'LRG':None,'QSO':None,'MWS_ANY':None}

    # Main index # TODO improve template handling target-based pages
    pagetext = template_index.render(pixels=pixels, exposures=exposures, bgs_pixels=target_pixels['BGS_ANY'], 
            elg_pixels=target_pixels['ELG'], lrg_pixels=target_pixels['LRG'], qso_pixels=target_pixels['QSO'],
            mws_pixels=target_pixels['MWS_ANY']) 
    indexfile = os.path.join(webdir,"index.html")
    with open(indexfile, "w") as fh:
        fh.write(pagetext)
        fh.close()
        st = os.stat(indexfile)
        os.chmod(indexfile, st.st_mode | stat.S_IROTH) # "chmod a+r"
        log.info("Main index done")
