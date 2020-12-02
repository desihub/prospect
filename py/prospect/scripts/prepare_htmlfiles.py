# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
==================================
prospect.scripts.prepare_htmlfiles
==================================

Write html index pages from existing pages/images produced by :mod:`~prospect.plotframes`.
"""


import os, glob, stat
import argparse
from pkg_resources import resource_filename
from desiutil.log import get_logger

from jinja2 import Environment, FileSystemLoader


def _parse():

    parser = argparse.ArgumentParser(description="Write html index pages")
    parser.add_argument('--webdir', help='Base directory for webpages', type=str, default=None)
    parser.add_argument('--template_dir', help='Template directory', type=str, default=resource_filename('prospect', 'templates'))
    parser.add_argument('--pixels', help='Pixel-based directory structure', action='store_true')
    parser.add_argument('--targets', help='Target-based directory structure', action='store_true')
    parser.add_argument('--exposures', help='Exposure-based directory structure', action='store_true')
    parser.add_argument('--with_thumbs', help='Also index thumbnail images', action='store_true')
    args = parser.parse_args()
    return args


def prepare_subdir(subdir, entry, template_index, template_vignette, target=None, do_expo=False, info=None, with_thumbs=False) :

    # TODO Needs restructure !!!

    pattern = entry
    if target is not None :
        pattern = target+"*"+entry

    spec_pages = glob.glob( subdir+"/specviewer_"+pattern+"_*.html" )
#     subsets = [ x[len(subdir+"/specviewer_"+pattern)+1:-5] for x in spec_pages ]
#     subsets.sort(key=int)
    subsets = [str(x+1) for x in range(len(spec_pages))]
    if with_thumbs : img_list = glob.glob( subdir+"/vignettes/*.png" )
    pp = [x for x in spec_pages if "_1.html" in x]
    pp=pp[0]
    basename = pp[len(subdir)+1:-7]
    if with_thumbs : nspec = len(img_list)
    else : nspec = 50 # TODO Ne pas laisser ca...
    if target is None and do_expo==False :
        pagetext = template_index.render(pixel=entry, subsets=subsets, nspec=nspec)
    elif target is None and do_expo==True :
        pagetext = template_index.render(expo=entry, subsets=subsets, nspec=nspec)
    else :
        pagetext = template_index.render(pixel=entry, subsets=subsets, nspec=nspec, target=target, info=info, basename=basename)

    with open( os.path.join(subdir,"index_"+entry+".html"), "w") as fh:
        fh.write(pagetext)
        fh.close()
    if with_thumbs :
        for subset in subsets :
            img_sublist = [ os.path.basename(x) for x in img_list if entry+"_"+subset in x ]
            pagetext = template_vignette.render(set=entry, i_subset=subset, n_subsets=len(subsets), imglist=img_sublist)
            with open( os.path.join(subdir,"vignettelist_"+entry+"_"+subset+".html"), "w") as fh:
                fh.write(pagetext)
                fh.close()
    for x in glob.glob(subdir+"/*.html") :
        st = os.stat(x)
        os.chmod(x, st.st_mode | stat.S_IROTH) # "chmod a+r "
    st = os.stat(subdir)
    os.chmod(subdir, st.st_mode | stat.S_IROTH | stat.S_IXOTH) # "chmod a+rx "
    if with_thumbs :
        thedir = os.path.join(subdir,"vignettes")
        st = os.stat(thedir)
        os.chmod(thedir, st.st_mode | stat.S_IROTH | stat.S_IXOTH) # "chmod a+rx "
        for x in glob.glob(subdir+"/vignettes/*.png") :
            st = os.stat(x)
            os.chmod(x, st.st_mode | stat.S_IROTH) # "chmod a+r "


def main():
    args = _parse()
    log = get_logger()

    webdir = args.webdir
    if webdir is None : webdir = os.environ["DESI_WWW"]+"/users/armengau/svdc2019c" # TMP, for test
    template_dir = args.template_dir
    # if template_dir is None : template_dir="../templates" # TMP, to define better.

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
            prepare_subdir(expo_dir, expo, template_expolist, template_vignettelist, do_expo=True, with_thumbs=args.with_thumbs)
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
#        target_dict = { # TODO locate elsewhere - no hardcode  ## for v1
#                 "BGS_ANY" : "bgs_targets",
#                 "ELG" : "elg_targets",
#                 "LRG" : "lrg_targets",
#                 "QSO" : "qso_targets",
#                 "MWS_ANY" : "mws_targets"}
        target_list = [ ["MWS_ANY", "mws", "All MWS targets"], ## for v2, still tmp
                        ["BGS_ANY", "bgs_bluesquare", "Blue square BGS : r>19.5"],
                        ["BGS_ANY", "bgs_greencircle", "Green circle BGS : r<19.5"],
                        ["LRG", "lrg", "All LRG targets"],
                        ["ELG", "elg_bluesquare", "Blue square ELGs : DeltaChi2 in [40 - 100]"],
                        ["ELG", "elg_greencircle", "Green circle ELGs : DeltaChi2>100"],
                        ["ELG", "elg_blackdiamond", "Black diamond ELGs : DeltaChi2<40"],
                        ["QSO", "qso_bluesquare", "Blue square QSOs : g>22.5"],
                        ["QSO", "qso_greencircle", "Green circle QSOs : g<22.5"] ]
        for i in range(len(target_list)) :
            target_dir = target_list[i][1]
            target_cat = target_list[i][0]
            target_pixels[target_dir] = os.listdir( os.path.join(webdir,target_dir) )
            for pix in target_pixels[target_dir] :
                pixel_dir = os.path.join(webdir,target_dir,pix)

                prepare_subdir(pixel_dir, pix, template_targetlist, template_vignettelist, target=target_cat, info=target_list[i][2] )
                log.info("Subdirectory done : "+pix)
        target_pixels['ELG']=None
        target_pixels['QSO']=None
        target_pixels['BGS_ANY']=None
    else : target_pixels={'BGS_ANY':None,'ELG':None,'LRG':None,'QSO':None,'MWS_ANY':None}

    # Main index # TODO improve template handling target-based pages
    pagetext = template_index.render(pixels=pixels, exposures=exposures, bgs_pixels=target_pixels['BGS_ANY'],
            elg_pixels=target_pixels['ELG'], lrg_pixels=target_pixels['lrg'], qso_pixels=target_pixels['QSO'], qsob_pix=target_pixels['qso_bluesquare'], qsog_pix=target_pixels['qso_greencircle'], elgg_pix=target_pixels['elg_greencircle'], elgb_pix=target_pixels['elg_bluesquare'], elgbb_pix=target_pixels['elg_blackdiamond'],bgsb_pix=target_pixels['bgs_bluesquare'], bgsg_pix=target_pixels['bgs_greencircle'],
            mws_pixels=target_pixels['mws'])
    indexfile = os.path.join(webdir,"index.html")
    with open(indexfile, "w") as fh:
        fh.write(pagetext)
        fh.close()
        st = os.stat(indexfile)
        os.chmod(indexfile, st.st_mode | stat.S_IROTH) # "chmod a+r"
        log.info("Main index done")
