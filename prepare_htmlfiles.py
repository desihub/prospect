# Script to have html files/index from templates + existing images/pages created by specviewer

import os, glob
from jinja2 import Environment, FileSystemLoader


webdir=os.environ['HOME']+"/prospect/website"

### Exposure-based : for each expo create an index of fiber subsets, and a vignette webpage for each fiber subset 
## Arborescence : webdir/exposures/expoN/ =>  specviewer_expoN_fibersetM.html ; vignettes/*.png

env = Environment(loader=FileSystemLoader('templates')) # path to templates
template_index = env.get_template('template_index.html')
template_expolist = env.get_template('template_expo_list.html')
template_pixellist = env.get_template('template_pixel_list.html')
template_vignettelist = env.get_template('template_vignettelist.html')

exposures=os.listdir(webdir+"/exposures")
for expo in exposures :
    basedir = webdir+"/exposures/"+expo
    pp=glob.glob(basedir+"/specviewer_"+expo+"_*.html")
    subsets = [int(x[x.find("fiberset")+8:-5]) for x in pp]
    subsets.sort()
    subsets = [str(x) for x in subsets]
    pp=glob.glob(basedir+"/vignettes/*.png")
    nspec = len(pp)
    pagetext=template_expolist.render(expo=expo, fiberlist=subsets, nspec=nspec)
    with open(basedir+"/index_"+expo+".html", "w") as fh:
        fh.write(pagetext)
        fh.close()
    for fiber in subsets :
        pp = glob.glob(basedir+"/vignettes/"+expo+"_fiberset"+fiber+"_*.png")
        vignettelist = [os.path.basename(x) for x in pp]
        pagetext = template_vignettelist.render(set=expo, i_subset=fiber, n_subsets=len(subsets), imglist=vignettelist)
        with open(basedir+"/vignettelist_"+expo+"_"+fiber+".html", "w") as fh:
            fh.write(pagetext)
            fh.close()

pixels = os.listdir(webdir+"/pixels")
for pix in pixels :
    basedir = webdir+"/pixels/"+pix
    pp=glob.glob(basedir+"/specviewer_"+pix+"_*.html")
    subsets = [int(x[x.rfind("_")+1:-5]) for x in pp]
    subsets.sort()
    subsets = [str(x) for x in subsets]
    pp=glob.glob(basedir+"/vignettes/*.png")
    nspec = len(pp)
    pagetext=template_pixellist.render(pixel=pix, subsets=subsets, nspec=nspec)
    with open(basedir+"/index_"+pix+".html", "w") as fh:
        fh.write(pagetext)
        fh.close()
    for subset in subsets :
        pp = glob.glob(basedir+"/vignettes/"+pix+"_"+subset+"_*.png")
        vignettelist = [os.path.basename(x) for x in pp]
        pagetext = template_vignettelist.render(set=pix, i_subset=subset, n_subsets=len(subsets), imglist=vignettelist)
        with open(basedir+"/vignettelist_"+pix+"_"+subset+".html", "w") as fh:
            fh.write(pagetext)
            fh.close()

pagetext = template_index.render(pixels=pixels, exposures=exposures)
with open(webdir+"/index.html", "w") as fh:
    fh.write(pagetext)
    fh.close()
