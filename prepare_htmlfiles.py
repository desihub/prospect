# Script to have html files/index from templates + existing images/pages created by specviewer

import os, glob
from jinja2 import Environment, FileSystemLoader


webdir=os.environ['HOME']+"/prospect/website"

### Exposure-based : for each expo create an index of fiber subsets, and a vignette webpage for each fiber subset 
## Arborescence : webdir/exposures/expoN/ =>  specviewer_expoN_fibersetM.html ; vignettes/*.png

env = Environment(loader=FileSystemLoader('templates')) # path to templates
template_expolist = env.get_template('template_expo_list.html')
template_vignettelist = env.get_template('template_vignettelist.html')

exposures=os.listdir(webdir+"/exposures")
for expo in exposures :
    pp=glob.glob(webdir+"/exposures/"+expo+"/specviewer_"+expo+"_fiberset*.html")
    fiberlist = [int(x[x.find("fiberset")+8:-5]) for x in pp]
    fiberlist.sort()
    fiberlist = [str(x) for x in fiberlist]
    pagetext=template_expolist.render(expo=expo, fiberlist=fiberlist)
    with open(webdir+"/exposures/"+expo+"/index_"+expo+".html", "w") as fh:
        fh.write(pagetext)
        fh.close()
    for fiber in fiberlist :
        pp = glob.glob(webdir+"/exposures/"+expo+"/vignettes/"+expo+"_fiberset"+fiber+"_*.png")
        vignettelist = [os.path.basename(x) for x in pp]
        pagetext = template_vignettelist.render(expo=expo, i_subset=fiber, n_subsets=len(fiberlist), imglist=vignettelist)
        with open(webdir+"/exposures/"+expo+"/vignettelist_"+expo+"_fiberset"+fiber+".html", "w") as fh:
            fh.write(pagetext)
            fh.close()

## TODO similar for pixels (also maybe try a clickable skymap in bokeh)
    