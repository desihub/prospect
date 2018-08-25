import numpy as np
import bokeh.plotting as bk

from bokeh.models import ColumnDataSource, CustomJS
from bokeh.models.widgets import Slider, Button
from bokeh.layouts import widgetbox
# from bokeh.layouts import row, column

import desispec.io
# bp.output_notebook()

framefiles = [
    'data/cframe-b0-00000020.fits',
    'data/cframe-r0-00000020.fits',
    'data/cframe-z0-00000020.fits',    
]
frames = list()
cds_spectra = list()

for filename in framefiles:
    fr = desispec.io.read_frame(filename)
    bad = (fr.ivar == 0.0) | (fr.mask != 0)
    fr.flux[bad] = np.nan
    frames.append(fr)

    cdsdata=dict(
        origwave=fr.wave.copy(),
        plotwave=fr.wave.copy(),
        plotflux=fr.flux[0],
        )

    for i in range(fr.nspec):
        key = 'origflux'+str(i)
        cdsdata[key] = fr.flux[i]
    
    cds_spectra.append(
        bk.ColumnDataSource(cdsdata, name=fr.meta['CAMERA'][0])
        )

# fig = bk.figure(height=450, width=800)
fig = bk.figure(height=350, width=800)
colors = dict(b='#1f77b4', r='#d62728', z='maroon')
for spec in cds_spectra:
    fig.line('plotwave', 'plotflux', source=spec, line_color=colors[spec.name])

# zslider_callback  = CustomJS(args=dict(spectrum=spectrum), code="""
#     var data = spectrum.data
#     var z = cb_obj.value
#     var origwave = data['origwave']
#     var plotwave = data['plotwave']
#     for (var i=0; i<wave.length; i++) {
#         wave[i] = obswave[i] / (1+z)
#     }
#     spectrum.change.emit()
#     """)
#
# zslider = Slider(start=0, end=2.0, value=0.0, step=0.01, title='Redshift')
# zslider.js_on_change('value', zslider_callback)

ifiberslider = Slider(start=0, end=fr.nspec-1, value=0, step=1, title='Fiber')
smootherslider = Slider(start=1, end=21, value=1, step=1, title='Smooth')

#-----
args = dict(
    spectra=cds_spectra,
    ifiberslider=ifiberslider,
    smootherslider=smootherslider
    )
update_plot = CustomJS(args=args, code="""
    ifiber = ifiberslider.value
    nsmooth = smootherslider.value
    for (var i=0; i<spectra.length; i++) {
        var data = spectra[i].data
        var plotflux = data['plotflux']
        var origflux = data['origflux'+ifiber]
        for (var j=0; j<plotflux.length; j++) {
            plotflux[j] = 0.0
            var n = 0
            for (var k=Math.max(0, j-nsmooth); k<Math.min(plotflux.length, j+nsmooth); k++) {
                fx = origflux[k]
                if(fx == fx) {
                    plotflux[j] = plotflux[j] + fx
                    n++
                }
            }
            plotflux[j] = plotflux[j] / n
        }
        spectra[i].change.emit()
    }
""")
smootherslider.js_on_change('value', update_plot)
ifiberslider.js_on_change('value', update_plot)

### bk.show(bk.Column(fig, widgetbox(zslider), widgetbox(ifiberslider)))
bk.show(bk.Column(fig, widgetbox(ifiberslider, width=800), widgetbox(smootherslider)))

#--- DEBUG ---
# import IPython
# IPython.embed()
#--- DEBUG ---