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
    frames.append(fr)

    cdsdata=dict(
        obswave=fr.wave.copy(),
        wave=fr.wave.copy(),
        flux=fr.flux[0],
        )

    for i in range(fr.nspec):
        key = 'flux'+str(i)
        cdsdata[key] = fr.flux[i]
    
    cds_spectra.append(
        bk.ColumnDataSource(cdsdata, name=fr.meta['CAMERA'][0])
        )

fig = bk.figure(height=450, width=800)
colors = dict(b='#1f77b4', r='#d62728', z='maroon')
for spec in cds_spectra:
    fig.line('wave', 'flux', source=spec, line_color=colors[spec.name])

# zslider_callback  = CustomJS(args=dict(spectrum=spectrum), code="""
#     var data = spectrum.data
#     var z = cb_obj.value
#     var wave = data['wave']
#     var obswave = data['obswave']
#     for (var i=0; i<wave.length; i++) {
#         wave[i] = obswave[i] / (1+z)
#     }
#     spectrum.change.emit()
#     """)
#
# zslider = Slider(start=0, end=2.0, value=0.0, step=0.01, title='Redshift')
# zslider.js_on_change('value', zslider_callback)

#-----
ifiberslider_callback  = CustomJS(args=dict(spectra=cds_spectra), code="""
    var ifiber = cb_obj.value
    for (var i=0; i<spectra.length; i++) {
        var data = spectra[i].data
        var flux = data['flux']
        var newflux = data['flux'+ifiber]
        for (var j=0; j<flux.length; j++) {
            flux[j] = newflux[j]
        }
        spectra[i].change.emit()
    }
    """)

ifiberslider = Slider(start=0, end=fr.nspec-1, value=0, step=1, title='Fiber')
ifiberslider.js_on_change('value', ifiberslider_callback)

#-----
smoother_callback = CustomJS(args=dict(spectra=cds_spectra, ifiberslider=ifiberslider), code="""
    ifiber = ifiberslider.value
    nsmooth = cb_obj.value
    for (var i=0; i<spectra.length; i++) {
        var data = spectra[i].data
        var flux = data['flux']
        var newflux = data['flux'+ifiber]
        for (var j=nsmooth; j<flux.length-nsmooth; j++) {
            flux[j] = 0.0
            for (var k=-nsmooth; k<=nsmooth; k++) {
                flux[j] = flux[j] + newflux[j+k] / (2*nsmooth+1)
            }
        }
        spectra[i].change.emit()
    }
""")
smootherslider = Slider(start=1, end=21, value=1, step=1, title='Smooth')
smootherslider.js_on_change('value', smoother_callback)

### bk.show(bk.Column(fig, widgetbox(zslider), widgetbox(ifiberslider)))
bk.show(bk.Column(fig, widgetbox(ifiberslider), widgetbox(smootherslider)))

#--- DEBUG ---
import IPython
IPython.embed()
#--- DEBUG ---