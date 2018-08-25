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

plot_width=800
fig = bk.figure(height=350, width=plot_width)
colors = dict(b='#1f77b4', r='#d62728', z='maroon')
for spec in cds_spectra:
    fig.line('plotwave', 'plotflux', source=spec, line_color=colors[spec.name])

zslider_callback  = CustomJS(
    args=dict(spectra=cds_spectra),
    code="""
    var z = cb_obj.value
    for (var i=0; i<spectra.length; i++) {
        var data = spectra[i].data
        var origwave = data['origwave']
        var plotwave = data['plotwave']
        for (var j=0; j<plotwave.length; j++) {
            plotwave[j] = origwave[j] / (1+z)
        }
        spectra[i].change.emit()
    }
    """)

zslider = Slider(start=0, end=2.0, value=0.0, step=0.01, title='Redshift')
zslider.js_on_change('value', zslider_callback)

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

#-----
#- Add navigation buttons
navigation_button_width = 30
prev_button = Button(label="<", width=navigation_button_width)
next_button = Button(label=">", width=navigation_button_width)

prev_callback = CustomJS(
    args=dict(ifiberslider=ifiberslider),
    code="""
    if(ifiberslider.value>0) {
        ifiberslider.value--
    }
    """)
next_callback = CustomJS(
    args=dict(ifiberslider=ifiberslider, nspec=fr.nspec),
    code="""
    if(ifiberslider.value<nspec+1) {
        ifiberslider.value++
    }
    """)

prev_button.js_on_event('button_click', prev_callback)
next_button.js_on_event('button_click', next_callback)

#-----
### bk.show(bk.Column(fig, widgetbox(zslider), widgetbox(ifiberslider)))
slider_width = plot_width - 2*navigation_button_width
navigator = bk.Row(prev_button, next_button, widgetbox(ifiberslider, width=slider_width))
bk.show(bk.Column(
    fig,
    navigator,
    widgetbox(smootherslider),
    widgetbox(zslider),
    ))

#--- DEBUG ---
# import IPython
# IPython.embed()
#--- DEBUG ---