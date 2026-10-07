"""Restyle exact saved runtime plot areas without reconstructing missing values.
The complete core runtime table is absent from both code archives. The five
high-resolution plot interiors are retained unchanged; every label/tick is TeX.
Pixel coordinates below locate existing axes and ticks, not estimated run values.
"""
from plot_style import *
from PIL import Image
import matplotlib.pyplot as plt
import json
HERE=Path(__file__).resolve().parent
src=HERE/'legacy-panels'/'fig-sim-gauss-3.png'
im=Image.open(src).convert('RGB')
# x0,y0,x1,y1; original pixel tick centers; labels transcribed from existing axes.
panels=[
 ((214,119,907,812),[291.5,425,558,691.5,824.5],['10','50','100','200','500'],[215,417.5,620],['1.5','1.0','0.5'],'Barycenter size'),
 ((1079,119,1772,812),[1175,1340,1505,1670],['50','100','200','500'],[231,348,465.5,582.5,699.5],['0.5','0.4','0.3','0.2','0.1'],'Target size'),
 ((1943,119,2636,812),[2040,2205,2370,2535],['2','4','8','16'],[188,394,600,805.5],['1.5','1.0','0.5','0.0'],'Measures'),
 ((194,1076,887,1768),[321.5,538,754.5],['2','10','50'],[1232,1410.5,1589,1767.5],['0.18','0.15','0.12','0.09'],'Dimension'),
 ((1059,1076,1752,1768),[1245,1560],['0.5','1'],[1097.5,1272,1446,1620],['1.2','1.0','0.8','0.6'],'Step size'),
]
fig,axs=plt.subplots(2,3,figsize=(.9*TEXTWIDTH_IN,4.03))
fig.subplots_adjust(left=.11,right=.994,bottom=.14,top=.93,wspace=.87,hspace=.89)
for j,(box,xt,xl,yt,yl,xlabel) in enumerate(panels):
 ax=axs.flat[j]; x0,y0,x1,y1=box
 ax.imshow(im.crop(box),extent=(0,1,0,1),origin='upper',aspect='auto',interpolation='none')
 ax.set_xticks([(v-x0)/(x1-x0) for v in xt],xl,rotation=40,ha='right')
 ax.set_yticks([(y1-v)/(y1-y0) for v in yt],yl)
 ax.set(xlim=(0,1),ylim=(0,1),xlabel=xlabel,ylabel='Time (s)')
 ax.tick_params(pad=2);ax.set_title(f'({chr(65+j)})',loc='left',pad=5)
axs.flat[5].set_axis_off()
save_figure(fig,'fig-sim-gauss-3')
(OUTPUT/'runtime-figure-provenance.json').write_text(json.dumps({'figure':2,'source':'legacy-panels/fig-sim-gauss-3.png','mode':'Exact existing plot interiors with replacement TeX labels','reason':'Full assembled_runtime_scaling-core.csv absent from available experiment outputs','numerical_values_reconstructed':False,'panels':panels},indent=2))
print('Runtime figure exported with unchanged plot interiors.')
