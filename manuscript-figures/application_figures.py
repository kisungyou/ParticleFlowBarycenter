"""Rebuild application figures at their manuscript display sizes.

Figure 7 is plotted from the original T4transport digit image, exported by
export_digit_image.R. Figure 9 PBMC and FASHION panels use saved processed data
and the original sampling/PCA code, exported by export_application_pca.R.
The complete raw inputs for the other published panels in Figures 8--10
are absent from the project and its submitted code archive. For those panels
the exact existing plotted areas are retained as raster layers, with new TeX
axes, titles and legends. No plotted numerical values or error bars are inferred
or replaced. PDF and PGF exports therefore contain vector typography and the
documented original raster layers, not a reconstruction of missing experiments.
"""
from pathlib import Path
import json
import numpy as np
from PIL import Image
from plot_style import TEXTWIDTH_IN, COLORS, OUTPUT, save_figure
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
mpl.rcParams['savefig.dpi'] = 600

ROOT = Path(__file__).resolve().parent
OLD = ROOT / 'legacy-panels'
DATA = ROOT / 'data'
DATA.mkdir(exist_ok=True)


def original(stem):
    return np.asarray(Image.open(OLD / (stem + '.png')).convert('RGB'))


def panel(ax, label):
    ax.set_title(label, loc='left', pad=5)


def figure7():
    img = np.loadtxt(DATA / 'digit8_image.csv', delimiter=',')
    threshold = float((DATA / 'digit8_threshold.txt').read_text())
    fig = plt.figure(figsize=(.9 * TEXTWIDTH_IN, 1.94))
    # The first two panels have no axes. The histogram and point cloud retain
    # their quantitative axes, all in the same Computer Modern type.
    axes = [fig.add_axes(x) for x in [
        [.018,.27,.205,.61], [.253,.27,.205,.61],
        [.555,.27,.175,.61], [.817,.27,.175,.61]]]
    for ax, data, lab in zip(axes[:2], [img, img>threshold], ['(A)','(B)']):
        ax.imshow(data, cmap='gray_r', interpolation='nearest')
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color('.45')
        panel(ax, lab)
    ax = axes[2]
    ax.hist(img.ravel(), bins=30, color='.78', edgecolor='.25', linewidth=.3)
    ax.axvline(threshold, color=COLORS[1], linestyle='--', linewidth=1.2)
    ax.set_xlabel('Intensity', labelpad=3)
    ax.set_ylabel('Count', labelpad=2)
    ax.set_xticks([0, .005, .01])
    ax.set_xticklabels(['0', '.005', '.01'])
    ax.set_yticks([0, 400, 800])
    panel(ax, '(C)')
    row, col = np.where(img > threshold)
    ax = axes[3]
    ax.scatter(2*col/27-1, 1-2*row/27, s=3, color=COLORS[0], linewidths=0)
    ax.set_xlim(-1.1,1.1)
    ax.set_ylim(-1.1,1.1)
    ax.set_xticks([-1,0,1])
    ax.set_yticks([-1,0,1])
    ax.set_xlabel('$x$', labelpad=3)
    ax.set_ylabel('$y$', labelpad=0)
    ax.set_aspect('equal', adjustable='box')
    panel(ax, '(D)')
    save_figure(fig, 'fig-real-digits-1')
    plt.close(fig)


def density_grayscale(rgb):
    # The original viridis scale has monotonically increasing luminance.
    # Direct grayscale preserves both the KDE and its antialiased contours.
    # Bandwidths and densities are not recomputed or inferred.
    return np.asarray(Image.fromarray(rgb).convert('L'))


def figure8():
    source = original('fig-real-digits-2')
    # Coordinates are exact panel-border locations in the 2700x1380 source.
    xbounds = [(504,1050),(1077,1623),(1650,2196)]
    ybounds = [(142,688),(808,1356)]
    width = .95*TEXTWIDTH_IN
    fig, axes = plt.subplots(2,3,figsize=(width,3.68))
    fig.subplots_adjust(left=.035,right=.995,bottom=.025,top=.90,hspace=.13,wspace=.08)
    for row,(y0,y1) in enumerate(ybounds):
        for col,(x0,x1) in enumerate(xbounds):
            ax=axes[row,col]
            raster=source[y0:y1,x0:x1].copy()
            if row == 1:
                raster=density_grayscale(raster)
                ax.imshow(raster,cmap='gray',vmin=0,vmax=255,interpolation='none')
            else:
                ax.imshow(raster,interpolation='none')
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_color('.45')
            if row==0:
                ax.set_title(r'$m=%d$' % [40,80,160][col],pad=5)
    fig.text(.005,.91,'(A)',ha='left',va='bottom')
    fig.text(.005,.455,'(B)',ha='left',va='bottom')
    save_figure(fig,'fig-real-digits-2')
    plt.close(fig)


def figure9():
    import csv
    source=original('fig-real-clustering-1')
    # Existing PCA views are retained, including the original observations.
    bounds=[(212,728,243,773),(901,1494,243,773),(1742,2350,243,773)]
    xticks=[([300,439,578,718],['-100','0','100','200']),
            ([1008,1169.5,1331,1492],['-2.5','0','2.5','5']),
            ([1915,2113,2311],['-10000','0','10000'])]
    yticks=[([366,505,644],['100','0','-100']),
            ([266,395,524,653],['4','2','0','-2']),
            ([398,595],['10000','0'])]
    fig,axes=plt.subplots(1,3,figsize=(.85*TEXTWIDTH_IN,2.12))
    fig.subplots_adjust(left=.103,right=.985,bottom=.28,top=.85,wspace=.46)
    for i,(ax,(x0,x1,y0,y1)) in enumerate(zip(axes,bounds)):
        if i != 1:
            name='pbmc' if i==0 else 'fashion'
            data=np.genfromtxt(DATA/('pca_'+name+'.csv'),delimiter=',',names=True)
            with (DATA/('pca_'+name+'_palette.csv')).open() as f:
                palette={int(row['label']):row['color'] for row in csv.DictReader(f)}
            colors=[palette[int(label)] for label in data['label']]
            # Rasterize only the dense point layer, at 600 dpi, to keep PGF
            # within pdfTeX's memory limit. Typography and axes remain vector.
            ax.scatter(data['x'],data['y'],s=1.8,c=colors,alpha=.75,linewidths=0,rasterized=True)
            for dim,setlim in [('x',ax.set_xlim),('y',ax.set_ylim)]:
                lo,hi=np.min(data[dim]),np.max(data[dim])
                pad=.05*(hi-lo)
                setlim(lo-pad,hi+pad)
            ax.set_xticks([-100,0,100,200] if i==0 else [-10000,0,10000])
            ax.set_yticks([-100,0,100] if i==0 else [0,10000])
            ax.grid(True,color='.92',linewidth=.45)
            ax.set_axisbelow(True)
            ax.set_aspect('equal',adjustable='box')
            ax.set_anchor('N')
        else:
            # NEWS plot areas preserve the original coordinates and overlaps.
            ax.imshow(source[y0:y1,x0:x1],extent=(0,1,0,1),aspect='auto',interpolation='none')
            xpos,xlab=xticks[i]
            ypos,ylab=yticks[i]
            ax.set_xticks([(x-x0)/(x1-x0) for x in xpos],xlab)
            ax.set_yticks([(y1-y)/(y1-y0) for y in ypos],ylab)
            ax.set_xlim(0,1)
            ax.set_ylim(0,1)
            ax.set_box_aspect((y1-y0)/(x1-x0))
            ax.set_anchor('N')
        ax.set_xlabel('PC1',labelpad=3)
        ax.set_ylabel('PC2',labelpad=1)
        ax.tick_params(labelsize=9)
        panel(ax,'(%s)' % chr(65+i))
    save_figure(fig,'fig-real-clustering-1')
    plt.close(fig)


def figure10():
    source=original('fig-real-clustering-2')
    rows=[(136,572,[(293,731),(944,1381),(1616,2054),(2234,2671)]),
          (833,1269,[(276,713),(948,1386),(1612,2049),(2251,2689)]),
          (1530,1965,[(282,720),(933,1370),(1583,2020),(2245,2682)])]
    # Existing tick centers locate the new typography exactly relative to the
    # preserved plotted image. These are display coordinates, not raw results.
    ytick_pixels=[
        [[146.5,229,312,395,478,560.5],[138,244.5,351,458,564.5],
         [152,274.5,397,520],[169,267,365.5,463.5,562]],
        [[901.5,1004,1107,1209],[891.5,960,1029,1097,1166,1234.5],
         [899,1027,1155],[866,941,1016,1091,1166,1241]],
        [[1563,1691.5,1820,1949],[1542,1671,1801,1930.5],
         [1551,1622,1693,1764.5,1835.5,1907],[1604,1684.5,1765,1846,1927]]]
    ytick_labels=[
        [['.035','.030','.025','.020','.015','.010'],['.07','.06','.05','.04','.03'],
         ['.035','.030','.025','.020'],['65','60','55','50','45']],
        [['.0050','.0045','.0040','.0035'],['.018','.017','.016','.015','.014','.013'],
         ['.00','-.02','-.04'],['260','250','240','230','220','210']],
        [['.36','.32','.28','.24'],['.52','.48','.44','.40'],
         ['.19','.18','.17','.16','.15','.14'],['1800','1700','1600','1500','1400']]]
    fig,axes=plt.subplots(3,4,figsize=(.95*TEXTWIDTH_IN,5.33))
    fig.subplots_adjust(left=.084,right=.995,bottom=.225,top=.91,hspace=.37,wspace=.43)
    for r,(y0,y1,xbounds) in enumerate(rows):
        for c,(ax,(x0,x1)) in enumerate(zip(axes[r],xbounds)):
            raster=source[y0:y1,x0:x1]
            ax.imshow(raster,extent=(0,1,0,1),aspect='auto',interpolation='none')
            ax.set_xlim(0,1)
            ax.set_ylim(0,1)
            ax.set_xticks([.055,.5,.945],['18','20','22'] if r==1 else ['8','10','12'])
            yy=[(y1-y)/(y1-y0) for y in ytick_pixels[r][c]]
            ax.set_yticks(yy,ytick_labels[r][c])
            ax.tick_params(labelsize=9,pad=2)
            if r==0:
                ax.set_title(['ARI','NMI','Silhouette','CH'][c],pad=5)
            if r==2:
                ax.set_xlabel('$k$',labelpad=2)
        # Compact row labels avoid repeating dataset names in every panel.
        axes[r,0].text(-.40,.50,['(A) PBMC','(B) 20NEWS','(C) FASHION'][r],
                       transform=axes[r,0].transAxes,rotation=90,
                       va='center',ha='center',fontsize=11)
    # Exact styles from the original ggplot (solid, dashed, dotted, dotdash,
    # longdash and twodash). Color is redundant with the line pattern.
    original_colors=['#0072B2','#E69F00','#009E73','#CC79A7','#56B4E9','#D55E00']
    ltys=['-',(0,(4,4)),(0,(1,3)),(0,(1,3,4,3)),(0,(7,3)),(0,(2,2,6,2))]
    names=[r'$k$-means',r'Spherical $k$-means',
           r'DVQ $S=5$, $\alpha=1$',r'DVQ $S=5$, $\alpha=0.5$',
           r'DVQ $S=10$, $\alpha=1$',r'DVQ $S=10$, $\alpha=0.5$']
    handles=[Line2D([],[],color=co,ls=ls,lw=1.15,marker='o',markersize=3.5,label=la)
             for co,ls,la in zip(original_colors,ltys,names)]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.53,.013),ncol=2,
               frameon=False,fontsize=10,handlelength=3.8,columnspacing=1.0,labelspacing=.5)
    save_figure(fig,'fig-real-clustering-2')
    plt.close(fig)


if __name__=='__main__':
    import sys
    selected = set(sys.argv[1:] or ['7','8','9','10'])
    for key,fn in [('7',figure7),('8',figure8),('9',figure9),('10',figure10)]:
        if key in selected:
            fn()
            print('Exported Figure', key, flush=True)
    manifest={
      'figure7':'Original T4transport first digit-8 image; original 256-bin Otsu rule.',
      'figure8':'Original plotted particle and KDE panels retained, with monotone grayscale density recoloring and TeX labels. Complete core prototype files are not present.',
      'figure9':'PBMC and FASHION regenerated from saved processed features with the original R sampling/PCA procedure. Original NEWS plotted area retained with TeX axes because its processed features are absent. No class-identification claim is made in the revised caption.',
      'figure10':'Original plotted means, error bars and six line patterns retained, with TeX axes and matching legend. Complete 162-run metric table is not present, so no numerical data were digitized or recomputed.',
      'formats':'PGF plus pdfLaTeX-compiled PDF. Figures 8, 9B and 10 contain preserved raster plot layers and vector typography. Figures 7, 9A and 9C are plotted from saved data.',
    }
    (OUTPUT/'application-figure-provenance.json').write_text(json.dumps(manifest,indent=2)+'\n')
