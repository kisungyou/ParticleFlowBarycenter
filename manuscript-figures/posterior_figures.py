from pathlib import Path
from plot_style import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE=Path(__file__).resolve().parent
DATA=HERE/'data'
REFRESH=HERE.parent/'revision-work'/'refresh'
df=pd.read_csv(REFRESH/'wasp_metrics.csv')
ref=pd.read_csv(REFRESH/'wasp_full_reference.csv').iloc[0]

# The same KDE grids used by the existing representative posterior illustration.
fig,axs=plt.subplots(1,4,figsize=(.95*TEXTWIDTH_IN,2.12))
fig.subplots_adjust(left=.018,right=.995,bottom=.27,top=.86,wspace=.16)
for j,name in [(0,'analytic'),(1,'full_mcmc'),(3,'barycenter')]:
    g=pd.read_csv(DATA/f'posterior_{name}.csv')
    x=np.sort(g.x.unique()); y=np.sort(g.y.unique());z=g.z.to_numpy().reshape((len(y),len(x)))
    axs[j].pcolormesh(x,y,z,cmap='viridis',shading='auto',rasterized=True)
    levels=np.linspace(z.min(),z.max(),8)[1:-1]
    axs[j].contour(x,y,z,levels=levels,colors='white',linewidths=.35)
patterns=['solid','dashed','dotted','dashdot',(0,(6,2,1,2,1,2))]
for i in range(1,6):
    g=pd.read_csv(DATA/f'posterior_subset_{i}.csv')
    x=np.sort(g.x.unique());y=np.sort(g.y.unique());z=g.z.to_numpy().reshape((len(y),len(x)))
    axs[2].contour(x,y,z,levels=np.linspace(0,z.max(),8)[1:-1],colors=[COLORS[i-1]],linestyles=[patterns[i-1]],linewidths=.65)
for j,ax in enumerate(axs):
    ax.set(xlim=(.63,1.13),ylim=(.77,1.27),aspect='equal')
    ax.set_axis_off();ax.set_title(f'({chr(65+j)})',loc='left',pad=5)
fig.legend([Line2D([0],[0],color=COLORS[i],ls=patterns[i],lw=1.15) for i in range(5)], [f'Subset {i}' for i in range(1,6)],loc='lower center',ncol=5,frameon=False,handlelength=2.0,columnspacing=.7,bbox_to_anchor=(.5,.01))
save_figure(fig,'fig-sim-normal-1');plt.close(fig)

markers={50:'o',100:'^',200:'s'}
# Categorical positions are dodged horizontally solely to reveal overlapping summaries.
# Means and standard deviations remain those of the saved runs, with all runs included.
def metric_panel(ax,metric,label,j):
    for si,m in enumerate([50,100,200]):
        for ai,alpha in enumerate([1,.5]):
            s=df[(df.support==m)&(df.alpha==alpha)].groupby('nsplit')[metric].agg(['mean','std']).reindex([5,10,20])
            offset=(2*si+ai-2.5)*.047
            x=np.arange(3)+offset
            ax.errorbar(x,s['mean'],yerr=s['std'],color=COLORS[si],ls='-' if alpha==1 else '--',
                        marker=markers[m],markersize=4.0 if alpha==1 else 4.8,
                        markerfacecolor=COLORS[si] if alpha==1 else 'white',markeredgewidth=.85,
                        linewidth=1.05,elinewidth=.65,capsize=1.6,zorder=3+ai)
    if metric!='runtime':ax.axhline(ref[metric],color='.2',ls=':',lw=1.05,zorder=1)
    ax.set_title(f'({chr(65+j)})',loc='left',pad=6)
    ax.set_xticks([0,1,2],['5','10','20']);ax.set_xlim(-.27,2.27)
    ax.set_xlabel('Subsets');ax.set_ylabel(label,labelpad=4)
    if metric in ['coef_inclusion','mean_inclusion']:
        ax.set_ylim(-.045,1.08);ax.set_yticks([0,.5,1])
    ax.tick_params(axis='both',pad=2)
    if metric not in ['coef_inclusion','mean_inclusion']: ax.yaxis.set_major_locator(plt.MaxNLocator(4))
    if metric=='mean_rmse':
        from matplotlib.ticker import FuncFormatter
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v,pos: f'{100*v:g}'))
        ax.set_ylabel(r'RMSE ($10^{-2}$)',labelpad=4)

def legend(fig):
    h=[Line2D([0],[0],color=COLORS[i],marker=markers[m],ls='None',markerfacecolor='white',ms=5) for i,m in enumerate([50,100,200])]
    h+=[Line2D([0],[0],color='.2',ls='-',marker='o',ms=4),Line2D([0],[0],color='.2',ls='--',marker='o',mfc='white',ms=4),Line2D([0],[0],color='.2',ls=':')]
    labs=[r'$m=50$',r'$m=100$',r'$m=200$',r'$\alpha=1$',r'$\alpha=0.5$','Full-data MCMC']
    fig.legend(h,labs,loc='lower center',bbox_to_anchor=(.5,.005),ncol=3,frameon=False,columnspacing=1.3,handlelength=2.4,labelspacing=.35)

fig,axs=plt.subplots(1,3,figsize=(.95*TEXTWIDTH_IN,2.75))
fig.subplots_adjust(left=.12,right=.995,bottom=.34,top=.87,wspace=.68)
for j,(v,l) in enumerate([('semi_w2',r'$W_2$ error'),('mean_error','Mean error'),('cov_error','Covariance error')]):metric_panel(axs[j],v,l,j)
legend(fig);save_figure(fig,'fig-sim-normal-2');plt.close(fig)

fig,axs=plt.subplots(1,4,figsize=(.99*TEXTWIDTH_IN,2.75))
fig.subplots_adjust(left=.12,right=.995,bottom=.34,top=.87,wspace=.87)
for j,(v,l) in enumerate([('mean_rmse','RMSE'),('mean_inclusion','Inclusion'),('coef_inclusion','Inclusion'),('runtime','Time (s)')]):metric_panel(axs[j],v,l,j)
legend(fig);save_figure(fig,'fig-sim-normal-3');plt.close(fig)
print('Posterior figures4,5,6 exported as PGF and PDF.')
