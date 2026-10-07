"""Replot Figures 1 and 3 from saved revision results, without refitting.

The fixed canvases match their manuscript inclusion widths. All data markers
remain at their true x coordinates. Nested open/filled markers and distinguishable
line patterns expose coincident series without jittering or changing summaries.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from plot_style import TEXTWIDTH_IN, COLORS, OUTPUT, save_figure
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, FuncFormatter

BASE=Path(__file__).resolve().parents[1]
WORK=BASE/'revision-work'
OUT=OUTPUT
OUT.mkdir(parents=True, exist_ok=True)

def summarize(frame,groups,columns):
    result=frame.groupby(groups)[columns].agg(['mean','std','count'])
    result.columns=['_'.join(col) for col in result.columns]
    return result.reset_index()

def series(ax,x,mean,sd,*,label,color,marker,ms,filled,ls,capsize,zorder):
    """SD, not SE/CI. No horizontal display displacement is applied."""
    return ax.errorbar(x,mean,yerr=sd,label=label,color=color,
       marker=marker,markersize=ms,markerfacecolor=color if filled else 'white',
       markeredgecolor=color,markeredgewidth=.8,linestyle=ls,linewidth=1.05,
       capsize=capsize,capthick=.65,elinewidth=.65,zorder=zorder)

def tidy(ax,panel,ylabel,xvalues):
    ax.set_title(panel,loc='left',pad=5)
    ax.set_xlabel('Support size',labelpad=3)
    ax.set_ylabel(ylabel,labelpad=4)
    ax.set_xticks(np.arange(len(xvalues)),[str(v) for v in xvalues])
    ax.tick_params(axis='both',pad=2)
    ax.set_xlim(-.28,len(xvalues)-.72)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4,min_n_ticks=3))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{x:g}'))

def figure_one():
    d=pd.read_csv(WORK/'refresh/gaussian_resolution.csv')
    columns=['semi_w2','mean_error','cov_error']
    s=summarize(d,['support','alpha'],columns)
    s.to_csv(OUT/'fig1_point_summaries.csv',index=False)
    support=[10,25,50,100,200,500];x=np.arange(len(support))
    fig,axs=plt.subplots(1,3,figsize=(.95*TEXTWIDTH_IN,2.85))
    fig.subplots_adjust(left=.092,right=.982,bottom=.265,top=.875,wspace=.43)
    for ax,col,panel,ylabel in zip(axs,columns,['(A)','(B)','(C)'],
             [r'$W_2$ error','Mean error','Bures error']):
        # Larger open squares are drawn before smaller filled circles so both
        # method symbols remain visible even when their ordinates coincide.
        for alpha,style in [(.5,dict(label=r'$\alpha=0.5$',color=COLORS[1],marker='s',ms=5.8,filled=False,ls=(0,(4,2)),capsize=3.,zorder=2)),
                             (1.,dict(label=r'$\alpha=1$',color=COLORS[0],marker='o',ms=2.8,filled=True,ls='-',capsize=1.5,zorder=3))]:
            g=s[s.alpha==alpha].set_index('support').loc[support]
            series(ax,x,g[col+'_mean'],g[col+'_std'],**style)
        tidy(ax,panel,ylabel,support)
        ax.set_ylim(bottom=0)
    handles,labels=axs[0].get_legend_handles_labels()
    fig.legend([handles[1],handles[0]],[labels[1],labels[0]],ncol=2,
       loc='lower center',bbox_to_anchor=(.5,.017),frameon=False,
       handlelength=2.5,columnspacing=2.,handletextpad=.6)
    save_figure(fig,'fig-sim-gauss-2')
    plt.close(fig)

def figure_three():
    d=pd.read_csv(WORK/'gaussian/paired_results.csv')
    columns=['semi_w2','cov_error','runtime_sec']
    s=summarize(d,['support','method'],columns)
    s.to_csv(OUT/'fig3_point_summaries.csv',index=False)
    supports=[50,100,200];x=np.arange(len(supports))
    styles=[
      ('r_exact',dict(label='R/T4 exact MM',color=COLORS[0],marker='s',ms=5.8,filled=False,ls='-',capsize=3.,zorder=3)),
      ('pot_exact',dict(label='Python/POT exact',color=COLORS[1],marker='o',ms=2.8,filled=True,ls=(0,(4,2)),capsize=1.5,zorder=4)),
      ('sinkhorn_0.1',dict(label=r'Sinkhorn, $\varepsilon=0.1s^2$',color=COLORS[2],marker='^',ms=4.7,filled=False,ls=':',capsize=2.,zorder=2)),
      ('sinkhorn_0.3',dict(label=r'Sinkhorn, $\varepsilon=0.3s^2$',color=COLORS[3],marker='D',ms=4.2,filled=False,ls='-.',capsize=2.,zorder=2)),
      ('sinkhorn_1',dict(label=r'Sinkhorn, $\varepsilon=s^2$',color=COLORS[4],marker='v',ms=4.5,filled=True,ls=(0,(6,2,1,2)),capsize=2.,zorder=2))]
    fig,axs=plt.subplots(1,3,figsize=(.99*TEXTWIDTH_IN,3.12))
    fig.subplots_adjust(left=.087,right=.993,bottom=.34,top=.89,wspace=.45)
    for ax,col,panel,ylabel in zip(axs,columns,['(A)','(B)','(C)'],
             [r'$W_2$ error','Bures error','Runtime (s)']):
        for method,style in styles:
            g=s[s.method==method].set_index('support').loc[supports]
            series(ax,x,g[col+'_mean'],g[col+'_std'],**style)
        tidy(ax,panel,ylabel,supports)
        if col=='runtime_sec':
            ax.set_yscale('log')
        else:ax.set_ylim(bottom=0)
    handles,labels=axs[0].get_legend_handles_labels()
    fig.legend(handles,labels,ncol=3,loc='lower center',bbox_to_anchor=(.5,.016),
       frameon=False,handlelength=2.2,columnspacing=1.15,handletextpad=.5,
       labelspacing=.5,borderaxespad=0.)
    save_figure(fig,'fig-sim-gauss-4')
    plt.close(fig)

def audit():
    d=pd.read_csv(WORK/'refresh/gaussian_resolution.csv')
    half=d[d.alpha==.5].set_index(['rep','support'])
    full=d[d.alpha==1].set_index(['rep','support'])
    info={'source':'revision-work/refresh/gaussian_resolution.csv',
          'number_of_runs':len(d),'all_runs_converged':bool((d.status=='converged').all()),
          'summaries':'arithmetic mean and sample SD (ddof=1), no changes to saved data',
          'horizontal_offsets':False,'textwidth_in':TEXTWIDTH_IN,
          'figure1_width_in':.95*TEXTWIDTH_IN,'figure3_width_in':.99*TEXTWIDTH_IN,
          'mean_identity':'Balanced transport plans imply sum_i v_i M_i = sum_n lambda_n sum_j w_nj x_nj. The full update reaches this empirical input mean immediately; damping approaches the same mean.'}
    for col in ['semi_w2','mean_error','cov_error','F','niter','runtime']:
        delta=half[col]-full[col]
        info[col]={'max_abs_paired_difference':float(delta.abs().max()),'mean_paired_difference':float(delta.mean()),
                   'half_mean':float(half[col].mean()),'full_mean':float(full[col].mean()),
                   'half_median':float(half[col].median()),'full_median':float(full[col].median())}
    info['support_array_audit']='data/gaussian_support_audit.csv, computed by reading saved RDS support arrays; no optimization rerun'
    (OUT/'gaussian_alpha_audit.json').write_text(json.dumps(info,indent=2))

if __name__=='__main__':
    audit();figure_one();figure_three()
