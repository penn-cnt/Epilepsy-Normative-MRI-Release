import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.lines as mlines
from matplotlib.gridspec import GridSpec
plt.rcParams['font.size']=8; plt.rcParams['font.family']='sans-serif'
import sys; sys.path.insert(0,'.')
from common import tract_order

long=pd.read_csv('harmonization_table_long.csv'); age=pd.read_csv('age_diag_2a.csv')
VERS=['ComBat as run','No ComBat','ComBat + age/sex, ref-only']
STY={'ComBat as run':('#2d6a9f','o'),'No ComBat':('#c4432b','s'),'ComBat + age/sex, ref-only':('#3f8f5f','^')}
comps=['Epilepsy vs. Functional Seizures','Focal vs. Generalized','Drug-Resistant vs. Drug-Responsive']
labels=[f'{t} {m}' for t in tract_order for m in ['FA','MD']]

fig=plt.figure(figsize=(7.6,8.2))
gs=GridSpec(2,2,hspace=0.30,wspace=0.46,figure=fig)
pos=[gs[0,0],gs[0,1],gs[1,0]]
_axlist=[]
for pi,comp in enumerate(comps):
    ax=fig.add_subplot(pos[pi]); d=long[long.Comparison==comp]
    for yi,lab in enumerate(labels):
        t,m=lab.rsplit(' ',1)
        if yi%2==0: ax.axhspan(yi-.5,yi+.5,color='#fafafa',zorder=0)
        for v in VERS:
            r=d[(d.Tract==t)&(d.Metric==m)&(d.Version==v)]
            if not len(r): continue
            r=r.iloc[0]; c,mk=STY[v]
            ax.scatter(r.AUC,yi,color=c,marker=mk,s=22,zorder=3,
                       edgecolor='black' if r.q<0.05 else 'none',linewidth=.9 if r.q<0.05 else 0)
    ax.axvline(0.5,ls='--',c='#aaa',lw=.9,zorder=1)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels,fontsize=5.8)
    ax.set_ylim(-.6,len(labels)-.4); ax.set_xlim(0.45,0.82)
    ax.set_title(comp,fontsize=8.0,loc='left',pad=4)
    _axlist.append(ax)
    ax.set_xlabel('AUC   (max(AUC, 1 − AUC))',fontsize=7.4)
    ax.spines[['top','right']].set_visible(False); ax.tick_params(labelsize=6.8)

ax=fig.add_subplot(gs[1,1])
bins=['22-35','36-50','51-65','66+']
for v in VERS:
    s=age[(age.version_label==v)&(age.scalar=='dti_md')].set_index('age_bin').reindex(bins)
    c,mk=STY[v]; ax.plot(range(len(bins)),s.median_z.values,marker=mk,color=c,lw=1.6,ms=6,label=v)
ax.axhline(0,ls='--',c='#aaa',lw=.9); ax.axvline(0.5,ls=':',c='#c4432b',lw=1.3)
ax.annotate('HCP-YA → HCP-Aging',xy=(0.56,ax.get_ylim()[1]),fontsize=6.4,ha='left',va='top',color='#c4432b')
ax.set_xticks(range(len(bins))); ax.set_xticklabels(bins,fontsize=7)
ax.set_xlabel('reference-subject age bin (years)',fontsize=7.4)
ax.set_ylabel('median MD z-score',fontsize=7.4)
ax.set_title('Reference age curve (MD)',fontsize=8.0,loc='left',pad=4)
_axlist.append(ax)
ax.spines[['top','right']].set_visible(False); ax.tick_params(labelsize=6.8)

h=[mlines.Line2D([],[],color=STY[v][0],marker=STY[v][1],ls='none',ms=6,label=v) for v in VERS]
h.append(mlines.Line2D([],[],color='none',marker='o',ls='none',ms=6,markeredgecolor='black',
                       markeredgewidth=.9,label='FDR q < 0.05'))
# letters in AXES coords so bbox_inches='tight' cannot shift them;
# y=1.035 sits on the single-line title's baseline (pad=4)
for _ltr,_ax,_dx in zip(['A.','B.','C.','D.'],_axlist,[-0.30,-0.30,-0.30,-0.17]):
    _ax.text(_dx,1.035,_ltr,transform=_ax.transAxes,fontsize=11,fontweight='bold',va='baseline',ha='left')
fig.legend(handles=h,loc='lower center',ncol=4,frameon=False,fontsize=7.0,bbox_to_anchor=(0.5,-0.012))
fig.subplots_adjust(top=0.97,bottom=0.075,left=0.135,right=0.985)
fig.savefig('supp_fig_s4_harmonization.pdf',bbox_inches='tight')
fig.savefig('supp_fig_s4_harmonization.png',dpi=120,bbox_inches='tight')
print('saved supp_fig_s4_harmonization.pdf')
