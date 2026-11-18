import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.lines as mlines
plt.rcParams['font.size']=8; plt.rcParams['font.family']='sans-serif'
import sys; sys.path.insert(0,'.')
from common import tract_order

long=pd.read_csv('harmonization_table_long.csv')
age=pd.read_csv('age_diag_2a.csv')
VERS=['ComBat as run','No ComBat','ComBat + age/sex, ref-only']
STY={'ComBat as run':('#2d6a9f','o'),'No ComBat':('#c4432b','s'),'ComBat + age/sex, ref-only':('#3f8f5f','^')}
comps=['Epilepsy vs. Functional Seizures','Focal vs. Generalized','Drug-Resistant vs. Drug-Responsive']
labels=[f'{t} {m}' for t in tract_order for m in ['FA','MD']]

fig=plt.figure(figsize=(13.5,7.4))
gs=fig.add_gridspec(1,4,width_ratios=[1,1,1,0.85],wspace=0.42)
for pi,comp in enumerate(comps):
    ax=fig.add_subplot(gs[0,pi]); d=long[long.Comparison==comp]
    for yi,lab in enumerate(labels):
        t,m=lab.rsplit(' ',1)
        for v in VERS:
            r=d[(d.Tract==t)&(d.Metric==m)&(d.Version==v)]
            if not len(r): continue
            r=r.iloc[0]; c,mk=STY[v]
            ax.scatter(r.AUC,yi,color=c,marker=mk,s=34,zorder=3,
                       edgecolor='black' if r.q<0.05 else 'none',linewidth=1.1 if r.q<0.05 else 0)
    ax.axvline(0.5,ls='--',c='#aaa',lw=.9,zorder=1)
    for yi in range(len(labels)):
        ax.axhline(yi,c='#f0f0f0',lw=.6,zorder=0)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels if pi==0 else ['']*len(labels),fontsize=6.6)
    ax.set_ylim(-.7,len(labels)-.3); ax.set_xlim(0.45,0.80)
    ax.set_xlabel('AUC  (max(AUC, 1−AUC))',fontsize=8)
    ax.set_title(f"{'ABC'[pi]}. {comp}",fontsize=8.6,loc='left')
    ax.spines[['top','right']].set_visible(False); ax.tick_params(labelsize=7)

ax=fig.add_subplot(gs[0,3])
bins=['22-35','36-50','51-65','66+']
for v in VERS:
    s=age[(age.version_label==v)&(age.scalar=='dti_md')].set_index('age_bin').reindex(bins)
    c,mk=STY[v]
    ax.plot(range(len(bins)),s.median_z.values,marker=mk,color=c,lw=1.5,ms=6,label=v)
ax.axhline(0,ls='--',c='#aaa',lw=.9); ax.axvline(0.5,ls=':',c='#c4432b',lw=1.2)
ax.annotate('HCP-YA → HCP-Aging',xy=(0.5,ax.get_ylim()[1]),fontsize=6.4,ha='center',va='top',color='#c4432b')
ax.set_xticks(range(len(bins))); ax.set_xticklabels(bins,fontsize=7.5)
ax.set_xlabel('reference-subject age bin (y)',fontsize=8)
ax.set_ylabel('median MD z-score',fontsize=8)
ax.set_title('D. Reference age curve (MD)',fontsize=8.6,loc='left')
ax.spines[['top','right']].set_visible(False); ax.tick_params(labelsize=7)

h=[mlines.Line2D([],[],color=STY[v][0],marker=STY[v][1],ls='none',ms=6,label=v) for v in VERS]
h.append(mlines.Line2D([],[],color='none',marker='o',ls='none',ms=6,markeredgecolor='black',
                       markeredgewidth=1.1,label='FDR q < 0.05'))
fig.legend(handles=h,loc='lower center',ncol=4,frameon=False,fontsize=7.6,bbox_to_anchor=(0.5,-0.015))
fig.suptitle('Supplementary Figure: tractometry discrimination under three harmonization schemes',
             fontsize=10,x=0.5,y=0.985)
fig.savefig('supp_fig_harmonization.pdf',bbox_inches='tight')
print("saved supp_fig_harmonization.pdf")
