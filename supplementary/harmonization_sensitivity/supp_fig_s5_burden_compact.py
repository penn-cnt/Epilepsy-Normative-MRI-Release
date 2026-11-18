"""Figure S5: whole-brain deviation burden + hippocampal positive control.
Reuses the composite_figures.ipynb conventions: cohort = 355, morph from centilebrain_zscores/, spelling fix applied,
AUC = max(AUC,1-AUC) via Mann-Whitney with 5000-bootstrap CI (seed 42), BH-FDR."""
import pandas as pd, numpy as np, sys, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests
sys.path.insert(0,'.')
from common import BASE, META, mwu_auc, bootstrap_auc_ci, fold_auc
plt.rcParams['font.size']=8; plt.rcParams['font.family']='sans-serif'

cb_dir=f'{BASE}/centilebrain_zscores'
ct=pd.concat([pd.read_csv(f'{cb_dir}/zscore_CorticalThickness_male.csv'),
              pd.read_csv(f'{cb_dir}/zscore_CorticalThickness_female.csv')],ignore_index=True)
sv=pd.concat([pd.read_csv(f'{cb_dir}/zscore_SubcorticalVolume_male.csv'),
              pd.read_csv(f'{cb_dir}/zscore_SubcorticalVolume_female.csv')],ignore_index=True)
ct=ct.rename(columns=lambda c:c.replace('_entorhil_','_entorhinal_').replace('_supramargil_','_supramarginal_'))
both=set(ct.SubjectID)&set(sv.SubjectID); ct=ct[ct.SubjectID.isin(both)]; sv=sv[sv.SubjectID.isin(both)]
sv_z=[c for c in sv.columns if c not in ['SITE','SubjectID','Vendor','FreeSurfer_Version','age']]
morph=ct.merge(sv[['SubjectID']+sv_z],on='SubjectID',how='inner')
meta=pd.read_csv(META).rename(columns={'record_id':'SubjectID'})
morph=morph.merge(meta[['SubjectID','group_1','group_2','lesion_mts','lesion_left','lesion_right']],
                  on='SubjectID',how='left').drop_duplicates('SubjectID',keep='first')
tract_subs=set(pd.read_csv(f'{BASE}/node_zscores.csv',usecols=['subject'])['subject'].unique())
morph=morph[morph.SubjectID.isin(tract_subs)].copy()
print(f"cohort n={len(morph)}")

roi=[c for c in ct.columns if c not in ['SITE','SubjectID','Vendor','FreeSurfer_Version','age','sex']]+sv_z
print(f"ROIs counted: {len(roi)}")
morph['burden']=(morph[roi].abs()>2).sum(axis=1)
epi=morph[morph.group_2.isin(['drug-resistant','drug-responsive'])]; fs=morph[morph.group_1=='pnee']
a,b=epi.burden.values,fs.burden.values
auc,p=mwu_auc(a,b); lo,hi=bootstrap_auc_ci(a,b); A,L,H=fold_auc(auc,lo,hi)
print(f"(A) burden AUC={A:.3f} ({L:.3f}-{H:.3f}) p={p:.4f} | epilepsy n={len(a)} median={np.median(a):.0f}, FS n={len(b)} median={np.median(b):.0f}")

mts=morph[morph.lesion_mts==1].copy()
mts['side']=np.where(mts.lesion_left==1,'L','R')
mts['scler']=np.where(mts.side=='L',mts.Lhippo,mts.Rhippo)
mts['contra']=np.where(mts.side=='L',mts.Rhippo,mts.Lhippo)
st,pw=wilcoxon(mts.scler,mts.contra)
print(f"(B) MTS n={len(mts)} (L={(mts.side=='L').sum()},R={(mts.side=='R').sum()}) "
      f"sclerotic median={mts.scler.median():.2f} contra={mts.contra.median():.2f} Wilcoxon p={pw:.4f}")
print(f"    sclerotic z<-2: {(mts.scler<-2).sum()}/{len(mts)}")

morph['hippo_mean']=morph[['Lhippo','Rhippo']].mean(axis=1)
morph['hippo_ai']=(morph.Lhippo-morph.Rhippo)/((morph.Lhippo.abs()+morph.Rhippo.abs())/2)
grp={'Epilepsy':epi.index,'Functional\nseizures':fs.index,'MTS':mts.index}
print("(C) by group:")
for g,idx in grp.items():
    d=morph.loc[idx]
    print(f"    {g.replace(chr(10),' ')}: n={len(d)} mean hippo z={d.hippo_mean.median():.3f} |AI| median={d.hippo_ai.abs().median():.3f}")

fig=plt.figure(figsize=(7.0,6.0))
from matplotlib.gridspec import GridSpec
_gs=GridSpec(2,4,height_ratios=[1.5,0.92],hspace=0.38,wspace=0.95,figure=fig)
axes=[fig.add_subplot(_gs[0,0:2]),fig.add_subplot(_gs[0,2:4]),fig.add_subplot(_gs[1,0:3])]
ax=axes[0]
for i,(lbl,d,c) in enumerate([('Epilepsy',a,'#2d6a9f'),('Functional\nseizures',b,'#c4432b')]):
    x=np.random.default_rng(42).normal(i,0.07,len(d))
    ax.scatter(x,d,s=15,alpha=.55,color=c,edgecolor='none')
    ax.hlines(np.median(d),i-.26,i+.26,color='black',lw=2,zorder=4)
ax.set_xticks([0,1]); ax.set_xticklabels(['Epilepsy','Functional\nseizures'])
ax.set_ylabel(f'regions with |z| > 2  (of {len(roi)})')
ax.set_title('Whole-brain deviation burden',fontsize=8.6,loc='left',pad=14)
ax.text(0.0,1.005,f'AUC {A:.2f} (95% CI {L:.2f}-{H:.2f})',transform=ax.transAxes,fontsize=7.4,va='bottom',ha='left',color='#444')

ax=axes[1]
for i,r in enumerate(mts.itertuples()):
    ax.plot([0,1],[r.scler,r.contra],color='#999',lw=.8,zorder=1)
ax.scatter([0]*len(mts),mts.scler,s=28,color='#c4432b',zorder=3,label='sclerotic')
ax.scatter([1]*len(mts),mts.contra,s=28,color='#2d6a9f',zorder=3,label='contralateral')
ax.axhline(-2,ls='--',c='#555',lw=1)
ax.set_xticks([0,1]); ax.set_xticklabels(['Sclerotic','Contralateral']); ax.set_xlim(-.35,1.35)
ax.set_ylabel('hippocampal volume z-score')
ax.set_title(f'MTS positive control (n={len(mts)})',fontsize=8.6,loc='left',pad=14)
ax.text(0.0,1.005,f'Wilcoxon p = {pw:.4f}',transform=ax.transAxes,fontsize=7.4,va='bottom',ha='left',color='#444')
ax.legend(fontsize=6.8,frameon=False,loc='upper left',bbox_to_anchor=(-0.02,1.0),handletextpad=.4)

ax=axes[2]; labs=['Epilepsy','Functional\nseizures','MTS']; h=0.35
for j,(metric,col,c) in enumerate([('mean hippocampal z','hippo_mean','#2d6a9f'),
                                   ('|asymmetry index|','hippo_ai','#3f8f5f')]):
    vals=[morph.loc[idx,col].abs().median() if j==1 else morph.loc[idx,col].median() for idx in grp.values()]
    ax.barh(np.arange(len(labs))+(j-0.5)*h,vals,h,color=c,label=metric)
ax.axvline(0,c='#333',lw=.8); ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs,fontsize=7.5)
ax.set_xlabel('median value'); ax.set_xlim(right=1.52)
ax.legend(fontsize=6.8,frameon=False,loc='upper right',bbox_to_anchor=(1.0,0.97),handletextpad=.4,labelspacing=.35)
ax.set_title('Hippocampal z and asymmetry by group',fontsize=8.6,loc='left')
for x in axes: x.spines[['top','right']].set_visible(False); x.tick_params(labelsize=7.5)
fig.canvas.draw()
_pA=axes[0].get_position(); _pB=axes[1].get_position(); _pC=axes[2].get_position()
# make C start where A starts and end where B ends
axes[2].set_position([_pA.x0,_pC.y0,(_pB.x1-_pA.x0),_pC.height])
fig.canvas.draw()
_xA=axes[0].get_position().x0; _xB=axes[1].get_position().x0
# Letters sit on the FIRST title line. A title line is ~0.026 fig-units at this
# size; a 2-line title therefore starts one extra line above the axes top.
for _ltr,_ax,_fx in zip('ABC',axes,[_xA,_xB,_xA]):
    _y=_ax.get_position().y1+(0.0335 if _ax is not axes[2] else 0.016)
    fig.text(_fx-0.055,_y,_ltr+'.',fontsize=11,fontweight='bold',va='baseline',ha='left')
fig.savefig('supp_fig_s5_burden_compact.pdf',bbox_inches='tight')
fig.savefig('supp_fig_s5_burden_compact.png',dpi=120,bbox_inches='tight')
pd.DataFrame([{'panel':'A','metric':'burden AUC','value':f"{A:.3f} ({L:.3f}-{H:.3f})",'p':round(p,4),
               'epilepsy_median':float(np.median(a)),'fs_median':float(np.median(b)),'n_roi':len(roi)},
              {'panel':'B','metric':'MTS sclerotic vs contra','value':f"{mts.scler.median():.2f} vs {mts.contra.median():.2f}",
               'p':round(pw,4),'n':len(mts)}]).to_csv('supp_fig_s5_burden_values.csv',index=False)
print("saved supp_fig_s5_burden_hippocampus.pdf")
