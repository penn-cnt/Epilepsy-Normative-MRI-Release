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

fig,axes=plt.subplots(3,1,figsize=(5.4,10.4))
ax=axes[0]
for i,(lbl,d,c) in enumerate([('Epilepsy',a,'#2d6a9f'),('Functional\nseizures',b,'#c4432b')]):
    x=np.random.default_rng(42).normal(i,0.07,len(d))
    ax.scatter(x,d,s=15,alpha=.55,color=c,edgecolor='none')
    ax.hlines(np.median(d),i-.26,i+.26,color='black',lw=2,zorder=4)
ax.set_xticks([0,1]); ax.set_xticklabels(['Epilepsy','Functional\nseizures'])
ax.set_ylabel(f'regions with |z| > 2  (of {len(roi)})')
ax.set_title(f'A. Whole-brain deviation burden\nAUC {A:.2f} (95% CI {L:.2f}-{H:.2f})',fontsize=8.6,loc='left')

ax=axes[1]
for i,r in enumerate(mts.itertuples()):
    ax.plot([0,1],[r.scler,r.contra],color='#999',lw=.8,zorder=1)
ax.scatter([0]*len(mts),mts.scler,s=28,color='#c4432b',zorder=3,label='sclerotic')
ax.scatter([1]*len(mts),mts.contra,s=28,color='#2d6a9f',zorder=3,label='contralateral')
ax.axhline(-2,ls='--',c='#555',lw=1)
ax.set_xticks([0,1]); ax.set_xticklabels(['Sclerotic','Contralateral']); ax.set_xlim(-.35,1.35)
ax.set_ylabel('hippocampal volume z-score')
ax.set_title(f'B. MTS positive control (n={len(mts)})\nWilcoxon p = {pw:.4f}',fontsize=8.6,loc='left')
ax.legend(fontsize=7,frameon=False,loc='lower left')

ax=axes[2]; labs=[l.replace(chr(10),' ') for l in grp]; h=0.35
for j,(metric,col,c) in enumerate([('mean hippocampal z','hippo_mean','#2d6a9f'),
                                   ('|asymmetry index|','hippo_ai','#3f8f5f')]):
    vals=[morph.loc[idx,col].abs().median() if j==1 else morph.loc[idx,col].median() for idx in grp.values()]
    ax.barh(np.arange(len(labs))+(j-0.5)*h,vals,h,color=c,label=metric)
ax.axvline(0,c='#333',lw=.8); ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs,fontsize=7.5)
ax.set_xlabel('median value'); ax.legend(fontsize=7,frameon=False,loc='lower right')
ax.set_title('C. Hippocampal z and asymmetry by group',fontsize=8.6,loc='left')
for x in axes: x.spines[['top','right']].set_visible(False); x.tick_params(labelsize=7.5)
fig.tight_layout(h_pad=2.2)
fig.savefig('supp_fig_s5_burden_portrait.pdf',bbox_inches='tight')
fig.savefig('supp_fig_s5_burden_portrait.png',dpi=120,bbox_inches='tight')
pd.DataFrame([{'panel':'A','metric':'burden AUC','value':f"{A:.3f} ({L:.3f}-{H:.3f})",'p':round(p,4),
               'epilepsy_median':float(np.median(a)),'fs_median':float(np.median(b)),'n_roi':len(roi)},
              {'panel':'B','metric':'MTS sclerotic vs contra','value':f"{mts.scler.median():.2f} vs {mts.contra.median():.2f}",
               'p':round(pw,4),'n':len(mts)}]).to_csv('supp_fig_s5_burden_values.csv',index=False)
print("saved supp_fig_s5_burden_hippocampus.pdf")
