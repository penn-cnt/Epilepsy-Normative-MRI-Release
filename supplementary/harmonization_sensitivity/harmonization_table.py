import pandas as pd, numpy as np, sys
sys.path.insert(0,'.')
from common import *
# Fornix = F_R only so all three versions are like-for-like (F_L imputed inside v1/v3 ComBat)
bilat=dict(bilateral_tracts); bilat['Fornix']=['F_R']
morph=morph_cohort()
VERS=[('ComBat as run',f'{BASE}/node_zscores.csv'),
      ('No ComBat','node_zscores_unharmonized.csv'),
      ('ComBat + age/sex, ref-only','node_zscores_v3.csv')]
allrows=[]; ns={}
for label,path in VERS:
    z=pd.read_csv(path,usecols=['subject','tract','scalar','z_score','cohort'],low_memory=False)
    z=z[z.cohort=='Clinical']
    wide=build_wide(z,morph,bilat); g=groups(wide); ns[label]={k:len(v) for k,v in g.items()}
    for ga,gb,la,lb in comparisons:
        fa,md=tract_auc_family(g[ga],g[gb],tract_order)
        for res,dom in [(fa,'FA'),(md,'MD')]:
            for i,t in enumerate(tract_order):
                a,lo,hi=fold_auc(res['auc'][i],res['ci_lo'][i],res['ci_hi'][i])
                allrows.append({'Comparison':f'{la} vs. {lb}','Tract':t,'Metric':dom,'Version':label,
                    'AUC (95% CI)':f"{a:.3f} ({lo:.3f}-{hi:.3f})",'q':round(res['q'][i],3),
                    'AUC':round(a,4),'Higher In':la if res['auc'][i]>=0.5 else lb})
    print(label,'done',ns[label])
long=pd.DataFrame(allrows); long.to_csv('harmonization_table_long.csv',index=False)
wide_t=long.pivot_table(index=['Comparison','Tract','Metric'],columns='Version',
                        values=['AUC (95% CI)','q'],aggfunc='first')
wide_t.columns=[f'{b} — {a}' for a,b in wide_t.columns]
order=['Epilepsy vs. Functional Seizures','Focal vs. Generalized','Drug-Resistant vs. Drug-Responsive']
wide_t=wide_t.reset_index()
wide_t['_o']=wide_t.Comparison.map({c:i for i,c in enumerate(order)})
wide_t['_t']=wide_t.Tract.map({t:i for i,t in enumerate(tract_order)})
wide_t=wide_t.sort_values(['_o','_t','Metric']).drop(columns=['_o','_t'])
cols=['Comparison','Tract','Metric']+[f'{v} — {m}' for v in [x[0] for x in VERS] for m in ['AUC (95% CI)','q']]
wide_t=wide_t[[c for c in cols if c in wide_t.columns]]
wide_t.to_csv('harmonization_table.csv',index=False)
pd.DataFrame(ns).to_csv('harmonization_table_group_ns.csv')
print("\nn per group per version:\n",pd.DataFrame(ns).to_string())
print(f"\nsaved harmonization_table.csv ({len(wide_t)} rows)")
for v in [x[0] for x in VERS]:
    s=(long[(long.Version==v)&(long.q<0.05)])
    print(f"  {v}: {len(s)}/54 FDR-significant")
