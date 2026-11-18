import pandas as pd, numpy as np
from common import *
morph=morph_cohort()
print(f"cohort n = {len(morph)}")
raw=pd.read_csv(f'{BASE}/node_zscores.csv')
clin=raw[raw['cohort']=='Clinical'].copy()
wide=build_wide(clin,morph)
g=groups(wide)
print("n per group:", {k:len(v) for k,v in g.items()})
rows=[]
for ga,gb,la,lb in comparisons:
    fa,md=tract_auc_family(g[ga],g[gb],tract_order)
    for res,dom in [(fa,'Tractometry (FA)'),(md,'Tractometry (MD)')]:
        for i,t in enumerate(tract_order):
            a,lo,hi=fold_auc(res['auc'][i],res['ci_lo'][i],res['ci_hi'][i])
            rows.append({'Comparison':f'{la} vs. {lb}','Domain':dom,'Region/Tract':t,
                'AUC (95% CI)':f"{a:.3f} ({lo:.3f}-{hi:.3f})",
                'Higher In':la if res['auc'][i]>=0.5 else lb,'q':round(res['q'][i],3)})
out=pd.DataFrame(rows); out.to_csv('step1_reproduced.csv',index=False)
tgt=pd.read_csv(f'{BASE}/table_auc_results.csv')
tgt=tgt[tgt.Domain.str.startswith('Tractometry')]
m=out.merge(tgt,on=['Comparison','Domain','Region/Tract'],suffixes=('_new','_tgt'))
m['match_auc']=m['AUC (95% CI)_new']==m['AUC (95% CI)_tgt']
m['match_q']=np.isclose(m.q_new,m.q_tgt,atol=0.0011)
m['match_dir']=m['Higher In_new']==m['Higher In_tgt']
print(f"\nrows compared: {len(m)}")
print(f"  AUC+CI string exact: {m.match_auc.sum()}/{len(m)}")
print(f"  q match (±0.001):    {m.match_q.sum()}/{len(m)}")
print(f"  direction match:     {m.match_dir.sum()}/{len(m)}")
bad=m[~(m.match_auc&m.match_q&m.match_dir)]
if len(bad): print("\nMISMATCHES:\n",bad[['Comparison','Domain','Region/Tract','AUC (95% CI)_new','AUC (95% CI)_tgt','q_new','q_tgt']].to_string(index=False))
else: print("\n*** EXACT REPRODUCTION ***")
