import pandas as pd, numpy as np
from common import *

# Fornix = F_R only, so harmonized and unharmonized are like-for-like
# (F_L is mean-imputed inside ComBat for 68 subjects; absent pre-ComBat)
bilat_noFL=dict(bilateral_tracts); bilat_noFL['Fornix']=['F_R']
morph=morph_cohort()

def run(zdf,label):
    wide=build_wide(zdf,morph,bilat_noFL); g=groups(wide)
    print(f"{label} n per group:", {k:len(v) for k,v in g.items()})
    rows=[]
    for ga,gb,la,lb in comparisons:
        fa,md=tract_auc_family(g[ga],g[gb],tract_order)
        for res,dom in [(fa,'FA'),(md,'MD')]:
            for i,t in enumerate(tract_order):
                a,lo,hi=fold_auc(res['auc'][i],res['ci_lo'][i],res['ci_hi'][i])
                rows.append({'Comparison':f'{la} vs. {lb}','Metric':dom,'Tract':t,
                    'AUC':a,'CI':f"({lo:.3f}-{hi:.3f})",'q':res['q'][i],
                    'Higher In':la if res['auc'][i]>=0.5 else lb})
    return pd.DataFrame(rows), {k:len(v) for k,v in g.items()}

print("=== harmonized ===")
h=pd.read_csv(f'{BASE}/node_zscores.csv',usecols=['subject','tract','scalar','z_score','cohort'])
h=h[h.cohort=='Clinical']
H,nh=run(h,'harmonized')
print("=== unharmonized ===")
u=pd.read_csv('node_zscores_unharmonized.csv',usecols=['subject','tract','scalar','z_score','cohort'])
u=u[u.cohort=='Clinical']
U,nu=run(u,'unharmonized')

m=H.merge(U,on=['Comparison','Metric','Tract'],suffixes=('_harm','_unharm'))
m['dAUC']=m.AUC_unharm-m.AUC_harm
m['sig_harm']=m.q_harm<0.05; m['sig_unharm']=m.q_unharm<0.05
m['FLAG']=np.where(m.dAUC.abs()>0.03,'AUC>0.03','')
m.loc[m.sig_harm!=m.sig_unharm,'FLAG']=(m.loc[m.sig_harm!=m.sig_unharm,'FLAG']+' SIG-CHANGE').str.strip()
m['Harmonized AUC (95% CI)']=m.apply(lambda r:f"{r.AUC_harm:.3f} {r.CI_harm}",axis=1)
m['Unharmonized AUC (95% CI)']=m.apply(lambda r:f"{r.AUC_unharm:.3f} {r.CI_unharm}",axis=1)
out=m[['Comparison','Tract','Metric','Harmonized AUC (95% CI)','q_harm',
       'Unharmonized AUC (95% CI)','q_unharm','dAUC','Higher In_harm','Higher In_unharm','FLAG']]
out=out.rename(columns={'q_harm':'q (harm)','q_unharm':'q (unharm)','dAUC':'delta AUC',
                        'Higher In_harm':'Higher In (harm)','Higher In_unharm':'Higher In (unharm)'})
out.round(3).to_csv('step3_harmonization_sensitivity.csv',index=False)
pd.Series(nh).to_csv('step4_n_harmonized.csv'); pd.Series(nu).to_csv('step4_n_unharmonized.csv')
print(f"\n=== rows: {len(out)} | flagged: {(out.FLAG!='').sum()} ===")
print(f"max |dAUC| = {m.dAUC.abs().max():.3f} | mean |dAUC| = {m.dAUC.abs().mean():.3f}")
print(f"sig harmonized: {m.sig_harm.sum()} | sig unharmonized: {m.sig_unharm.sum()}")
if (out.FLAG!='').any(): print("\nFLAGGED:\n",out[out.FLAG!=''].round(3).to_string(index=False))
