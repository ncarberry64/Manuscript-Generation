"""Preserve both published alternatives for the one ambiguous nominal row."""
from pathlib import Path
import pandas as pd
import lightcurve_stage as lc
def main():
 tasks=lc.build_tasks(); row,item=next(t for t in tasks if t[0]['uid']=='P0015')
 idx=pd.read_csv(lc.OUT/'RAW_LIGHTCURVE_INDEX.csv',keep_default_na=False)
 candidates=idx.loc[(idx.group=='PS1s_LOWZ_JRK07_NOCFA_DS15')&(idx.SNID.map(lc.norm).eq(lc.norm(row['CID']))|idx.IAUC.map(lc.norm).eq(lc.norm(row['CID'])))]
 output=[]
 for number,candidate in enumerate(candidates.to_dict('records'),1):
  r=row.copy(); r['uid']=row['uid']+f'_ALT{number}'; b,j,s=lc.reconstruct((r,candidate)); output.append(b)
  pd.DataFrame(j).to_csv(lc.OUT/(r['uid']+'_JACKKNIFE.csv'),index=False); pd.DataFrame(s).to_csv(lc.OUT/(r['uid']+'_SYSTEMATICS.csv'),index=False)
 pd.DataFrame(output).to_csv(lc.OUT/'AMBIGUOUS_SOURCE_ALTERNATIVES.csv',index=False)
 print(pd.DataFrame(output)[['uid','CID','status','source']].to_string(index=False))
if __name__=='__main__': main()
