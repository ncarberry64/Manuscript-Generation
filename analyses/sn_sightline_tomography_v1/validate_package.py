"""Scientific invariants and provenance checks; no model optimization."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
import frozen_analysis as fa
ROOT=Path(__file__).resolve().parent
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
 return h.hexdigest()
def main():
 checks={}; previous={}
 for name in ['BHSM_FULL_COVARIANCE_RESULTS','BHSM_PV_INTERPRETATION_TESTS']:
  path=ROOT.parent/name
  if not (path/'RESULT_MANIFEST.json').exists(): previous[name]={'status':'not locally available; original manifest bundled'}; continue
  manifest=json.loads((path/'RESULT_MANIFEST.json').read_text()); failures=[f for f,h in manifest.items() if not (path/f).exists() or digest(path/f)!=h]
  previous[name]=dict(status='verified unchanged',files=len(manifest),failures=failures); assert not failures
  (ROOT/'audit'/(name+'_ORIGINAL_MANIFEST.json')).write_text(json.dumps(manifest,indent=2))
 inputs=json.loads((ROOT/'INPUT_PROVENANCE.json').read_text()); bad=[x['path'] for x in inputs if digest(ROOT/x['path'])!=x['sha256']]; assert not bad
 checks['input_hashes_verified']=len(inputs); checks['previous_packages']=previous
 p=pd.read_csv(ROOT/'interpretation/PANTHEON_NUISANCE_FEATURES.csv'); assert not any('pv' in c.lower() or 'common_delta' in c.lower() for c in p.columns)
 checks['no_fitted_velocity_or_common_amplitude_nuisance']=True
 analysis=pd.read_csv(ROOT/'interpretation/ANALYSIS_ROWS_WITH_SOURCE_FEATURES.csv.gz'); lc=pd.read_csv(ROOT/'lightcurves/LIGHTCURVE_RECONSTRUCTION.csv').set_index('uid')
 for c in ['delta_c','delta_x1','delta_t0']:
  value=analysis.uid.map(lc[c].where(lc.status=='ok')).to_numpy(); check=analysis['lc_'+c].to_numpy(); assert np.allclose(value,check,equal_nan=True)
 checks['source_features_match_final_reconstruction']=True
 fold=pd.read_csv(ROOT/'interpretation/FOLD_MEMBERSHIP.csv.gz'); ids=analysis.set_index('uid').sn_block
 for name,g in fold.groupby('test'):
  train=set(g.loc[g.role=='train','uid'].map(ids)); test=set(g.loc[g.role=='test','uid'].map(ids)); assert not train&test
 checks['folds_without_same_SN_training_test_overlap']=fold.test.nunique()
 results=pd.read_csv(ROOT/'interpretation/LIKELIHOODS.csv'); assert np.max(abs(results.logL_frozen-results.logL_null+results.delta_chi2/2))<1e-8
 checks['null_frozen_loglikelihood_identity']=True
 mocks=np.load(ROOT/'interpretation/GAUSSIAN_NULL_MOCKS.npz'); z=[]; A=-0.04321623435997208
 for row in results.itertuples():
  vals=mocks[row.test]; sd=2*abs(A)*np.sqrt(row.template_information)
  if sd>1e-10: z.append(abs(vals.mean()-A*A*row.template_information)/(sd/np.sqrt(len(vals))))
 checks['max_mock_mean_standard_error_deviation']=max(z); assert max(z)<5
 state=np.load(ROOT/'interpretation/P_SKY_CROSSFIT_JOINT_OPERATORS.npz'); B=state['B']; C=state['C']; projection=B@C@B.T; error=float(np.max(abs(projection@projection-projection))); checks['crossfit_projected_covariance_idempotence_error']=error; assert error<1e-5
 r=results.set_index('test'); assert abs(r.loc['P_original_heldout_release_only','delta_chi2']-2.984509)<1e-5; assert abs(r.loc['DES_original_heldout_release_only','delta_chi2']-.240298)<1e-5
 checks['previous_frozen_likelihoods_reproduced']=True
 fa.savejson(ROOT/'VALIDATION.json',checks); print(json.dumps(checks,indent=2))
if __name__=='__main__': main()
