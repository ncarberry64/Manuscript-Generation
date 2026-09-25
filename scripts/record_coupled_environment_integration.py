"""Record completed evidence; never run or bypass the full test suite."""
import hashlib
import json
from pathlib import Path
import re
import subprocess
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]
BASE='e32b13ad86d1b66889a4a0483fed9ef20ae84561'

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    replay=ROOT/'artifacts/provenance/R1_COUPLED_ENVIRONMENT_INTEGRATION_REPLAY.json'
    a=json.loads(replay.read_text())
    xml=ROOT/'artifacts/coupled_environment_full_pytest.xml'
    collection=ROOT/'artifacts/coupled_environment_collection.txt'
    expected={s.strip() for s in collection.read_text().splitlines() if '::' in s}
    tests={'status':'RUNNING_OR_PENDING','expected_count':len(expected)}
    if xml.exists():
        cases=ET.parse(xml).findall('.//testcase')
        failures=[x.get('name') for x in cases if x.find('failure') is not None or x.find('error') is not None]
        skips=[x.get('name') for x in cases if x.find('skipped') is not None]
        actual={x.get('classname','').replace('.','/')+'.py::'+x.get('name','') for x in cases}
        identities_match=actual==expected
        allowed=['test_dae_preserves_00_0i_constraints_from_index_consistent_seed']
        ok=identities_match and not failures and skips==allowed
        tests.update(status='PASS' if ok else 'FAIL',count=len(cases),passed=len(cases)-len(failures)-len(skips),
                     failures=failures,skipped=skips,collection_identities_match=identities_match,sha256=digest(xml))
    log=ROOT/'manuscript/main.log'
    problems=re.findall(r'^.*(?:Overfull|undefined|Warning|! LaTeX Error).*$' ,log.read_text(errors='replace'),flags=re.M)
    protected=subprocess.run(['git','-C',str(ROOT),'diff','--quiet',BASE,'--','code','preregistration','analyses'],check=False).returncode==0
    ranks=all(r['rank_rel_1e-8']==r['q2_Pi2_rank_rel_1e-8']==2 and r['matter_delta_v_determinant']!=0 and r['q2_Pi2_matter_determinant']!=0 for r in a['rows'])
    invariants=protected and ranks and not a['retuning'] and not a['observational_environment_selection']
    passed=invariants and not problems and tests['status']=='PASS'
    status='PASS' if passed else ('PENDING_FULL_SUITE' if tests['status']=='RUNNING_OR_PENDING' and invariants and not problems else 'FAIL')
    files=['preregistration/prediction_manifest.json','manuscript/main.pdf','manuscript/main.tex',
           'manuscript/sections/04ajb_coupled_environmental_state.tex','scripts/audit_r1_coupled_environmental_state.py',
           'tests/test_coupled_environmental_state.py','artifacts/provenance/R1_COUPLED_ENVIRONMENT_INTEGRATION_REPLAY.json',
           'artifacts/coupled_environment_focused_pytest.xml','artifacts/coupled_environment_collection.txt']
    receipt=dict(base_commit=BASE,branch='codex/r1-coupled-environment-integration',
        statuses=dict(COUPLED_ENVIRONMENTAL_TEMPORAL_RANK=2,MATTER_ONLY_TEMPORAL_RANK=2,
            SPATIAL_SINGLE_PROFILE='CONDITIONAL_RANK_LE_1',ISOLATED_TEMPORAL_ATTRACTOR_REQUIRED=False,
            ANCHOR_INITIAL_STATE_DERIVED_FROM_EARLY_UNIVERSE=False,BHSM_R1_NORMALIZATION_BRIDGE='OPEN',
            RETUNING=False,MANUSCRIPT_INTEGRATION=status),
        full_suite=tests,latex_log_problems=problems,protected_scientific_paths_unchanged=protected,
        visual_review_pages=[2,5,11,12,13,30],remote_main_modified=False,
        hashes={p:digest(ROOT/p) for p in files})
    output=ROOT/'artifacts/provenance/R1_COUPLED_ENVIRONMENT_INTEGRATION_RECEIPT.json'
    output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n',encoding='utf8')
    print(json.dumps(receipt['statuses']))

if __name__=='__main__':main()
