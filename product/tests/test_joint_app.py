"""Exercise native joint input, saved-state verification and portable exports."""
import base64
import copy
import hashlib
from io import BytesIO
import json
import time
import subprocess
import sys
from pathlib import Path
import zipfile
from urllib.error import HTTPError

import numpy as np
import pytest

from test_local_app import service, request
from test_joint_mass import synthetic_document


@pytest.fixture(scope='module')
def joint_result():
    from bhtracker.joint_mass import JointModel, trajectory, fit_joint
    d=synthetic_document()
    d['metadata']['synthetic']=True
    d['target']='SYNTHETIC-JOINT-APP-CONTROL'
    raw=b'Independent synthetic regression control; not a measured event.\n'
    d['native_files']={'control.txt':{'sha256':hashlib.sha256(raw).hexdigest(),'raw_bytes_base64':base64.b64encode(raw).decode()}}
    model=JointModel(d)
    x=np.array([5760.,np.log(240.),-.08,1.7,.12,20.,.05])
    q=trajectory(model.tp,model.decode(x),model.ep)
    u2=np.sum(q*q,axis=1)
    light=20.-2.5*np.log10(1+.05*((u2+2)/np.sqrt(u2*(u2+4))-1))
    qa=trajectory(model.ta,model.decode(x),model.ea)
    position=np.array([3.,-4.])+((model.ta-model.ref)/365.25)[:,None]*np.array([-2.,3.])-.13*(model.ea-model.eref)+4.8*qa/(np.sum(qa*qa,axis=1)+2)[:,None]
    for row,value in zip(d['photometry'],light):row['magnitude']=float(value)
    for row,value in zip(d['astrometry'],position):row.update(north_mas=float(value[0]),east_mas=float(value[1]))
    fitted=fit_joint(d,initial_states=[x+np.array([.1,.001,.001,.001,.001,.001,.0001])],max_nfev=120)
    assert fitted['objective_deviance']<1e-14
    return fitted


def test_joint_saved_result_native_export_and_generic_replay(service,joint_result,tmp_path):
    for path in ['/joint','/joint/','/joint.js']:
        status,body,_=request(service,path)
        assert status==200 and len(body)>100
    for path in ['/','/photometry']:
        _,body,_=request(service,path)
        assert b'href="/joint"' in body
    dataset=joint_result['dataset_input']
    _,body,_=request(service,'/api/joint/validate',{'dataset':dataset})
    assert json.loads(body)['summary']['photometry_rows']==150
    _,body,_=request(service,'/api/joint/result',{'result':joint_result})
    checked=json.loads(body)['replay']
    assert checked['passed'] and checked['same_physical_parameters']
    _,body,_=request(service,'/api/report',{'result':joint_result})
    report=json.loads(body)['html']
    assert 'BHSM' in report and 'Conditional solver parameters' in report
    assert 'does not certify a black hole' in report
    _,body,headers=request(service,'/api/export',{'result':joint_result})
    assert headers['Content-Type']=='application/zip'
    with zipfile.ZipFile(BytesIO(body)) as archive:
        assert json.loads(archive.read('dataset.json'))==dataset
        assert len(archive.read('photometry.csv').splitlines())==151
        assert len(archive.read('astrometry.csv').splitlines())==25
        assert archive.read('native_source_01.bin')==b'Independent synthetic regression control; not a measured event.\n'
        assert json.loads(archive.read('manifest.json'))['replay_method']=='FORWARD_STATE_REPLAY'
        provenance=json.loads(archive.read('software_provenance.json'))
        assert provenance['export_environment']['project_citation']['repository_url']=='https://github.com/ncarberry64/Berger-Hopf-Standard-Model'
        interpretation=json.loads(archive.read('bhsm_interpretation_sources.json'))
        assert interpretation['framework']=='BHSM' and interpretation['source_provenance']['sources']
    bundle=tmp_path/'joint.zip';bundle.write_bytes(body)
    from bhtracker.report import replay_analysis
    replay=replay_analysis(bundle,tmp_path/'replayed')
    assert replay['passed'] and replay['same_native_residuals']
    assert replay['max_astrometric_difference_mas']<1e-12


def test_archived_v1_lossless_browser_json_preserves_whole_float_hash(service,joint_result):
    original=json.dumps(joint_result,allow_nan=False)
    # JavaScript JSON.parse/stringify drops the .0 on a whole-valued JD.
    # Original v1 scientific hashing is preserved; transport must retain text.
    lossy=copy.deepcopy(joint_result)
    lossy['dataset_input']['reference_jd_tdb']=int(lossy['dataset_input']['reference_jd_tdb'])
    with pytest.raises(HTTPError) as error:
        request(service,'/api/joint/result',{'result':lossy})
    assert 'hash mismatch' in json.loads(error.value.read())['error']
    _,body,_=request(service,'/api/joint/result',{'result_json':original})
    checked=json.loads(body)
    assert checked['replay']['passed']
    assert json.loads(checked['dataset_json'])==joint_result['dataset_input']
    _,body,_=request(service,'/api/report',{'result_json':original})
    assert 'BHSM' in json.loads(body)['html']


def test_joint_sourcepack_rejects_unreviewed_bytes_before_import(service):
    with pytest.raises(HTTPError) as error:
        request(service,'/api/joint/import-benchmark',{'content_base64':base64.b64encode(b'other data').decode()})
    assert error.value.code==400
    assert 'source pack hash differs' in json.loads(error.value.read())['error']


@pytest.mark.parametrize('change',['mass','residual','native_bytes','scientific_source'])
def test_joint_verification_rejects_tampered_scientific_summary(service,joint_result,change):
    altered=copy.deepcopy(joint_result)
    if change=='mass':altered['parameters']['mass_msun']*=2
    elif change=='residual':altered['astrometric_residuals_NE_mas'][0][0]+=.1
    elif change=='native_bytes':altered['dataset_input']['native_files']['control.txt']['raw_bytes_base64']=base64.b64encode(b'altered').decode()
    else:altered['scientific_source_sha256']='0'*64
    with pytest.raises(HTTPError) as error:
        request(service,'/api/joint/result',{'result':altered})
    assert error.value.code==400


def test_joint_http_fit_and_invalid_controls(service,joint_result):
    data={'dataset':joint_result['dataset_input'],'scientific_version':'v1','options':{'starts':2,'max_nfev':120}}
    _,body,_=request(service,'/api/joint/fit',data)
    job=json.loads(body)['job_id']
    deadline=time.monotonic()+50
    while time.monotonic()<deadline:
        _,body,_=request(service,'/api/jobs/'+job);output=json.loads(body)
        if output['status'] not in ('pending','running'):break
        time.sleep(.03)
    assert output['status']=='complete',output
    fitted=json.loads(output['result_json'])
    assert fitted['dataset_input']==joint_result['dataset_input']
    _,body,_=request(service,'/api/joint/replay',{'result_json':output['result_json']})
    assert json.loads(body)['replay']['passed']
    for changes in [{'scientific_version':'other'},{'options':{'starts':True}}, {'options':{'gp':'yes'}}, {'options':{'invented_covariance':True}}]:
        with pytest.raises(HTTPError) as error:request(service,'/api/joint/fit',{**data,**changes})
        assert error.value.code==400


def test_joint_export_rejects_parameter_tamper_even_with_rehashed_manifest(joint_result,tmp_path):
    from bhtracker.report import export_result,replay_analysis
    export_result(joint_result,tmp_path/'original')
    root=tmp_path/'original'
    result=json.loads((root/'result.json').read_text(encoding='utf-8'))
    result['parameters']['pi_E']*=1.2
    (root/'result.json').write_text(json.dumps(result),encoding='utf-8')
    manifest=json.loads((root/'manifest.json').read_text(encoding='utf-8'))
    manifest['files']['result.json']=hashlib.sha256((root/'result.json').read_bytes()).hexdigest()
    (root/'manifest.json').write_text(json.dumps(manifest),encoding='utf-8')
    with pytest.raises(ValueError,match='physical parameters'):
        replay_analysis(root,tmp_path/'replayed')


def test_python_module_cli_propagates_failed_joint_replay_status(joint_result,tmp_path):
    altered=copy.deepcopy(joint_result)
    altered['photometric_predictions_magnitude'][0]+=.01
    saved=tmp_path/'bad-forward-state.json'
    saved.write_text(json.dumps(altered),encoding='utf-8')
    checked=subprocess.run([sys.executable,'-m','bhtracker','joint','--scientific-version','v1','replay',str(saved)],
                           cwd=Path(__file__).resolve().parents[1],capture_output=True,text=True,timeout=25)
    assert checked.returncode==1,checked.stdout+checked.stderr
    assert json.loads(checked.stdout)['passed'] is False
