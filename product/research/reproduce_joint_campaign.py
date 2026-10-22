#!/usr/bin/env python
"""Portable reproduction of the internally frozen exploratory2023 joint campaign.

python product/research/reproduce_joint_campaign.py --sourcepack SOURCE.zip --output NEWDIR
Use --verify-only to verify/import the exact observations without refitting.
"""
import argparse,datetime,hashlib,json,os,platform,shutil,stat,subprocess,sys,zipfile
from pathlib import Path,PurePosixPath,PureWindowsPath

SOURCEPACK_SHA256='71d36a03e4ba275e91bea62f4784a5193f34ad32d3777b637453f91151bafa99'
REFERENCE_INPUT_SHA256='d4fc2bfc9320c7eab8a33b93517f68f035a5d498ffb35ff548a638648143af96'
SCIENCE_SHA256='6eaf050481fc23aebf4b929e49fdab57e324ee80eb08e1547b69f056269eeaaa'


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write_same_or_new(path,data):
    path=Path(path)
    if path.exists():
        if path.read_bytes()!=data:raise ValueError(f'Refusing to replace different existing bytes: {path}')
    else:path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
def write_json(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def source_directory(source,out):
    if source.is_dir():return source
    if digest(source)!=SOURCEPACK_SHA256:raise ValueError('Sourcepack ZIP hash differs from the frozen artifact')
    target=(out/'sourcepack').resolve();target.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(source) as archive:
        for item in archive.infolist():
            name=item.filename;part=PurePosixPath(name)
            if part.is_absolute() or PureWindowsPath(name).drive or '..' in part.parts or '\\' in name or stat.S_ISLNK(item.external_attr>>16):raise ValueError('Unsafe sourcepack archive entry')
            destination=(target/Path(*part.parts)).resolve()
            if not destination.is_relative_to(target):raise ValueError('Sourcepack entry escapes requested output')
            if not item.is_dir():write_same_or_new(destination,archive.read(item))
    return target/'portable_2023_benchmark_v1'


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--sourcepack',required=True,type=Path);parser.add_argument('--output',required=True,type=Path);parser.add_argument('--repository',type=Path,default=Path(__file__).resolve().parents[2]);parser.add_argument('--verify-only',action='store_true');parser.add_argument('--stages',default='baselines,profiles,predictions,controls,summary',help='Comma-separated stages; import/source verification always runs')
    a=parser.parse_args();repo=a.repository.resolve();out=a.output.resolve();source=a.sourcepack.resolve();out.mkdir(parents=True,exist_ok=True);sys.path.insert(0,str(repo/'product'))
    from bhtracker.joint_mass_v2 import JointData
    from bhtracker.joint_benchmark_v2 import import_2023_benchmark
    from bhtracker.provenance import project_citation,VERSION
    if digest(repo/'product/bhtracker/joint_mass_v2.py')!=SCIENCE_SHA256:raise ValueError('Scientific source differs; use the frozen version or declare a new analysis version')
    resource=Path(__file__).resolve().parent/'joint_campaign_drivers';drivers=json.loads((resource/'ORIGINAL_DRIVER_MANIFEST.json').read_text())
    for name,record in drivers['files'].items():
        if digest(resource/name)!=record['sha256']:raise ValueError(f'Original driver hash mismatch: {name}')
        write_same_or_new(out/name,(resource/name).read_bytes())
    write_same_or_new(out/'ORIGINAL_DRIVER_MANIFEST.json',(resource/'ORIGINAL_DRIVER_MANIFEST.json').read_bytes())
    raw=source_directory(source,out);data=import_2023_benchmark(raw)
    if data.sha256!=REFERENCE_INPUT_SHA256:raise ValueError('Derived input differs from frozen study. Use the recorded Python/NumPy/SciPy/Astropy environment and inspect derived clock metadata; do not silently reuse reference results.')
    inputfile=out/'JOINT_2023_NATIVE_INPUT_v1.json';write_same_or_new(inputfile,(json.dumps(data.document,separators=(',',':'))+'\n').encode())
    environment={'python':platform.python_version(),**{name:__import__(name).__version__ for name in ['numpy','scipy','astropy']}}
    receiptfile=out/'PORTABLE_REPRODUCTION_MANIFEST.json'
    if receiptfile.exists():
        previous=json.loads(receiptfile.read_text())
        if previous['environment']!=environment:raise ValueError('Resume environment differs; use a new output directory to retain prior results')
    receipt={'schema':'bhtracker.portable_joint_reproduction.v1','sourcepack_sha256':SOURCEPACK_SHA256 if source.is_file() else None,'source_directory':str(raw),'scientific_source_sha256':SCIENCE_SHA256,'input_sha256':data.sha256,'hash_algorithm':'bhtracker.canonical_numeric_json.v2','project_citation':project_citation(),'product_version':VERSION,'environment':environment,'driver_manifest':drivers,'commands':[],'no_automatic_mass_or_BH_acceptance':True}
    if receiptfile.exists():receipt['commands']=json.loads(receiptfile.read_text()).get('commands',[])
    write_json(receiptfile,receipt)
    if a.verify_only:print(json.dumps({'verified':True,'photometric_rows':len(data.document['photometry']),'astrometric_rows':len(data.document['astrometry']),'input_sha256':data.sha256}));return
    stages=set(a.stages.split(','));env=dict(os.environ);env['PYTHONPATH']=str(repo/'product')+os.pathsep+env.get('PYTHONPATH','')
    def run(label,args,expected=None,complete=None):
        if expected is not None and expected.exists():
            existing=json.loads(expected.read_text())
            if existing.get('scientific_source_sha256',SCIENCE_SHA256)!=SCIENCE_SHA256:raise ValueError('Existing result has different scientific source')
            if existing.get('input_sha256',data.sha256)!=data.sha256 or existing.get('full_input_sha256',data.sha256)!=data.sha256:raise ValueError('Existing result has different native joint input')
            if 'base_result' in existing and existing['base_result']!=baseline:raise ValueError('Existing control summary refers to a different baseline')
            if 'reference_objective_deviance' in existing and abs(existing['reference_objective_deviance']-json.loads((out/baseline).read_text())['objective_deviance'])>1e-8:raise ValueError('Existing profile has a different reference minimum')
            if complete is None or complete(existing):return
            snapshot=out/'preserved_resume_snapshots'/datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');snapshot.mkdir(parents=True)
            for f in out.glob(expected.stem.replace('_PARAMETRIC_CONTROLS20','')+'*'):
                if f.is_file():shutil.copy2(f,snapshot/f.name)
        record={'label':label,'argv':[str(x) for x in args],'cwd':str(repo),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};receipt['commands'].append(record);write_json(receiptfile,receipt)
        print(label,flush=True)
        with (out/(label+'.log')).open('w',encoding='utf-8') as log:completed=subprocess.run([str(x) for x in args],cwd=repo,env=env,stdout=log,stderr=subprocess.STDOUT)
        record['exit_code']=completed.returncode;record['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();write_json(receiptfile,receipt)
        if completed.returncode:raise RuntimeError(f'{label} failed; inspect its retained log')
    candidates=[('G0','JOINT_V2_G0_DIAGONAL16',[]),('G1','JOINT_V2_G1_MARGINAL16',['--gp']),('G2','JOINT_V2_G2_CALIBRATED16',['--gp','--options','{"noise_model":"calibrated_marginal_ou"}'])]
    if 'baselines' in stages:
        for code,stem,options in candidates:
            run(stem,[sys.executable,'-m','bhtracker.joint_workflow_v2','fit',inputfile,out/(stem+'.json'),*options],out/(stem+'.json'))
            run(stem+'_CONTINUATION',[sys.executable,out/'continue_v2.py',stem+'.json'],out/(stem+'_CONTINUED.json'))
    for code,stem,options in candidates[1:]:
        baseline=stem+'_CONTINUED.json';profile='JOINT_V2_'+code+'_MASS_PROFILE31.json'
        if 'profiles' in stages:
            run(code+'_PROFILE31',[sys.executable,'-m','bhtracker.joint_workflow_v2','profile',inputfile,out/baseline,out/profile],out/profile)
            refined='JOINT_V2_'+code+'_MASS_REFINED_SUPPORT.json';run(code+'_PROFILE_REFINEMENT',[sys.executable,out/'refine_v2_profiles.py',baseline,profile,refined],out/refined)
        if 'predictions' in stages:
            run(code+'_ORIGINAL_UTC_CONDITIONAL_DIAGNOSTIC',[sys.executable,out/'run_v2_blocked_prediction.py',baseline],out/(Path(baseline).stem+'_BLOCKED_PREDICTION.json'))
            run(code+'_V2_1_NIGHT_PREDICTION',[sys.executable,out/'run_v2_1_blocked_prediction.py',baseline],out/(Path(baseline).stem+'_V2_1_BLOCKED_PREDICTION.json'))
        if 'controls' in stages:
            file=out/(Path(baseline).stem+'_PARAMETRIC_CONTROLS20.json');run(code+'_PARAMETRIC20',[sys.executable,out/'run_v2_parametric_controls.py',baseline],file,lambda r:r.get('control_count')==20)
    if 'summary' in stages:
        table=[]
        for code,stem,options in candidates:
            result=json.loads((out/(stem+'_CONTINUED.json')).read_text());table.append({'candidate':code,'solver_status':result['status'],'conditional_mass_point_msun':result['parameters']['mass_msun'],'accepted_mass_measurement':False,'failed_structure_bins':result['residual_structure']['failed_bins'],'active_search_bounds':result['nonlinear_active_bounds'],'numerical_convergence':result['numerical_convergence']})
        write_json(out/'PORTABLE_CAMPAIGN_SUMMARY.json',{'schema':'bhtracker.portable_campaign_summary.v1','candidates':table,'interpretation':'No automatic physical acceptance. Apply frozen profile, noise, predictive, timing/calibration and darkness gates; retain every inadequate result.','all_source_rows_preserved':15533,'new_discoveries_claimed':[]})

if __name__=='__main__':main()
