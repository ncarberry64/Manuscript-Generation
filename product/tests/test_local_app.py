import json
import threading
import time
from urllib.request import Request,urlopen
from urllib.error import HTTPError
import pytest
from bhtracker.server import make_server

@pytest.fixture
def service():
    server=make_server(0);thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
    yield f'http://127.0.0.1:{server.server_port}'
    server.shutdown();server.server_close();server.executor.shutdown(wait=True)

def request(base,path,body=None,headers=None):
    data=None if body is None else json.dumps(body).encode()
    req=Request(base+path,data=data,headers={'Content-Type':'application/json',**(headers or {})})
    with urlopen(req,timeout=20) as response:return response.status,response.read(),response.headers

def test_app_full_flow(service):
    for path in ['/','/app.js','/style.css']:
        status,body,headers=request(service,path);assert status==200 and len(body)>100
    _,body,_=request(service,'/api/examples');examples=json.loads(body);assert len(examples)>=3
    _,body,_=request(service,'/api/examples/nearby_motion');dataset=json.loads(body)
    _,body,_=request(service,'/api/validate',{'dataset':dataset});assert json.loads(body)['valid']
    status,body,_=request(service,'/api/fit',{'dataset':dataset,'family':'gr','nuisance':'shear','starts':1})
    assert status==202;job=json.loads(body)['job_id']
    deadline=time.monotonic()+20
    while time.monotonic()<deadline:
        _,body,_=request(service,'/api/jobs/'+job);out=json.loads(body)
        if out['status']!='running':break
        time.sleep(.02)
    assert out['status']=='complete',out
    result=out['result'];assert result['status']=='TRACKED'
    _,body,_=request(service,'/api/report',{'result':result});report=json.loads(body)['html']
    assert '<html' in report and 'TRACKED' in report and 'known' in report
    assert result['data_summary']['provenance']['kind']=='synthetic'

def test_bad_covariance_rejected(service):
    _,body,_=request(service,'/api/examples/nearby_motion');dataset=json.loads(body)
    dataset['covariance'][0][0]=0
    with pytest.raises(HTTPError) as error:request(service,'/api/validate',{'dataset':dataset})
    assert error.value.code==400

def test_cross_origin_and_path_rejected(service):
    with pytest.raises(HTTPError) as error:request(service,'/api/health',headers={'Origin':'https://example.com'})
    assert error.value.code==403
    with pytest.raises(HTTPError) as error:request(service,'/../pyproject.toml')
    assert error.value.code==404
