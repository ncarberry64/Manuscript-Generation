"""Loopback-only app. Measurements stay in memory; no remote service used."""
from concurrent.futures import ThreadPoolExecutor
from http.server import BaseHTTPRequestHandler,ThreadingHTTPServer
from pathlib import Path
import json
import threading
import uuid
import webbrowser
from .data import load_dataset
from .examples import example_names,make_example
from .inference import fit,profile
from .report import report_html

WEB=Path(__file__).resolve().parent/'web'

def make_server(port=8765):
    pool=ThreadPoolExecutor(max_workers=1)
    jobs={};lock=threading.Lock()
    class Handler(BaseHTTPRequestHandler):
        def log_message(self,*args):pass
        def send_json(self,value,status=200):
            blob=json.dumps(value,allow_nan=False).encode()
            self.send_response(status);self.send_header('Content-Type','application/json; charset=utf-8')
            self.send_header('Content-Length',str(len(blob)));self.send_header('Cache-Control','no-store')
            self.send_header('X-Content-Type-Options','nosniff');self.end_headers();self.wfile.write(blob)
        def local_request(self):
            hosts={f'127.0.0.1:{self.server.server_port}',f'localhost:{self.server.server_port}'}
            if self.headers.get('Host') not in hosts:return False
            origin=self.headers.get('Origin')
            return origin is None or origin in {f'http://{h}' for h in hosts}
        def do_GET(self):
            if not self.local_request():self.send_json({'error':'Local origin required'},403);return
            if self.path=='/api/health':self.send_json({'status':'ready','product':'BHSM Blackhole Tracker','version':'0.1.1'});return
            if self.path=='/api/examples':self.send_json(example_names());return
            if self.path.startswith('/api/examples/'):
                try:self.send_json(make_example(self.path.rsplit('/',1)[1])[0])
                except ValueError as error:self.send_json({'error':str(error)},404)
                return
            if self.path.startswith('/api/jobs/'):
                name=self.path.rsplit('/',1)[1]
                with lock:item=jobs.get(name)
                self.send_json(item or {'error':'Unknown job'},200 if item else 404);return
            path={'/':'index.html','/index.html':'index.html','/app.js':'app.js','/style.css':'style.css'}.get(self.path)
            if not path:self.send_json({'error':'Not found'},404);return
            blob=(WEB/path).read_bytes();self.send_response(200)
            self.send_header('Content-Type',{'html':'text/html; charset=utf-8','js':'text/javascript','css':'text/css'}[path.rsplit('.',1)[1]])
            self.send_header('Content-Length',str(len(blob)));self.send_header('Cache-Control','no-store')
            self.send_header('Content-Security-Policy',"default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; connect-src 'self'; img-src 'self' data:; frame-ancestors 'none'")
            self.end_headers();self.wfile.write(blob)
        def do_POST(self):
            if not self.local_request():self.send_json({'error':'Local origin required'},403);return
            try:
                size=int(self.headers.get('Content-Length','0'))
                if not 0<size<=10_000_000:raise ValueError('Request must contain at most 10 MB of JSON')
                data=json.loads(self.rfile.read(size))
                if self.path=='/api/validate':
                    dataset=load_dataset(data['dataset']);self.send_json({'valid':True,'summary':dataset.summary()});return
                if self.path=='/api/report':self.send_json({'html':report_html(data['result'])});return
                if self.path=='/api/fit':
                    dataset=load_dataset(data['dataset'])
                    with lock:
                        if len(jobs)>=64:raise ValueError('Session limit reached; restart the local app to clear in-memory results')
                        name=uuid.uuid4().hex;jobs[name]={'status':'running'}
                    def work():
                        try:
                            result=fit(dataset,family=data.get('family','gr'),nuisance=data.get('nuisance','shear'),starts=data.get('starts',2))
                            item={'status':'complete','result':result}
                        except Exception as error:item={'status':'failed','error':str(error)}
                        with lock:jobs[name]=item
                    pool.submit(work);self.send_json({'job_id':name},202);return
                self.send_json({'error':'Not found'},404)
            except (ValueError,KeyError,TypeError) as error:self.send_json({'error':str(error)},400)
    server=ThreadingHTTPServer(('127.0.0.1',port),Handler)
    server.executor=pool
    return server

def serve(port=8765,open_browser=False):
    server=make_server(port)
    url=f'http://127.0.0.1:{server.server_port}/'
    print('Blackhole Tracker ready at '+url,flush=True)
    print('Local-only research app. Ctrl+C stops the server; inputs/results stay in this process.',flush=True)
    if open_browser:webbrowser.open(url)
    try:server.serve_forever()
    except KeyboardInterrupt:pass
    finally:server.server_close();server.executor.shutdown(wait=False,cancel_futures=True)
