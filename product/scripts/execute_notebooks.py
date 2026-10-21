from pathlib import Path
import hashlib
import json
import nbformat
from nbclient import NotebookClient
ROOT=Path(__file__).resolve().parents[1]
receipts=[]
for path in sorted((ROOT/'notebooks').glob('*.ipynb')):
    print('Executing '+path.name,flush=True)
    book=nbformat.read(path,as_version=4)
    NotebookClient(book,timeout=180,kernel_name='python3',resources={'metadata':{'path':str(ROOT)}}).execute()
    nbformat.write(book,path)
    receipts.append({'notebook':path.name,'status':'EXECUTED_NO_ERRORS','code_cells':sum(c.cell_type=='code' for c in book.cells),
                     'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
out=ROOT/'validation';out.mkdir(exist_ok=True)
(out/'notebook_receipt.json').write_text(json.dumps(receipts,indent=2)+'\n')
print(json.dumps(receipts,indent=2))
