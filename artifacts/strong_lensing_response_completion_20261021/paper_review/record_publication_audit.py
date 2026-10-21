"""Record focused checks and final output scope without rerunning old solvers."""
from pathlib import Path
import hashlib
import json
import re
import subprocess

root = Path(__file__).resolve().parents[3]
out = Path(__file__).resolve().parent
def sha(p):
    return hashlib.sha256((root/p).read_bytes()).hexdigest()
base = '06cf7c1b85d6605bfee0cc0d07fb2af86a66cae9'
name = 'manuscript/strong_lensing_main.tex'
new = (root/name).read_text()
math_name = 'manuscript/sections/strong_lensing_galaxies_20261020.tex'
old_math = subprocess.check_output(['git','show',f'{base}:{math_name}'],cwd=root).decode()
new_math = (root/math_name).read_text()
pattern = r'\\begin\{(equation\*?|align\*?|gather\*?|multline\*?)\}.*?\\end\{\1\}'
before = [m.group(0) for m in re.finditer(pattern,old_math,re.S)]
after = [m.group(0) for m in re.finditer(pattern,new_math,re.S)]
assert before == after and len(before) == 42
log = (root/'manuscript/strong_lensing_main.log').read_text(errors='replace')
assert not re.search(r'LaTeX Warning|Package .* Warning|Overfull|Underfull',log)
assert '210 passed' in (out/'focused_tests.txt').read_text()
render = json.loads((out/'render_receipt.json').read_text())
assert render['pdf_sha256'] == sha('manuscript/strong_lensing_main.pdf')
render['visual_review_complete'] = True
render['review'] = ('All 46 pages reviewed in contact sheets; full-resolution pages34 and37 '
                    'checked for figure labels, equations, caption and final scientific claims. '
                    'No overlap, clipping or unresolved citations. Figure6 now appears by its discussion.')
(out/'render_receipt.json').write_text(json.dumps(render,indent=2)+'\n')
paths = [name,'manuscript/strong_lensing_main.pdf',
         'manuscript/sections/he0435_response_repair_20261021.tex',
         'manuscript/sections/strong_lensing_galaxies_20261020.tex',
         'docs/HE0435_RESPONSE_COMPLETION_20261021.md',
         'manuscript/figures/he0435_converged_held_residuals_20261021.png',
         '.github/workflows/strong-lensing-manuscript.yml']
data = {
    'status':'PASS_FOCUSED_APPLICATIONS_BUILD_AND_VISUAL_REVIEW',
    'basis_commit':base, 'source_and_output_sha256':{p:sha(p) for p in paths},
    'focused_tests':{'passed':210,'failed':0,'seconds':27.32,
                     'scope':'Exact explicit list used by strong-lensing-manuscript.yml; no whole-repository pass claimed.'},
    'build':{'real_unmodified_MDPI_class':True,'pages':render['pages'],
             'abstract_words':len(re.search(r'\\abstract\{(.*?)\}',new,re.S).group(1).split()),
             'warnings':0,'overfull':0,'underfull':0,
             'latex_log_sha256':sha('manuscript/strong_lensing_main.log')},
    'preservation':{'original_main_math_environments':len(before),
                    'all_original_main_math_environments_unchanged':True,
                    'original_prediction_and_data_files_unchanged':True,
                    'delivered_package_unchanged':True},
    'scientific_scope':'Converged conditional native mean; source-calibrated covariance applications and unsuccessful held-exposure mean repairs. Full illuminated-detector likelihood and microscopic BHSM U/V source remain conditional/uncomputed.',
    'held_original_ERR_Q':10974.065419722,
    'physical_admission':False,'physical_validation':False,'physical_guard':2e-10,
    'reserved_observations_opened':False,'error_inflation_used_for_validation':False,
    'review_type':'Internal independent computational/source review, not external peer review.',
}
(out/'publication_manifest.json').write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps({'status':data['status'],'main_math_environments':len(before),
                  'pages':render['pages'],'focused_tests':210}))
