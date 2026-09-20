"""Prepare author-review delivery and compact archive; exclude third-party raw data."""
from pathlib import Path
import csv, hashlib, json, shutil, subprocess, zipfile

ROOT=Path(__file__).resolve().parents[1]
OUT=Path.home()/'Downloads/UNIVERSE_COSMOLOGY_SUBMISSION_FINAL'
def copy(src,dst):
    dst.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(src,dst)
def hashes(folder):
    return {p.relative_to(folder).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(folder.rglob('*')) if p.is_file() and p.name!='SHA256.json'}
def archive(folder,dest):
    (folder/'SHA256.json').write_text(json.dumps(hashes(folder),indent=2)+'\n')
    with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted(folder.rglob('*')):
            if p.is_file(): z.write(p,p.relative_to(folder).as_posix())
    with zipfile.ZipFile(dest) as z: assert z.testzip() is None

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    copy(ROOT/'manuscript/main.pdf',OUT/'UNIVERSE_COSMOLOGY_FINAL.pdf')
    copy(ROOT/'submission/universe-source.zip',OUT/'universe-source.zip')
    copy(ROOT/'manuscript/references.bib',OUT/'references.bib')
    for pattern in ['submission/*.md','docs/UNIVERSE_CLAIM_PROVENANCE.md',
                    'docs/UNIVERSE_NUMERICAL_CLAIM_INDEX.md','docs/REVIEWER_RISK_AUDIT.md',
                    'docs/UNIVERSE_STATISTICAL_LANGUAGE_AUDIT.md','docs/UNIVERSE_BIBLIOGRAPHY_VERIFICATION.md']:
        for p in ROOT.glob(pattern): copy(p,OUT/p.name)
    for p in (ROOT/'manuscript/figures').iterdir():
        if p.is_file(): copy(p,OUT/'figures'/p.name)
    for p in [ROOT/'scripts/universe_publication_figures.py', ROOT/'code/figures.py']:
        copy(p,OUT/'figure_scripts'/p.name)
    bundle=OUT/'UNIVERSE_SUBMISSION_REPRODUCIBILITY_PACKAGE'
    bundle.mkdir(exist_ok=True)
    with zipfile.ZipFile(ROOT/'submission/universe-source.zip') as z:
        for name in z.namelist():
            dest=bundle/'manuscript'/name
            dest.resolve().relative_to(bundle.resolve())
            dest.parent.mkdir(parents=True,exist_ok=True); dest.write_bytes(z.read(name))
    for pattern in ['code/*.py','tests/*.py','preregistration/*','artifacts/*.json',
                    'artifacts/*pytest.xml','requirements*.txt','scripts/*.py',
                    'docs/*.md','submission/*.md']:
        for p in ROOT.glob(pattern):
            if p.is_file(): copy(p,bundle/p.relative_to(ROOT))
    base=ROOT/'analyses/sn_sightline_tomography_v1'
    # Author-generated scripts, provenance and aggregate comparison tables only.
    for pattern in ['*.py','*.json','*.md','requirements.txt','interpretation/*.json',
                    'interpretation/LIKELIHOODS.csv','interpretation/*TOMOGRAPHY_SCORES.csv',
                    'interpretation/TOMOGRAPHY_ALL_SUBSETS.csv','inputs/**/*PROVENANCE*.json']:
        for p in base.glob(pattern):
            if p.is_file(): copy(p,bundle/p.relative_to(ROOT))
    for p in [ROOT/'manuscript/build/final-build.log',ROOT/'manuscript/build/final-second-pass.log']:
        copy(p,bundle/'validation'/p.name)
    metadata={'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
              'status':'Prepared for author confirmation; not deposited or submitted',
              'raw_data':'Excluded; original local full package remains unchanged',
              'license':'Author copyright retained; no new redistribution license granted'}
    for dest in [OUT/'RELEASE_METADATA.json',bundle/'RELEASE_METADATA.json']:
        dest.write_text(json.dumps(metadata,indent=2)+'\n')
    readme='''# Universe final reproducibility package

Read submission/SUBMISSION_READINESS.md and docs/UNIVERSE_CLAIM_PROVENANCE.md.
This compact archive contains manuscript sources, figures, author code, tests,
frozen manifests, aggregate derived tables and acquisition metadata. Third-party
raw observations, calibration models and field grids are excluded. Copyright in
author material is retained; no new license is granted by preparation of this bundle.

## Reproduce the saved-data publication figures

From this package root, install requirements-validation.txt, then run:
`python scripts/universe_publication_figures.py`
This reads the included frozen likelihood table and does no fit. Existing R1
figures can be reproduced with `python code/figures.py` (reference computation).

## Build

From manuscript: `latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex`.
Run a second `pdflatex -interaction=nonstopmode -halt-on-error main.tex` pass.
The source-only ZIP is independently self-contained. Python 3.14 on Windows and
MiKTeX were used; dependencies and the earlier complete-suite evidence are retained.

## Inexpensive checks

`python -m pytest -q tests/test_bhsm_make_or_break_prediction.py tests/test_bhsm_optical_low_rank_covariance.py tests/test_bhsm_optical_transfer.py -p no:cacheprovider`

The final-document verifier expects the complete repository sightline snapshot and
PDF text extraction, so its hash check is not runnable on this intentionally reduced
archive alone. Its saved report is in artifacts. This is not a complete observational
rerun package; obtain the original data and frozen operators for that purpose.

## Data acquisition and limitations

Pantheon+: https://github.com/PantheonPlusSH0ES/DataRelease
Pinned commit c447f0fea703fcd0fff57de5000947b5ca81286b.
DES: https://github.com/des-science/DES-SN5YR
Pinned commit c9a4fcafc4cbd19bd750dee47fc76194a45c181f.
Independent fields: https://github.com/rlilow/2MRS-NeuralNet
2MRS catalog: https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJS/199/26

Use analyses/sn_sightline_tomography_v1/INPUT_PROVENANCE.json for exact per-file
URLs and SHA256 values. Preserve original release terms; public access alone is
not a redistribution license. The preparation script also depends on frozen
BHSM_PV_INTERPRETATION_TESTS inputs/operators and earlier package files; these
must be recovered from the original local package, not replaced with current data.
The complete local archive is C:/Users/carbe/Downloads/BHSM_SN_SIGHTLINE_TOMOGRAPHY_V1/
BHSM_SN_SIGHTLINE_TOMOGRAPHY_V1.zip, hash
db704faa3bf6f164017af10ed80028f2799b9f7a1dc4854b88c7108a7dfe8988.
Historical foreground input checksums and the original localization statistic
remain missing; no byte-reproducible historical observational rerun is claimed.

The copied original snapshot manifests describe larger packages. Only this bundle's
root SHA256.json defines the delivered compact contents. Follow
submission/ARCHIVE_RELEASE_INSTRUCTIONS.md before any public deposit.
'''
    (bundle/'README.md').write_text(readme,encoding='utf8')
    archive(bundle,OUT/'UNIVERSE_SUBMISSION_REPRODUCIBILITY_PACKAGE.zip')
    (OUT/'SHA256.json').write_text(json.dumps(hashes(OUT),indent=2)+'\n')
    print(f'Prepared {OUT}; {len(hashes(bundle))} reproducibility files')

if __name__=='__main__': main()
