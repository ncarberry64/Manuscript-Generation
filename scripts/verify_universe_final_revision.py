"""Cheap final-document and frozen-result checks. Does not refit science."""
from pathlib import Path
import csv, hashlib, json, math, re
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
def digest(path, normalized=False):
    b=path.read_bytes()
    return hashlib.sha256(b.replace(b'\r\n',b'\n') if normalized else b).hexdigest()

def main():
    snapshot=json.loads((ROOT/'artifacts/SUBMISSION_SCIENTIFIC_SNAPSHOT_v1.json').read_text())
    for name,sha in snapshot['files'].items():
        assert digest(ROOT/name,True)==sha, f'Frozen science changed: {name}'
    base=ROOT/'analyses/sn_sightline_tomography_v1'
    manifest=json.loads((base/'GIT_SNAPSHOT_MANIFEST.json').read_text())
    for name,sha in manifest.items():
        assert digest(base/name)==sha, f'Frozen sightline output changed: {name}'
    with (base/'interpretation/LIKELIHOODS.csv').open() as f:
        rows={r['test']:r for r in csv.DictReader(f)}
    expected={'P_original_heldout_release_only':(1208,2.985),
              'P_original_training_to_heldout':(1208,-.883),
              'P_training_to_heldout_fixed_external_PV_sensitivity':(1208,-4.909),
              'DES_original_heldout_release_only':(1467,.240),
              'DES_Pantheon_transfer':(1467,1.321)}
    for name,(n,v) in expected.items():
        row=rows[name]; assert int(row['rows'])==n
        assert round(float(row['delta_chi2']),3)==v
        for key in ('chi2_null','chi2_frozen','delta_chi2','logL_null','logL_frozen','mock_p_frozen_direction'):
            assert math.isfinite(float(row[key]))
        assert abs(float(row['chi2_frozen'])-float(row['chi2_null'])-float(row['delta_chi2']))<1e-9
        assert abs(2*(float(row['logL_null'])-float(row['logL_frozen']))-float(row['delta_chi2']))<1e-9
    maintex=(ROOT/'manuscript/main.tex').read_text(encoding='utf8')
    paths=[ROOT/'manuscript/main.tex']+[ROOT/'manuscript'/(s+'.tex') for s in re.findall(r'\\input\{([^}]+)\}',maintex)]
    source='\n'.join(p.read_text(encoding='utf8') for p in paths)
    labels=re.findall(r'\\label\{([^}]+)\}',source)
    assert len(labels)==len(set(labels)), 'Duplicate label'
    for ref in re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',source): assert ref in labels,ref
    for fig in re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}',source): assert (ROOT/'manuscript'/fig).exists(),fig
    bib=(ROOT/'manuscript/references.bib').read_text(encoding='utf8')
    keys=set(re.findall(r'@\w+\{([^,]+),',bib))
    for group in re.findall(r'\\cite\w*\{([^}]+)\}',source):
        for key in group.split(','): assert key.strip() in keys,key
    pdftext=(ROOT/'manuscript/build/final.txt').read_text(encoding='utf8')
    assert not re.search(r'TODO|\?\?|\bNaN\b|\bInf\b',pdftext)
    log=(ROOT/'manuscript/main.log').read_text(encoding='utf8',errors='replace')
    assert not re.search(r'undefined|multiply defined|Overfull|^!',log,re.M)
    count=0
    for name in ['universe_final_focused_pytest.xml','universe_final_additional_pytest.xml']:
        suite=ET.parse(ROOT/'artifacts'/name).getroot()
        cases=suite.findall('.//testcase'); count+=len(cases)
        assert not suite.findall('.//failure') and not suite.findall('.//error')
    abstract=re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}',maintex,re.S).group(1)
    assert 220<=len(abstract.split())<=260
    result={'scientific_snapshot_files_verified':len(snapshot['files']),
            'sightline_snapshot_files_verified':len(manifest),'fresh_tests_passed':count,
            'abstract_whitespace_words':len(abstract.split()),'primary_saved_comparisons_verified':len(expected),
            'pages':int(re.search(r'Output written on main.pdf \((\d+) pages',log).group(1)),
            'pdf_sha256':digest(ROOT/'manuscript/main.pdf'),
            'source_checks':'references, citations, figures, labels, no PDF TODO/??/NaN/Inf, no overflow',
            'limits':'No new full numerical or raw-observation rerun. Historical full suite is separate evidence.'}
    (ROOT/'artifacts/UNIVERSE_FINAL_REVISION_VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
