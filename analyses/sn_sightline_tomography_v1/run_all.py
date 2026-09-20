"""Offline numerical reproduction from the bundled, hash-pinned public inputs."""
from pathlib import Path
import subprocess,sys
ROOT=Path(__file__).resolve().parent
for stage in ['sightline_stage.py','lightcurve_stage.py','source_audit.py','epoch_jackknife.py','ambiguous_sources.py','enrich_sightlines.py','interpretation_stage.py','validate_package.py','write_report.py']:
 print('Running',stage,flush=True)
 with (ROOT/(stage.replace('.py','')+'_REPRODUCTION.log')).open('w') as log:
  subprocess.run([sys.executable,str(ROOT/stage)],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
print('Complete. Read FINAL_REPORT.md; package_outputs.py refreshes the distributable archive.')
