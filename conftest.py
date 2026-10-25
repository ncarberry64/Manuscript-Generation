"""Expose the existing numerical modules consistently to pytest and CI."""
from pathlib import Path
import sys

CODE = str(Path(__file__).resolve().parent / "code")
if CODE not in sys.path:
    sys.path.insert(0, CODE)
