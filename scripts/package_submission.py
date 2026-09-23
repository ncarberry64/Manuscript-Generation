"""Deterministic generic-LaTeX source ZIP from the active manuscript inputs.

Run after building main.pdf. Exact source hashes avoid embedding a
self-referential final Git SHA. No observational inputs are included.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import zipfile

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "manuscript"


def package(destination):
    files = {}

    def include(relative):
        path = (MANUSCRIPT / relative).resolve()
        path.relative_to(MANUSCRIPT.resolve())
        key = path.relative_to(MANUSCRIPT).as_posix()
        if key in files:
            return
        files[key] = path.read_bytes()
        if path.suffix == ".tex":
            text = files[key].decode("utf-8-sig")
            for child in re.findall(r"\\input\{([^}]+)\}", text):
                include(child if child.endswith(".tex") else child + ".tex")
            for child in re.findall(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}", text):
                include(child)

    include("main.tex")
    include("references.bib")
    include("main.pdf")
    files["README.txt"] = (
        "Universe generic-LaTeX source package for author/referee review.\n"
        "Build: latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex\n"
        "Use full TeX Live or MiKTeX with Latin Modern and natbib.\n"
        "Author declarations await confirmation; archive deposit instructions are supplied separately.\n"
        "Not an author-approved submission; raw survey data are not included.\n"
        "SHA256.json identifies exact included source and PDF bytes.\n"
    ).encode()
    files["SHA256.json"] = (json.dumps({name: hashlib.sha256(data).hexdigest()
                                       for name, data in sorted(files.items())}, indent=2) + "\n").encode()
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(destination, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, data in sorted(files.items()):
            info = zipfile.ZipInfo(name, date_time=(2026, 9, 19, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(info, data)
    with zipfile.ZipFile(destination) as archive:
        assert archive.testzip() is None
        manifest = json.loads(archive.read("SHA256.json"))
        assert all(hashlib.sha256(archive.read(n)).hexdigest() == digest for n, digest in manifest.items())
    print(f"Verified {len(files)} archive entries: {destination}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "submission/universe-source.zip")
    args = parser.parse_args()
    package(args.output)
