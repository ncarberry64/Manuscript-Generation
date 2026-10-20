#!/usr/bin/env python3
"""Build the Galaxies strong-lensing manuscript with the unmodified MDPI class.

The official class names its EPS logos explicitly. Convert those original assets
to the filenames expected by epstopdf before running pdfLaTeX, or use the supplied
converted official assets when Ghostscript is absent. No class stub is used.
"""

from pathlib import Path
import shutil
import subprocess


def main() -> None:
    manuscript = Path(__file__).resolve().parents[1] / "manuscript"
    ghostscript = next((path for name in ("gs", "gswin64c", "gswin32c")
                        if (path := shutil.which(name))), None)
    latexmk = shutil.which("latexmk")
    if not latexmk:
        raise SystemExit("This build requires latexmk on PATH and the real MDPI class dependencies.")
    engine = shutil.which("pdflatex")
    engine_version = subprocess.run([engine, "--version"], capture_output=True, text=True,
                                    check=True).stdout if engine else ""
    # MiKTeX otherwise may wait for a graphical missing-package prompt even
    # in nonstop mode. This flag applies to this build, not global settings.
    engine_options = (["-pdflatex=pdflatex -enable-installer %O %S"]
                      if "MiKTeX" in engine_version else [])

    for source in sorted((manuscript / "Definitions").glob("*.eps")):
        target = source.with_name(source.stem + "-eps-converted-to.pdf")
        if not ghostscript:
            if target.is_file():
                print(f"Using supplied converted official asset: {target.name}", flush=True)
                continue
            raise SystemExit(f"Ghostscript is required to convert the missing asset {target.name}.")
        print(f"Converting official asset: {source.name}", flush=True)
        subprocess.run(
            [
                ghostscript,
                "-q",
                "-dSAFER",
                "-dBATCH",
                "-dNOPAUSE",
                "-sDEVICE=pdfwrite",
                "-dEPSCrop",
                "-dAutoRotatePages=/None",
                f"-sOutputFile={target}",
                str(source),
            ],
            check=True,
        )

    print("Building manuscript/strong_lensing_main.tex", flush=True)
    subprocess.run(
        [
            latexmk,
            *engine_options,
            "-pdf",
            "-g",
            "-interaction=nonstopmode",
            "-halt-on-error",
            "-file-line-error",
            "strong_lensing_main.tex",
        ],
        cwd=manuscript,
        check=True,
    )


if __name__ == "__main__":
    main()
