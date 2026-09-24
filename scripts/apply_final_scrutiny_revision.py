"""One-shot, asserted editorial revision; scientific owners are unchanged."""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
def read(p): return (ROOT / p).read_text(encoding="utf-8")
def write(p, s): (ROOT / p).write_text(s, encoding="utf-8", newline="\n")
def replace(p, old, new):
    s = read(p)
    assert old in s, (p, old)
    write(p, s.replace(old, new))

TITLE = "Geometry Before Fields: Constraint-Reduced Cosmology and Optical Predictions in the Berger--Hopf Framework"
oldtitle = "Constraint-Reduced Cosmology on a Closed Hypersphere: Berger--Hopf Optical Transfer and Falsifiable Cross-Observable Predictions"
abstract = r"""We develop a conditional cosmological response chain on a closed $S^3$
parent, resolving constraints and gauge redundancy before assigning physical
degrees of freedom. The regular one-scalar quotient removes the exceptional
$n=1$ sector; the first ordinary physical scalar carrier is $n=2$.
One selected normalized profile has two temporal coordinates; its pure
observer-sky dipole is proportional to $\sin(2\psi)\hat{\mathbf n}\cdot\hat{\mathbf p}$.
The frozen R1 action-native dust--radiation quadratic reduction has a positive
audited kinetic matrix; no ghost or high-$k$ scalar gradient instability is
detected within the audited linear/quadratic domain. Environmental transfer
has temporal rank two at eight sampled post-anchor epochs, without selecting
a spatial profile or axis. Schur reduction supplies the metric response,
while Jacobi/Sachs propagation transports optical means and covariances.
For $m$ coherent global amplitudes and fixed response,
$\operatorname{rank}C_{\rm topo}\le m$; local noise alone fixes no channel count.
Two prospective rejection gates test a shared deterministic state and a
minimal coherent one-amplitude covariance, each at $p<0.01$ under a
preregistered survey design. A frozen retrospective study reconstructs source
diagnostics for 3,521 supernova rows and adds independent foreground information.
With $\Delta\chi^2=\chi^2_{\rm frozen}-\chi^2_{\rm null}$, conditional held-out
Pantheon+ source/environment adjustment gives $-0.883$, whereas DES transfer
gives $+1.321$. A fixed external peculiar-velocity mean gives $-4.909$ for
Pantheon+, but incomplete joint uncertainty makes this sensitivity only.
Aggregate shell tests do not localize residual correlations; no coherent
topographic signal is established. Absolute source and microscopic
normalization, physical spatial selection, joint observational covariance,
nonlinear-fluid stability and ultraviolet completion remain open.
No scientific parameter is retuned.
"""
s = read("manuscript/main.tex").replace(oldtitle, TITLE)
s = re.sub(r"(?<=\\begin\{abstract\})[\s\S]*?(?=\\end\{abstract\})", lambda _: "\n"+abstract, s)
write("manuscript/main.tex", s)
for p in ["submission/UNIVERSE_SPECIAL_ISSUE_COVER_LETTER.md", "submission/UNIVERSE_PLANNED_PAPER_ABSTRACT_250_WORDS.md"]:
    s = read(p).replace(oldtitle, TITLE)
    write(p, s)
write("submission/UNIVERSE_PLANNED_PAPER_ABSTRACT_250_WORDS.md", "# "+TITLE+"\n\n"+abstract+"\nThe abstract is synchronized with manuscript/main.tex; mathematical notation is retained.\n")

p = "manuscript/sections/10_discussion.tex"
replace(p, "The physical selection and population mechanism remain open.", "The coupled environmental transfer generates both temporal components at\neight audited epochs (Section~\\ref{sec:coupled-environmental-state}); this\ndoes not select a common spatial profile, its population, or an observer axis.")
replace(p, "The propagator uses gravity--scalar reduced coefficients. Dust is transported\non that response with zero induced matter response at the anchor. This\nis not a fully backreacting matter+radiation stability theorem.", "The retained reference observable kernels use gravity--scalar reduced\ncoefficients and induced dust transport with zero response at the anchor.\nThe separate action-native coupled owner includes dust and radiation\nbackreaction at quadratic order. Its audited kinetic block is positive and\nits principal-symbol scan detects no high-$k$ scalar gradient instability.\nThese results do not silently replace the reference observable kernels.")
replace(p, "Full matter--radiation backreaction, an independently selected normalized\nstate direction, complete joint observational covariance and an independent\nprospective survey design also remain necessary.", "Nonlinear-fluid and ultraviolet control, independently reconstructed\nenvironmental anchor data and a normalized spatial profile, complete joint\nobservational covariance and an independent prospective survey design remain\nnecessary. Exact radial BAO and survey-ready RSD projections remain open.")
replace(p, "First, the current common-state basis is the canonical\n$X_2=(q_2,\\Pi_2)$.", "First, the common-state basis is the reference-normalized\n$X_2=(q_2,\\Pi_2)$. In the coupled system $\\Pi_2$ is not generally the\ncanonical momentum $p_\\zeta$, as the environmental section makes explicit.")
replace(p, "comparison with canonical state\ninference.", "comparison with reference-normalized state\ninference.")

p = "manuscript/sections/04aj_r1_n2_reduced_propagator.tex"
replace(p, "The current $\\mathcal F_{S,2}$ implementation is explicitly the gravity--scalar expression; the complete matter+radiation scalar stability/response system is not yet included.  Consequently this result is not yet the directional $f\\sigma_8$, lensing, or luminosity-distance observable kernel.", "This retained reference $\\mathcal F_{S,2}$ implementation is the gravity--scalar expression. The separate action-native coupled dust--radiation owner supplies the quadratic backreaction and audited stability results reported below; it does not retroactively change this reference kernel. This equation alone is not a directional $f\\sigma_8$, lensing, or luminosity-distance observable kernel.")
p = "manuscript/sections/04ak_r1_n2_growth_lensing_transfer.tex"
s = read(p); pos = s.index("\n", s.index(r"\section"))
s = s[:pos+1] + "\nThese retained reference rows use the stated gravity--scalar response and\ninduced dust transport. They are not silently replaced by the coupled\nenvironmental propagator or its reference-normalized momentum coordinate.\n" + s[pos+1:]; write(p,s)
for p in ["manuscript/sections/04b_matter_sourced_n2.tex", "manuscript/sections/04c_bianchi_dictionary_n2.tex"]:
    s=read(p); pos=s.index("\n",s.index(r"\section"))
    s=s[:pos+1]+r"""
\paragraph{Historical compatibility audit.}
The Gambino--Pace/curved effective equations below are retained as a
held-out compatibility calculation with a documented constraint/identity
defect. They do not define the current R1 dynamics. The action-native
ADM gravity--scalar plus ideal dust--radiation reduction owns the coupled
evolution; no failed identity is patched to obtain the results in this paper.
"""+s[pos+1:]; write(p,s)
p="manuscript/sections/01b_claim_status_and_falsifiability.tex"
replace(p,"DIAGNOSTIC & Raw",r"""NUMERICALLY AUDITED & Frozen R1 action-native dust--radiation quadratic
backreaction: positive audited kinetic matrix and no detected high-$k$
scalar gradient instability. Nonlinear-fluid stability, ultraviolet and
microscopic completion are not established.\\[0.5em]
DIAGNOSTIC & Raw""")
p="manuscript/sections/07_predictions_falsification.tex"
replace(p,"The bridge predicts selective", "The retained illustrative harmonic-filter realization predicts selective")
replace(p,"At large \\(n\\),", "Within the retained illustrative polynomial filter, at large \\(n\\),")
replace(p,"The model therefore predicts a strong separation", "This conditional filter realization therefore predicts a strong separation")
replace(p,"The minimal bridge is disfavored if", "The specified conditional realization (including its illustrative filter for the harmonic and high-$n$ tests) is disfavored if")

def pair(a,b,caption,label):
    return r"""
\begin{figure}[p]
\centering
\includegraphics[width=0.82\linewidth]{figures/"""+a+r"""}\\[0.5em]
\includegraphics[width=0.82\linewidth]{figures/"""+b+r"""}
\caption{"""+caption+r"""}
\label{"""+label+r"""}
\end{figure}
"""
p="manuscript/sections/02_hyperspherical_parent.tex"
write(p,read(p)+pair("n2_observer_shell_radial_basis.pdf","n2_observer_sky_dipole.pdf",r"Exact $n=2$ observer-shell radial factors (top) and the pure $A_{4i}$ angular dipole (bottom). The plotted axis is centered solely for visualization; no empirical RA/Dec, amplitude or redshift law is used. A general trace-free $4\times4$ tensor has nine independent components; the pure dipole is a conditional subspace.","fig:n2-exact-geometry"))
p="manuscript/sections/04ajb_coupled_environmental_state.tex"
write(p,read(p)+pair("coupled_environment_singular_values.pdf","coupled_environment_matter_determinant.pdf",r"Environmental transfer in the $(\zeta,\dot\zeta)$ output basis: singular values of $\widehat U_{XE}$ (top) and its matter-density/velocity subblock determinant (bottom). Points are the eight audited post-anchor epochs; lines only guide the eye. Both temporal components are accessible there. Singular-value magnitudes depend on coordinate normalization; rank is unchanged by the positive conversion to $(q_2,\Pi_2)$. The block vanishes at the anchor, which is not plotted. This is not spatial profile or axis selection.","fig:coupled-environment"))
p="manuscript/sections/09_results_R1.tex"
s=read(p).replace("and near-standard high-\\(n\\) structure growth.","and near-standard high-\\(n\\) structure growth within the retained illustrative filter realization.")
s+=pair("r1_omegaT.png","r1_growth_ratio.png",r"Conditional frozen R1 background: scalar energy fraction (top) and growth ratio (bottom), relative to matched curved $\Lambda$CDM with identical early normalization. These are reference calculations, not fits or survey-window RSD likelihoods. Near-standard high-$n$ recovery is conditional on the illustrative retained filter.","fig:r1-background-growth")
s+=r"""
\subsection{Audited action-native quadratic stability}
The minimally coupled ideal dust--radiation action is reduced with the
lapse/shift Schur complement. In the retained numerical audit on
$0\le z\le2.1$, its reduced kinetic matrix is positive; the extension
through $n=1280$ has minimum eigenvalue $1.21879873\times10^{-8}$.
The principal-symbol extrapolation gives positive propagating speed
squares, with radiation approaching $1/3$ and the non-radiation branch
bounded below by approximately $0.910373$ at the sampled epochs.
The positive normalized growth tail decreases as the harmonic sequence is
extended. No ghost or high-$k$ scalar gradient instability is detected in
this audited linear/quadratic domain. These numerical checks are not a
proof of nonlinear-fluid stability, a microscopic normalization, or an
ultraviolet completion. The detailed dispersion audit is retained in
the appendix; the gravity--scalar reference kernels above remain distinct.
"""
s+=pair("r1_principal_characteristic_speeds.pdf","r1_highk_growth_convergence.pdf",r"Frozen R1 coupled principal-symbol audit. Top: extrapolated characteristic speed squares at the audited redshifts, including the radiation $1/3$ reference. Bottom: the positive normalized growth tail at the two redshift endpoints decreases across the retained high-$n$ sequence. A positive finite-$n$ tail is not erased; convergence is numerical evidence within the audited quadratic domain, not an all-scale stability theorem.","fig:r1-principal-stability")
write(p,s)
p="manuscript/sections/09a_reference_figures.tex"
s=read(p)
for name in ["r1_omegaT","r1_growth_ratio"]:
    s=re.sub(r"\\begin\{figure\}[^\n]*\n(?:(?!\\end\{figure\})[\s\S])*?figures/"+name+r"\.png(?:(?!\\end\{figure\})[\s\S])*?\\end\{figure\}\n", "", s)
write(p,s.replace("\\clearpage\n\n\\begin{figure}","\n\\begin{figure}"))

p="code/final_scientific_figures.py"
s=read(p).replace('import matplotlib.pyplot as plt','import matplotlib.pyplot as plt\nimport hashlib\nplt.rcParams.update({"font.size": 12, "axes.titlesize": 12, "legend.fontsize": 10})')
s=s.replace('bbox_inches="tight")','bbox_inches="tight", metadata={"CreationDate": None, "ModDate": None})',1)
s=s.replace('(7.2, 4.8)','(6.4, 3.6)').replace('(7.2, 4.2)','(6.4, 3.5)')
s=s.replace('Exact observer-shell factors of the physical $n=2$ scalar harmonic','Exact $n=2$ observer-shell radial factors')
s=s.replace('s_1(U_{XE})',r's_1(\widehat U_{XE})').replace('s_2(U_{XE})',r's_2(\widehat U_{XE})').replace(r'\det U_{XE}',r'\det \widehat U_{XE}')
s=s.replace(r'Environmental $\rightarrow$ topographic temporal transfer is rank two','Sampled environmental temporal transfer')
s=s.replace('Matter density and velocity alone span both temporal components','Sampled matter density/velocity subblock')
s=s.replace('High-k positive-growth tail decreases with harmonic scale','High-$k$ positive-growth convergence')
s=s.replace('    fig, ax = plt.subplots(figsize=(6.4, 3.6))\n    ax.plot(z, s1', '    assert len(rows) == 8 and np.all(np.isfinite([s1, s2, det]))\n    assert np.all(s2 > 0) and np.all(det > 0)\n    fig, ax = plt.subplots(figsize=(6.4, 3.6))\n    ax.plot(z, s1')
s=s.replace('    print(f"Final scrutiny figures written to {OUT}")', '''    inputs = [COUPLED, STABILITY, Path(__file__)]
    outputs = [OUT / (stem + ext) for stem in [
        "n2_observer_shell_radial_basis", "n2_observer_sky_dipole",
        "coupled_environment_singular_values", "coupled_environment_matter_determinant",
        "r1_principal_characteristic_speeds", "r1_highk_growth_convergence"] for ext in [".pdf", ".png"]]
    receipt = {"retuning": False, "coupled_output_basis": ["zeta", "zeta_dot"],
               "rank_scope": "eight sampled post-anchor epochs; anchor block is zero",
               "analytic_axis": "visualization only; no empirical sky direction",
               "sha256": {str(p.relative_to(ROOT)).replace("\\\\", "/"): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs + outputs}}
    (ROOT / "artifacts/provenance/FINAL_SCIENTIFIC_FIGURES_RECEIPT.json").write_text(json.dumps(receipt, indent=2)+"\\n")
    print(f"Final scrutiny figures written to {OUT}")''')
write(p,s)

p="docs/REVIEWER_RISK_AUDIT.md"
replace(p,"Discussion: no complete backreacting matter+radiation stability theorem. Reference gravity-scalar response only. | artifacts/R1_n2_growth_lensing_transfer.json", "Frozen R1 action-native quadratic dust/radiation backreaction: positive audited kinetic matrix and no detected high-k scalar gradient instability. Nonlinear/UV limits remain open; legacy observable kernels remain reference calculations. | artifacts/R1_ACTION_NATIVE_GRADIENT_DISPERSION_CONVERGENCE_V1.json")
replace(p,"full matter/radiation, state direction", "nonlinear-fluid/UV control, environmental anchor data and spatial profile")
p="docs/HOSTILE_REFEREE_AUDIT.md"
replace(p,"| Two state dimensions are not action-owned | OPEN |", "| Two temporal coordinates do not select a spatial profile | PARTIALLY CLOSED | The action-native environmental transfer has rank two at eight audited epochs, independently of spatial selection.")
replace(p,"and complete matter/radiation backreaction remain open", "and nonlinear-fluid control remain open; quadratic coupled backreaction is action-native and audited")
replace(p,"Discussion now separates illustrative polynomial filter from a controlled large-n result of the exact complement. UV completion and full stability are OPEN.", "The polynomial harmonic filter is illustrative. The separate frozen R1 principal-symbol scan detects no ghost/high-k scalar gradient instability within its audited quadratic domain. A controlled exact-complement UV completion and nonlinear stability remain OPEN.")
p="docs/UNIVERSE_CLAIM_PROVENANCE.md"
replace(p,"Gravity-scalar reference; full matter/radiation stability open", "Retained gravity-scalar reference rows; separate coupled quadratic stability audit does not replace these observable kernels")
s=read(p)+"\n## Final scrutiny additions (2026-09-24)\n\n- Coupled temporal rank and basis conversion: `scripts/audit_r1_coupled_environmental_state.py`, `artifacts/provenance/R1_COUPLED_ENVIRONMENT_INTEGRATION_REPLAY.json`. Rank two only at eight audited post-anchor epochs; no spatial selection.\n- Quadratic stability: `scripts/audit_action_native_gradient_dispersion_v1.py`, `artifacts/R1_ACTION_NATIVE_GRADIENT_DISPERSION_CONVERGENCE_V1.json`. Positive audited kinetic matrix and no detected high-k scalar gradient instability; nonlinear/UV completion open.\n- Exact geometry and new audit figures: `code/final_scientific_figures.py`, `artifacts/provenance/FINAL_SCIENTIFIC_FIGURES_RECEIPT.json`. Exact analytic factors and committed numerical artifacts; no observational axis in the basis figure.\n- Earlier fresh/historical test descriptions above refer to their original validation passes. Current validation and hashes are recorded in `docs/FINAL_SCRUTINY_VALIDATION_20260924.md`.\n"
write(p,s)
p="README.md"
s=read(p).replace("# Geometry Before Fields: Manuscript-Generation", "# "+TITLE)
s=s.replace("| REJECTED / NOT ESTABLISHED |", "| NUMERICALLY AUDITED | Frozen R1 action-native quadratic dust/radiation backreaction; positive audited kinetic matrix and no detected high-k scalar gradient instability. Environmental transfer has temporal rank two at eight post-anchor epochs; no spatial profile/axis selection follows. |\n| REJECTED / NOT ESTABLISHED |")
s=s.replace("python -m pytest -q", "$env:PYTHONPATH='code'\npython -m pytest -q").replace("python code/figures.py", "python code/figures.py\npython code/final_scientific_figures.py")
s=s.replace("## Submission checkpoint", "## Historical submission checkpoint")
s+="\nCurrent review branch: `review/cosmology-final-scrutiny-20260924`, based on public `083e14654d1b81b3853484b48e132bd5bc41567c`. See [final scrutiny validation](docs/FINAL_SCRUTINY_VALIDATION_20260924.md). No merge or submission without author approval.\n"
write(p,s)
print("Final scrutiny editorial and figure revision applied; no scientific owner edited.")
