# External observational inputs

This directory documents inputs; it does not contain survey data. No unrelated BHSM repository is modified or required to run the local numerical tests.

| Input | Primary source / retained provenance | Local requirements and status |
|---|---|---|
| Pantheon+ | [Dataset paper](https://arxiv.org/abs/2112.03863), [cosmological analysis](https://arxiv.org/abs/2202.04077) | Foreground code needs compatible sky coordinates, redshift, magnitude/residual uncertainties and available controls. Exact cleaned historical table and hash missing; 277 selected objects in artifacts. |
| 2MRS density catalog | [Huchra et al. release](https://arxiv.org/abs/1108.0669) | Historical RA/Dec/z CSV; exact hash missing. Number-density result used 38600 objects. |
| 2MRS photometry | [VizieR release description](https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/ApJS/199/26?format=html&tex=true), table3 | RAJ2000, DEJ2000, z=cz/299792.458, Kcmag renamed Kmag; 0<=z<0.05. Temporary derived CSV not bundled; 38601 galaxies. See `docs/r1_2mrs_luminosity_catalog_provenance.md`. |
| BHSM upstream action | `docs/bhsm_native_cosmology_bridge.md` | Local note cites upstream PR numbers; immutable source exports and full microscopic normalization are not supplied. Treat imported ownership as conditional. |

For each recovered table record source release/version, retrieval date, SHA-256 of raw and derived files, selection script/commit, column units, frame, masks and row counts. Preserve original inputs. A later re-download does not automatically reproduce the historically used table.

No future survey targets have been selected or unblinded for the new two-gate protocol.

## Retrospective SN sightline V1 inputs

The new [sightline package](../analyses/sn_sightline_tomography_v1/GIT_SNAPSHOT.md) is distinct from the historical foreground tests above. It pins Pantheon+ commit `c447f0fea703fcd0fff57de5000947b5ca81286b` and DES commit `c9a4fcafc4cbd19bd750dee47fc76194a45c181f`, hashes raw photometry/models and official HD covariance inputs, and records independently downloaded 2MRS-NeuralNet fields and VizieR galaxy counts. Complete raw/derived inputs are retained in its local archive; compact provenance is committed here. A public permanent archive has not been created. These inputs do not silently replace any older frozen catalogue.
