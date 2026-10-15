# Actual single-run H2 results

Run: https://github.com/flyoung588/ah-quantum-world-e01c/actions/runs/38733218903
Dispatched source commit: `34809b009995ab59d8e9ee4b132cd0b7bf3c640d`
Workflow conclusion: failure
Scoped outcome: **E01_C_NO_GO_INCOMPLETE_AND_NUMERICAL_QA**
High precision outcome: **NOT_ESTABLISHED_INCOMPLETE**

A workflow or finite-basis numerical PASS does not establish complete-basis BO accuracy. This is known H2 reproduction, not an original discovery.

| Basis | Re (bohr) | Re (angstrom) | Emin (Eh) | same-basis De (Eh) | error at 1.4011 (microEh) |
|---|---:|---:|---:|---:|---:|
| cc-pVDZ | 1.43788022 | 0.76089344 | -1.163672981196 | 0.1651161744 | 11061.179 |
| cc-pVTZ | 1.40335368 | 0.74262279 | -1.172336686763 | 0.1727170642 | 2140.189 |

## Actual matched-distance basis comparison

These are computed values, not all numerically qualified values. Fixed R = 1.4011 bohr; independent BO reference -1.1744759314002167 Eh.

| Basis | computed total energy (Eh) | BO deviation (mEh) | residual norm (Eh) | frozen residual gate |
|---|---:|---:|---:|---|
| cc-pVDZ | -1.163414752076 | 11.061179 | 2.806e-15 | PASS |
| cc-pVTZ | -1.172335742664 | 2.140189 | 6.378e-07 | FAIL |
| cc-pVQZ | -1.173796348769 | 0.679583 | 8.908e-07 | FAIL |

## Error interpretation

- Method: full two-electron finite-orbital-space CI; this removes correlation truncation within that basis, not basis incompleteness.
- Basis: report DZ/TZ/QZ shifts and actual independent BO errors, not a guaranteed CBS extrapolation. Because numerical/root QA fails for parts of this run, do not attribute every discrepancy solely to basis error.
- Numerical: all CI residuals, spin values, tolerance changes and alternative orbital/solver checks are retained. Local minimum step sensitivity is separate.
- Physical model: fixed infinite-mass nuclei, nonrelativistic BO; no relativity, QED, nonadiabatic corrections or ZPE. De is not experimental D0.
- Endpoint: primary De uses two separate same-basis H atoms. The -1 Eh exact-endpoint diagnostic is distinct. No counterpoise correction applied.
- Reference distance 1.4011 bohr is rounded near equilibrium, not an independently known exact Re; its depth is a reference-point proxy.

## Resources and versions

Python 3.12.15; PySCF 2.14.0; ASE 3.29.0; NumPy 2.3.5; SciPy 1.16.3; h5py 3.15.1.
Computation wall time 660.629 s; process user CPU 654.777 s; peak RSS 815.13 MiB.
All numerical thread pools are checked against one thread; runner CPU affinity is restricted to one core.

## Provenance and cost boundary

One manual workflow_dispatch, attempt 1. A prior invalid-YAML configuration record had zero jobs; it was not rerun.
Only non-sensitive H2 code/results are public. No XREM production repository or Mac/Hetzner computation was used.
Standard public ubuntu-24.04 runner; no artifact uploads, cache, Packages, LFS, self-hosted or larger runner.
Raw scientific output files were recovered from logs and SHA256 verified, not regenerated.
No automatic continuation or workflow rerun. Prior Hetzner/Mac NO-GO decisions remain unchanged.

## Incomplete run

Only actually completed minima and points appear above. Missing QZ/minimum or validation outputs are NOT_TESTED, not replaced by extrapolation.
All 48 fixed-grid energies were computed, but QZ minimum and its cross-check were not completed. See execution_error.json and recovered_curve_from_actual_evaluations.json.
TZ/QZ residuals exceed the frozen 1e-7 Eh threshold; converged=True is not a numerical PASS.
DZ at 20 bohr is -0.519859847896 Eh versus independent BO -1.000000106740 Eh, and is inconsistent with the same-basis separated-atom limit near -0.998556806839 Eh. This is a ground-root plausibility failure, not a legitimate high-accuracy ground-state point.
Initial-guess/root selection is a suspected cause, not a confirmed diagnosis from an additional calculation. No repair or workflow rerun was performed.
