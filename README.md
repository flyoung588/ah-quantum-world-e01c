# AH Quantum World E01-C — public H2 reproduction

One manual GitHub-hosted standard `ubuntu-24.04` job, Python 3.12, single-thread,
15-minute job timeout. No schedule, push trigger, workflow rerun, self-hosted or
larger runner. Public non-sensitive reproduction only; not an original discovery.

## Physics

All-electron, nonrelativistic, infinite-nuclear-mass Born–Oppenheimer singlet H2.
Full configuration interaction in cc-pVDZ, cc-pVTZ, cc-pVQZ. Energies include
nuclear repulsion. Methods, thresholds, grid and reference are fixed in
`protocol.json` and `reference.json` before execution. No measured physical
properties, storage performance, atomic hydrogen inventory or useful energy
system claims follow from this molecule reproduction.

Independent reference: K. Pachucki, PRA 82, 032509 (2010),
[arXiv:1007.0322](https://arxiv.org/html/1007.0322), Eq.2 and the potential table.
Numerical FCI convergence is not complete-basis convergence or agreement with
the high-precision BO reference. De is not experimental D0.

## Reproduction

After separate resource approval: use an independent Python 3.12 environment,
install `requirements.txt`, then `python run_h2.py`. Original run dependency
freeze and platform details are saved with results. This instruction does not
authorize or initiate another Actions dispatch.

## Zero-additional-fee boundary

[GitHub standard hosted runners](https://docs.github.com/en/actions/reference/runners/github-hosted-runners)
are free for public repositories. This workflow never invokes Actions artifact
uploads, caches, Packages, LFS or paid custom/larger runners. Only ordinary logs
are emitted; retrieved results become small ordinary Git text files. This does
not assert that the owner's unrelated activity has no charges.

Only `workflow_dispatch` exists; do not automatically repeat a failed run.
All computations run on GitHub, not Mac or Hetzner. Their prior NO-GO decisions
are unchanged. Environment/script completion and physical error tests are
reported separately.

## Actual outcome

**NO-GO: incomplete validation and failed numerical/root QA.** The single
scientific job stopped at its predeclared 660-second scientific ceiling, inside
the 15-minute job limit. All 48 fixed-grid points were computed; QZ minimum and
cross-check were not completed. TZ/QZ residual thresholds failed, and the DZ
20-bohr point fails ground-root plausibility despite a converged solver flag.

See [actual results and limits](RESULTS.md), the original `results/` data,
`results/ground_root_anomaly.json`, and `results/execution_error.json`.
`python audit_results.py` performs only saved-data hashing/arithmetic checks,
not another scientific calculation. No workflow rerun was initiated.
