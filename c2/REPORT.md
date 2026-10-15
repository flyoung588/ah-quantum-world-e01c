# E01-C2 — actual incomplete numerical certification

**E01-C2 = NO-GO**
NUMERICAL_SOLVER_QUALIFIED = INCOMPLETE
H2_PHYSICAL_BASELINE = NOT_QUALIFIED
ORIGINAL_DISCOVERY_CAPABILITY = NOT_TESTED
Historical E01-C = NO-GO unchanged; compute environment remains qualified,
historical numerical failure and unestablished high-precision BO remain unchanged.

## Actual run and partial results

Run: https://github.com/flyoung588/ah-quantum-world-e01c/actions/runs/38741037586
Source SHA: `16441c0acfc0855df245c54f8ee880d22d8c947d`. One manual job, attempt 1, no rerun.

| Test | Actual status | Evidence retained | Not established |
|---|---|---|---|
| A / DZ / 20 bohr | PARTIAL_NOT_CERTIFIED | complete product-space eigenvalues; original pspace spectrum and singlet projections; two CI vectors | final root/spin/residual gates; full A certification |
| B / TZ / 1.4011 bohr | NOT_RUN | none from C2 | all requested C2 certification |
| C / QZ / 1.4011 bohr | NOT_RUN | none from C2 | all requested C2 certification |

The independent full product-space minimum is -0.9985568133616696 Eh. The first two
original pspace energies are -0.9985568133616707 and
-0.9985568133616702 Eh; their singlet projection norms are
0.0453301266450693 and 0.9989720614803713.
The original filter rejects both, then accepts index 2 at
-0.5198598478956922 Eh — matching the historical wrong energy.
Independent full/pspace spectra agree within
5.773e-15 Eh.
Thus a concrete degeneracy/purity-filter root-selection mechanism is supported;
it is not sufficient to blame only Davidson initial guesses. DZ fits inside
the default dense pspace shortcut. No missing Coulomb term is indicated by this
matched full-space spectrum. This is not a blanket certification of PySCF.

The full-space minimum differs from the historical same-basis two-atom endpoint
by -6.523e-09 Eh,
and from the independent BO reference by
1.443293378e-03 Eh. These support the asymptotic
diagnostic, NOT a completed singlet/residual certification. Saved-vector norm
and spin scalar checks are separately in `results/partial_saved_data_QA.json`.

## Failure introduced by this implementation

After A's vector and spectra files were saved, `dump(A_result.json,row)` raised
TypeError: a NumPy bool in the gate dictionary is not JSON serializable.
The Python exception string calls its type `bool`; the static expressions
comparing Python floats with NumPy eigenvalues produce NumPy scalar booleans.
This is a result-persistence programming defect, not a physical or convergence
failure. It nevertheless invalidates claiming completion: A's metrics were not
saved, and B/C never started. The exception, full traceback and original run logs
are preserved. No lost scalar metrics have been fabricated or recomputed.

The static preflight parsed Python/JSON/YAML and checked workflow/resource
boundaries but missed NumPy scalar serialization. That preflight was inadequate
for persistence QA; it is not renamed a numerical PASS. A minimal normalization
patch is supplied as **UNEXECUTED**; no scientific source was altered after the
run and no new scientific job is authorized or started in this tranche.

## Original audit and attempted minimal numerical repair

Original Hamiltonian, integral permutation, canonical AO orthogonalization and
electronic residual convention were audited in `STATIC_AUDIT.txt`. Original
`run_h2.py`, protocol, reference and all historical outputs remain byte-identical.
New code uses RHF MOs, explicit symmetric independent seed, direct_spin1,
pspace_size=0, conv_tol=1e-14 and conv_tol_residual=1e-10. The frozen acceptance
gate remains 1e-7. Independent assembly is
H_ij,kl=h_ik delta_jl+delta_ik h_jl+(ik|jl), projected into symmetric/antisymmetric
spatial sectors and solved with LAPACK. It shares the AO integral engine;
this is algorithmic, not independent physical/integral-engine validation.

PySCF 2.14 source shows the default Davidson residual limit sqrt(conv_tol).
This explains why 1e-12 energy control need not pass the 1e-7 residual QA.
But the attempted repair's B/C effectiveness is **NOT_TESTED**, not assumed.

## Physical error ceiling

BO comparison: [Pachucki, PRA82 032509 (2010)](https://arxiv.org/html/1007.0322),
total nonrelativistic infinite-mass BO potential including 1/R, matched R.
Finite Gaussian basis incompleteness is distinct from numerical residual and
method implementation error. Old TZ/QZ offsets are historical only, not new
C2 data. No CBS, microhartree, Re, De or experimental D0 conclusion follows.
Nuclear dynamics, nonadiabatic/adiabatic, relativistic/QED and ZPE effects omitted.

## Actual resources / permissions

Python 3.12.15, PySCF 2.14.0, ASE 3.29.0,
NumPy 2.3.5, SciPy 1.16.3; full freeze retained.
Scientific wall 1.876s; CPU user
0.492s, system 0.079s;
peak RSS 107.38 MiB. Whole job 28s.
One core/thread, hard 6 GiB address-space cap, 720s scientific deadline, 15min
job limit. No local or Hetzner electronic calculation; no production changes.
No artifacts/cache/Packages/LFS; small results through ordinary logs and Git.
[Public standard runner free policy](https://docs.github.com/en/billing/concepts/product-billing/github-actions)
was freshly verified. Zero new purchases, no claim about unrelated account costs.
Provider-total tokens and complete Mac coordination CPU usage are UNKNOWN.
No independent reviewer; all interpretation in the same coordinating context.

## NEXT_SCIENTIFIC_CAPABILITY — list only, no continuation

1. First, complete durable numeric certification with tested result persistence;
   this tranche has NOT certified even the three-point H2 baseline.
2. Transfer/error-domain evidence beyond two electrons: spin/root tracking,
   size consistency and independently checked many-electron approximations.
3. Coupled electron/nuclear dynamics, finite temperature and metastable reaction
   paths; fixed BO points cannot establish lifetimes or phases.
4. Grounded hosts/loading/boundaries, macroscopic inventory and reference/free
   energy accounting including constraint/environment burdens.
5. Useful release pathways and falsifiable predictions versus strong alternative
   mechanisms; no storage/output capability inherited from H2.
6. Independent implementations/integral engines, uncertainty propagation and
   reality-facing benchmarks/experiments for expanded domains.
7. Bounded unknown-hypothesis generation and discriminating test selection;
   numerical functionality alone never proves autonomous original discovery.

STOP. No rerun, automatic development, new Local Game or production adoption.
