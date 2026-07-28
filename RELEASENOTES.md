# splikit — Release Logbook

A chronological logbook of notable changes. Each entry is dated and signed by the
person who made the change, with a short summary of what was done and why. For the
canonical, CRAN-facing changelog see `NEWS.md`.

---

## 2.3.3 — 2026-07-27

**Summary — Exportable pseudo-correlation null results.**

- Replaced the pseudo-correlation result's `null_draws` matrix attribute with a
  structured S3 result containing `statistics`, a long `null_distribution`, and
  computation `metadata`.
- The long table has one row per retained event and permutation and preserves
  both identifiers for reliable export and reconstruction.
- Per-event empirical p-values and adjusted p-values are unchanged; the pooled
  null column is descriptive and is not used as a global event-level null.
- Added independent `matrix_source` and `matrix_file` controls to
  `make_gene_count()`, including direct support for raw or filtered matrices
  and STARsolo files such as `UniqueAndMult-EM.mtx`. The defaults now read
  `raw/UniqueAndMult-EM.mtx` when available and fall back to `raw/matrix.mtx`.

---

## 2.3.2 — 2026-06-09

**Author:** Arsham Mikaeili Namini

**Summary — Repeated permutation null for pseudo correlation.**

- Extended `get_pseudo_correlation()` (and the R6 method
  `SplikitObject$getPseudoCorrelation()`) to build the null model from multiple
  cell-permutations instead of a single one.
- Added two arguments:
  - `permutation_count` (default `100`): number of times the cells (columns) of
    `ZDB_matrix` are permuted and the pseudo correlation recomputed, producing a
    per-event empirical null. `permutation_count = 1` reproduces the original
    single-permutation behaviour.
  - `permutation_seed` (default `NULL`): seeds the permutations for reproducible
    nulls and empirical p-values. The caller's global RNG stream is saved and
    restored, so a seed here does not disturb downstream randomness.
- The observed correlation is computed once; only the null pass scales with
  `permutation_count`.
- New output columns: `null_sd`, `n_perm_valid`, `emp_pvalue` (two-sided,
  `(b + 1) / (m + 1)` corrected), `emp_padj` (Benjamini-Hochberg). The full
  per-event null draws are attached as `attr(result, "null_draws")`.
- `null_distribution` is now the mean of the per-event null draws (equal to the
  single draw when `permutation_count = 1`); `pseudo_correlation` is unchanged.
- Tests: added a dedicated permutation test (validity of empirical p-values,
  seed reproducibility, global-RNG isolation, `permutation_count = 1`
  equivalence, input validation); the full-data R6 test runs with a small
  `permutation_count` to stay fast.
- Version bumped to 2.3.2; `NEWS.md` updated.

— Arsham Mikaeili Namini
