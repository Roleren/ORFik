# Transcript mapping: core and R bindings

Measured on this server, 2026-10-10, ORFik 1.29.15. The baseline is the
previously optimized pointer implementation, not the older unbounded mapper.
Both sides use `-O2`. Prototypes were benchmarked before integration; final
numbers below are from the installed implementation. Timings are medians of
five alternating paired batches, with collection outside the timed region.

## What changed

`inst/include/ORFik/transcript_mapping_core.h` is a standalone C++11 header.
It contains validated exon catalogues, bounded coordinate searches, mapping
plans and interval merging. It has no R/Rcpp dependencies. `IntView` borrows
integer arrays; its owner must remain alive. Coordinates are inclusive and
one-based; catalogue group indices are zero-based. A Python binding could
borrow compatible integer arrays, translate coordinate conventions and wrap
errors; a Python binding is not implemented here.

`src/transcript_mapping_R.cpp` owns R-specific input decoding, flag/name
validation, factor/Rle handling and output factories. It borrows coordinate
storage, indexes each reference once, handles both strands in one pass and
constructs complete Bioconductor objects from cloned empty prototypes.
From-transcript mapping no longer duplicates reference exons for every query.
To-transcript interval merging no longer constructs an intermediate S4 result
and sends it through R `reduce()`. Explicit optional sorting still delegates
to ORFik's existing sorter. R entry points are thin bindings; legacy numeric
native signatures and schemas remain available.

## Full calls

Actual ATF4/FMR1/BOD1L1 transcript structures have 3/17/26 exons. Repeat this
three-transcript set to obtain larger inputs. To-transcript inputs are the
transcripts themselves, with default interval reduction. From-transcript
inputs are terminal 11-nt windows with `removeEmpty=TRUE`. Times include
preparation, validation, mapping, reduction and complete output construction;
the empty prototypes are initialized before timing.

| Call | Transcript count | Before ms | After ms | Speedup |
|---|---:|---:|---:|---:|
| To | 3 | 17.640 | 0.060 | 294× |
| To | 300 | 20.733 | 0.467 | 44.4× |
| To | 6,000 | 104.800 | 9.800 | 10.7× |
| From | 3 | 16.400 | 0.040 | 410× |
| From | 300 | 16.567 | 0.200 | 82.8× |
| From | 6,000 | 29.600 | 3.000 | 9.9× |

These are specific mapping workloads, not whole-pipeline speedup estimates.
The large ratios on tiny inputs reflect removal of R dispatch/construction
costs. Repeated references and default exon reduction favor this implementation;
different query shapes should be measured separately.

## Rcpp versus borrowed and owned arrays

A separate data-independent benchmark uses 100,000 valid last-exon queries.
Each mode uses the same catalogue/search logic. Rcpp mode reads and writes
through vector element access; borrowed mode uses integer pointers; owned
modes copy the inputs into std::vectors, optionally also mapping into owned
output vectors before exporting to R. Below, milliseconds per call:

| Exons | Borrowed pointers | Rcpp element loop | Copied inputs | Copied inputs and outputs |
|---|---:|---:|---:|---:|
| 3 | 1.26 | 1.78 | 1.33 | 1.40 |
| 20 | 1.35 | 1.88 | 1.41 | 1.48 |
| 200 | 1.73 | 2.24 | 1.82 | 1.87 |

Borrowed access is approximately 1.3–1.4× faster than the element loop in
these cases. Copying adds about 4–6% for inputs and 8–11% including outputs.
Rcpp remains useful for protected ownership; transferring everything into
std::vectors is not a speed optimization for R callers.

The legacy to-transcript numeric entry point now uses the common catalogue
and has a small extra preparation cost: 100,000-query cases at 3/20/200 exons
changed from 0.80/1.06/1.50 ms to 1.04/1.16/1.56 ms. This regression is retained
for a shared validated core; full public calls show the gains above.

## Verification and traps for future work

- Existing genomic-range and boundary suites: 339 assertions passed. New
  factory, immutability, malformed-input and reduction tests: 36 passed.
- 492 full-object comparisons, including randomized inputs, plus scalar, empty, metadata,
  unmapped-reduction checks passed against independent original native code.
- dominant_cell_states tests with local RiboCrypt loaded: 228 passed.
- `test_core.cpp` compiled and ran without R/Rcpp under UndefinedBehaviorSanitizer.
- Compare complete objects and run `validObject(..., complete=TRUE)`, not just
  start/end coordinates. Prototypes must remain immutable and outputs independent.
- Do not time load_all's default `-O0` debug build against sourceCpp's `-O2`.
- Saved R .Call wrappers resolve the current native library. Compile old native
  code under independent symbols, as run.R does, for a genuine baseline.
- Hold every wrapped output in a protected Rcpp vector before allocating the
  next one. An early owned-array benchmark returned two unprotected raw `wrap()`
  results and crashed during GC; the corrected benchmark uses protected holders.
- Named references currently produce NA output seqlengths in the established
  wrapper. That quirk is preserved for equivalence; changing it is a separate
  correctness decision. Flag validation is now scalar, logical and non-missing.

## Reproduce

From the ORFik repository root:

```bash
MAKEFLAGS=-j2 R CMD INSTALL --preclean .
env -u LC_ALL R --vanilla -q -f benchmark/transcript_mapping/run.R
env -u LC_ALL R --vanilla -q -f benchmark/transcript_mapping/storage.R
g++ -std=c++11 -O2 -Wall -Wextra -Werror -fsanitize=undefined \
  -fno-sanitize-recover=all benchmark/transcript_mapping/test_core.cpp -o /tmp/orfik_mapping_core_test
/tmp/orfik_mapping_core_test
```

run.R also runs validate.R. `DCS_COVERAGE_CACHE` selects the external gene
caches; data are not packaged here. `ORFIK_MAPPING_OUTPUT` selects the output
folder, otherwise ignored `results/` is used. `ORFIK_MAPPING_VALIDATE_ONLY=1`
skips timings. storage.R needs no annotation data or ORFik package loaded.
Small final summaries are retained beside this note; raw timings, build logs
and prototype evidence remain in dominant_cell_states's ignored
`results/orfik_mapping_bindings_2026-10-10/`.

Remaining costs are catalogue construction, coordinate-end buffers, Rle
expansion and per-group sorting/reduction on large inputs. An explicit reusable
catalogue API could amortize those costs across repeated calls. It would need
clear reference lifetime/invalidation rules; no implicit reference cache was
introduced here. The object factories depend on Bioconductor slot layouts,
so rerun complete-object tests after dependency updates.
