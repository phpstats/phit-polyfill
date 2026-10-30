# PHP 8.3's JIT computes wrong values and crashes

Recorded 2026-10-17 on PHP 8.3.33 (Laravel Herd, macOS arm64) and PHP 8.3.35
(the official `php:8.3-cli` image on Linux arm64 and on x86-64 under
emulation). Nothing has been reported upstream. PHP 8.3 receives security
fixes only, and its JIT is the DynASM backend that PHP 8.4 replaced with the
IR-based one, so fixes there are unlikely.

## Recommendation

On PHP 8.3, run the polyfill with OPcache and without the JIT
(`opcache.jit=disable` or `opcache.jit_buffer_size=0`, which is PHP 8.3's
default). The polyfill refuses to run under the JIT on every version. `opcache.jit=function` (`1205`) passed every check below, but the
tracing JIT and the function JIT with hot-counter or first-call triggers did
not. PHP 8.4 and 8.5 are not affected by the problems described here.

## What happened

Each test file was run on its own under each JIT setting, because the
failures depend on which code becomes hot first. The full suite in one
process passed under the tracing JIT on macOS and hid them.

- **Tracing JIT, wrong Vᵀ.** `Svd::assembleU` and `Svd::assembleVt` assigned
  `$shift = $upper ? 0 : 1` before a loop. Under the tracing JIT, PHP printed
  "Undefined variable $shift" inside the loop of `assembleU`, and `assembleVt`
  returned a wrong matrix without a message. `DenseTest` failed when run
  alone, and 13 of 90 SVD, solve, and eigen cases differed from the run
  without the JIT on all three platforms. Assigning `$shift` from a cast
  (`(int) ! $upper`) removed every difference.

- **Tracing JIT, wrong beta CDF.** After `Phit\dist_beta_quantile` has run on
  beta(1e5, 1e5), `Phit\dist_beta_cdf(0.5, 1e5, 1e5)` returns 0.4999999999988
  instead of 0.49999999999999994. `DistributionsStressTest` fails under every
  hot-counter setting tried. A root trace that starts where `rlog1()` returns
  inside `Cdflib::betaAsym` and inlines `errorFcScaled()` is involved, but
  splitting that expression and moving `betaAsym`'s inner sums into separate
  functions did not remove the error, and limiting the number of traces
  changes which traces compile. The cause was not found.

- **Function JIT with hot-counter or first-call triggers, crashes.** With
  `opcache.jit=1215` or `1235`, PHP crashes with a segmentation fault or bus
  error in `DenseTest`, `DistributionsTest`, and `DistributionsStressTest`.
  `opcache.jit_bisect_limit` points at compiling
  `Quantile::continuous`. `BobyqaTest` also fails one test under `1235`.

## What the polyfill does

- Every `Phit\` function the polyfill defines throws code 1007 while a JIT
  compiles the process, on every PHP version, so neither the wrong beta CDF
  nor the crashes below can reach a caller of the `Phit\` functions. See
  [function-jit.md](function-jit.md#what-the-polyfill-does).
- `Svd` assigns `$shift` from a cast, which fixes the SVD case above.
- The bobyqa reference is split in two tests, and the bit-exact one skips
  itself wherever a JIT is compiling; see
  [function-jit.md](function-jit.md). That covers the `BobyqaTest` failure
  above on every version.
- CI runs PHP 8.3 without OPcache, with the optimizer alone, with the tracing
  JIT and with `opcache.jit=1205` as blocking jobs. Only the `1235` job does
  not block, because of the crash below. The tracing job blocked from
  2026-10-21: the wrong beta CDF above has never appeared on a runner, so it
  is strict until it does.
- The README recommends no JIT on PHP 8.3.

## What CI has actually seen, 2026-10-21

Both runs of commit `a2fc92b` reproduced the `1235` behaviour on the runner,
which until then had only been recorded locally. One of the twelve bobyqa
reference cases failed, the rest passed, and the process then died with
`Segmentation fault (core dumped)`, exit 139, while running a file after
`BobyqaTest` whose header never printed. The first run crashed earlier, inside
`BobyqaTest`, before any of it was reported. So both halves of this note hold:
the function JIT with hot-counter triggers computes wrong values here and also
crashes, and the crash is why that one job stays informational. The bobyqa
skip cannot help it, because the test still runs the solver.

PHP 8.4 is no longer free of this either. The same solver miscompiles under the
*tracing* JIT on PHP 8.4.26 on GitHub's x86-64 runner, materially — one case's
objective moves by 18% — which is recorded with the build string, the flags and
the bit patterns in [function-jit.md](function-jit.md). This note's claim that
"PHP 8.4 and 8.5 are not affected by the problems described here" was written
about the problems above, and 8.4 now has one of its own.

## Finding affected code

- Run test files one at a time; the order of compilation decides which traces
  exist.
- `opcache.jit_debug=0x5000` prints each trace's start and when it is
  compiled; `0x45000` adds the trace's bytecode.
- `opcache.jit_max_root_traces=N` and `opcache.jit_max_side_traces=N` limit
  how many traces compile, but a different limit also changes which traces
  compile, so a binary search on them is not reliable.
