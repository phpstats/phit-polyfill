# The function JIT miscompiles code entered partway through a function

Recorded 2026-10-16 on PHP 8.5.10 (macOS arm64) and 8.5.11 (GitHub's Ubuntu
x86-64 and macOS runners, and the official `php:8.5-cli` image under x86-64
emulation). The upstream issue is
[php-src #24214](https://github.com/php/php-src/issues/24214), opened
2026-10-09 and still open. Nothing from this project has been added to it.

PHP 8.4.25 (Herd, macOS arm64) has the same bug. With the code from before the
fixes below, `opcache.jit=1235` sent the NelderMead test closure into the same
endless loop in three of four hot-counter settings and printed a `live var`
line each time.

This note used to say that the fixed code passes on PHP 8.4 under every JIT
setting. That was true of every machine it had been run on, and it is false on
GitHub's x86-64 runner: PHP 8.4.26 there miscompiles the bobyqa solver under
the tracing JIT, silently and without a `live var` line. See
[the tracing JIT on PHP 8.4](#the-tracing-jit-on-php-84) below.
PHP 8.3 has a different JIT with different problems; see
[php83-jit.md](php83-jit.md).

## Recommendation

Run the polyfill without the JIT. `opcache.jit=disable` with OPcache's
optimizer is the only configuration that has never produced a wrong number
here. Since 2026-10-30 the polyfill enforces it: every `Phit\` function it
defines throws code 1007 while a JIT compiles the process, on every PHP
version (see below). The rest of this section explains why.

**Do not fit anything through the bobyqa optimizer under the tracing JIT on
PHP 8.4.** On GitHub's x86-64 runner that configuration returns a different
minimum, not a differently rounded one: one reference case gives 3.3254 where
the crate and every other mode give 4.0388. `Phit\fit_lmer` and
`Phit\fit_glmer` default to that optimizer, so a mixed model fitted that way
could be wrong by a wide margin and say nothing about it. The tracing
JIT was recommended here until 2026-10-21.

Avoid a function JIT with hot-counter triggers too, such as `opcache.jit=1235`,
the value that setup-php's php.ini sets: under that setting the polyfill has
produced wrong singular values and an endless loop in test code, as described
below. PHP's `function` mode (`opcache.jit=function`, which is `1205`) compiles
every function when its file loads and so has no on-stack-replacement entries,
which is the mechanism behind this note; it is the least bad JIT here, and
still not one the numbers have been verified under on every build.

## The bug

The function JIT compiles a whole function at once. The compiled code can be
entered from the interpreter partway through the function:

- at a loop header, when a loop becomes hot before its function does (an
  on-stack-replacement, or OSR, entry). Only the hot-counter trigger
  (`opcache.jit` digit T = 3, as in `1235`) creates these.
- after a call that may have run in the interpreter, such as a call resolved
  by name at run time (a continuation entry).

At such an entry the JIT can restore only values that are bound to a slot of
the PHP call frame. Some values live only in a CPU register. In some of these
cases PHP prints

```
ENTRY n (block b start s) - live var v
```

to stderr while compiling, where `n` is the opcode index of the entry point.
In others it prints nothing. In both cases the compiled code reads whatever
the register holds. The maintainer's analysis in #24214 traces one cause to an
assignment such as `$n = count($a)` and the spill-conflict check added for
GH-16821.

## What happened here

CI jobs that enabled OPcache ran `opcache.jit=1235`, because setup-php sets it.

1. **An endless loop in a test closure.** The cost function in
   `tests/Unit/NelderMeadTest.php` computed `$n = count($y)` and looped
   `for ($i = 0; $i < $n; $i++)`. Entered through an OSR entry, `$n` was
   garbage, the loop ran past the array, and the test emitted warnings until
   memory ran out. This was the "hang" after `NanSafetyTest`. It depended on
   hot-counter timing and address-space layout, so it was intermittent.

2. **Wrong singular values, silently.** `Svd::reflectRows` applies a
   Householder reflection in nested loops. When an inner loop over rows
   became hot first, a value of the outer loop was not restored and one entry
   of Vᵀ came out as zero. No message was printed. `DenseTest` failed in
   about one run in ten with the default hot counters, and deterministically
   with `opcache.jit_hot_loop=10`.

3. **Harmless entries in `Cdflib::betaFrac`.** The JIT reported five
   continuation entries after the by-name calls `fdiv()`, `abs()`, and
   `is_nan()`. Those calls resolve to internal functions, which never return
   through the interpreter, so the entries were not taken and results were
   correct.

Fully qualifying every internal-function call in `src/` (`\abs()`) was tried as
a general fix and reverted (commit e902431, reverted in 1b7b95f). It removed
the by-name entries, but it changed the generated code everywhere. On arm64
it exposed three more miscompiles (in `betaFrac`, `Dense::multiply`, and
`LeastSquares::singularValues`), and on x86-64 the tracing JIT and
`opcache.jit=1205` then computed wrong incomplete-beta values in
`Cdflib::betaGrat`, which the unqualified code does not.

## What the polyfill does

- Every `Phit\` function the polyfill defines throws an exception with code
  1007 (`UNSUPPORTED`) when OPcache reports a JIT compiling the process, on
  every PHP version. `bin/generate-bootstrap` puts `Internal\Jit::refuse()`
  first in each function, so the check runs before any argument is read and
  no polyfill arithmetic runs under the JIT. It caches its verdict until
  `opcache.enable`, `opcache.enable_cli`, or `opcache.jit` changes, and costs
  about 60 ns a call. If `opcache.restrict_api` hides OPcache's status, the
  check reads the ini settings instead. Until 2026-10-30 only the mixed
  models that use an optimizer refused, and the other miscompiles in this
  note and in [php83-jit.md](php83-jit.md) had no guard.
  `Polyfill::implementation()` returns the implementation without the check,
  so the CI JIT jobs still run the suite under the JIT and record miscompiles.
  The mixed-model tests that fit a model skip themselves under a JIT, and
  `tests/Unit/JitTest.php` checks the refusal in a child PHP under the
  tracing JIT in every job.
  Calling `ini_set('opcache.jit', 'off')` satisfies the check: on PHP 8.4.25
  it returns the old mode and `jit.on` reads false afterwards. Whether it also
  stops traces that are already compiled has not been shown, so the message
  does not suggest it. `bin/jit-off-experiment` runs the bobyqa reference with
  the JIT on, switches it off, and runs it again; the `jit` CI jobs run it as
  an informational step. On 2026-10-28 it was inconclusive on 8.3, 8.4.27, and
  8.5.12, because the miscompile did not reproduce with the JIT on.
- The NelderMead test closure computes its loop bound without assigning a
  temporary and iterates the quadratic case with `foreach`.
- `Svd::reflectRows` moves its loops over rows into separate functions, so an
  OSR entry there has only parameters and local variables live. The arithmetic
  is unchanged.
- `betaFrac` calls `\fdiv()`, `\abs()`, and `\is_nan()` fully qualified and
  stores the new ratio in a fresh variable before copying the old one
  (`$next = \fdiv(...); $r0 = $r; $r = $next;`). Without the reordering the
  qualified version computed wrong values on arm64, because the JIT let the
  copy `$r0` see the call result written into `$r`.
- `tests/Unit/BobyqaTest.php` checks the reference in two tests. "it agrees
  with the crate on the status and the finiteness" asserts the status, its
  label and the finiteness flag, which survive a different search path, and
  runs in every configuration. "it reproduces the crate bit for bit" asserts
  the evaluation count, the objective, the solution, the restart count and the
  CRC of every evaluated point, which do not survive it, and skips itself when
  OPcache reports a JIT compiling the process. The ini value alone is not
  consulted, because `opcache.jit` keeps its setting when OPcache is off and
  no JIT runs. Both tests run the solver, so a JIT still executes it.
- CI runs the suite without OPcache, with the optimizer alone, with the
  tracing JIT, with `opcache.jit=1205` and with `1235`, on Linux and macOS.
  Every one of them blocks on a test failure and on any `- live var` line,
  except PHP 8.3 with `1235`, which segfaults for a reason of its own; see
  [php83-jit.md](php83-jit.md).

## Finding affected code

- `opcache.jit_debug=0x1000000` prints each function's IR before compiling it,
  headed `JIT$<function>: ; after folding`. The most recent header before an
  `ENTRY ... live var` line names the function. The printing changes timing,
  so it does not always reproduce what a plain run does.
- `opcache.jit_bisect_limit=N` compiles only the first `N` functions and names
  the first one it skips. A binary search on `N` finds a function that
  computes a wrong result without any message. The compile order varies
  between runs of the full suite, so the search needs a standalone script
  whose order is fixed.
- `opcache.jit_hot_loop` set to small values (1, 3, or 10) changes which loops
  become hot first and turns some intermittent OSR failures into
  deterministic ones.

## Verification

The suite (552 tests, and 656 with the extension's parity tests) was run on
arm64 and on x86-64 under emulation with OPcache off, the optimizer alone, the
tracing JIT, `opcache.jit=1205`, `1215`, and `1235` with default, low, and
OSR-first hot counters. Every run passed without a `live var` line, including
the `1235` runs. A sweep of 80 large inputs (regressions up to 120 columns,
formula and mixed models, prediction, distribution grids, samplers, and draw
kernels) gave results identical to the run without OPcache under all of those
settings on arm64, and under the tracing JIT and `1205` on x86-64. Passing
`1235` runs do not rule out other OSR miscompiles, which depend on which loops
become hot first.

## The tracing JIT on PHP 8.4

Recorded 2026-10-21 from two runs of the same commit in CI, and reported here
so it can be taken upstream.

- Build: `PHP 8.4.26 (cli) (built: Sep 23 2026 11:52:51) (NTS)`, Zend Engine
  v4.4.26, Zend OPcache v8.4.26, `ubuntu-latest` x86-64, installed by
  `shivammathur/setup-php@v2`.
- Flags: `-d opcache.enable_cli=1 -d opcache.file_update_protection=0
  -d opcache.jit=tracing -d opcache.jit_buffer_size=64M`, over setup-php's
  php.ini (`opcache.jit=1235`, `opcache.jit_buffer_size=256M`,
  `opcache.memory_consumption=128`, `opcache.optimization_level=0x7FFEBFFF`).
- Commit `a2fc92b`, which changes no arithmetic: the diff from the last green
  commit renames PHP parameters, regenerates `bootstrap.php` from them, adds a
  parity test and adds a lint job.

`Tests\Unit\BobyqaTest > it reproduces the crate bit for bit` failed ten of
its twelve reference cases on the first run and nine on the second. The nine
are the same cases with the same numbers both times; the tenth,
`rosenbrock-2-default`, failed once and passed once. `sphere-default`, which
runs first, and `all-infinite`, which returns before doing real arithmetic,
passed both times. No `live var` line was printed, and every other PHP 8.4 job
reproduced the crate exactly: no OPcache, the optimizer alone, `1235` and
`1205` on Linux, and the optimizer alone and `1205` on macOS.

The failures are a divergent search, not one wrong digit. Six cases differ
first in the evaluation count: 155 against the expected 127, then 171/168,
53/52, 83/89, 22/21 and 396/477. Three differ first in the objective's bit
pattern, expected against actual:

```
401027c5914b4d5a  400a9a70c4c3ab06
3ff9822f48da29d0  3ff98649ce928846
bfd48fb21e668740  bfd48fb21e6625d5
```

One differs first in the solution vector:

```
expected  bfe0c1bad0042674  3fbe4530710fc989
actual    bfe0c1bacf1aa622  3fbe453070ded7b2
```

In decimal, with the relative difference:

| expected | actual | relative |
|---|---|---|
| 4.0388395979128635 | 3.325410401555362 | 1.8e-01 |
| 1.5942833753384669 | 1.595285231517225 | 6.3e-04 |
| -0.3212704941635032 | -0.3212704941621188 | 4.3e-12 |

**The divergence is material, not a rounding difference.** The first case ends
18% away from the crate's minimum, and the second 0.06% away, which is far
outside double precision. The third, and both components of the solution vector
(3.2e-09 and 3.8e-10), are immaterial on their own; they are the cases where
the wrong arithmetic had not yet moved the search. So the skip in the test
suite is not about strictness. It marks the configurations where this
optimizer's answers cannot be trusted at all, which is why the recommendation
above tells users to stay off the tracing JIT on PHP 8.4 rather than to expect
a few wrong bits.

It does not reproduce off that runner. Seven full-suite runs at the same
commit with the same flags passed, with all twelve cases exact: macOS arm64 on
PHP 8.4.25; `php:8.4-cli` 8.4.26 under `--platform linux/amd64`, plain and
with `jit_hot_loop`, `jit_hot_func`, `jit_hot_return` and `jit_hot_side_exit`
all forced to 1; and `ubuntu:24.04` with `ppa:ondrej/php` PHP 8.4.26, the same
package line as the runner, plain and with `jit_hot_side_exit=1`,
`jit_max_root_traces=40`, and `jit_hot_loop=8` with `jit_hot_func=16`.
