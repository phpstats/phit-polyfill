# Third-party notices

phit-polyfill is written in PHP from the code of the
[phit](https://github.com/phpstats/phit) extension, which is released under
the MIT license. Parts of phit's code follow the Rust crates and the Rust
standard library routines listed below, and the polyfill follows the same
code operation for operation. Each entry gives the upstream project, its
version, its license, its copyright line as written in its license file, a
copy of that license in [`LICENSES/`](LICENSES), and the files in this package
that derive from it. Where an upstream is dual-licensed, this package uses it
under the MIT license.

## bobyqa 0.3.2

- Project: <https://github.com/pawlenartowicz/bobyqa>
- License: MIT OR Apache-2.0, used under MIT:
  [`LICENSES/bobyqa-MIT.txt`](LICENSES/bobyqa-MIT.txt)
- Copyright: Copyright (c) 2026 Paweł Lenartowicz
- Files: `src/Internal/Bobyqa.php` and `src/Internal/Bobyqa/`

The bobyqa crate derives from the BOBYQA routine of PRIMA, which its own
notices say remains under PRIMA's license:

- Project: PRIMA, <https://github.com/libprima/prima>, the version the crate
  names (v0.7.2+, commit `1d76fb88aeffb427cd17ed1e9d0d3b34f414913f`)
- License: BSD-3-Clause:
  [`LICENSES/PRIMA-BSD-3-Clause.txt`](LICENSES/PRIMA-BSD-3-Clause.txt)
- Copyright: Copyright (c) 2020--2026, Zaikun ZHANG ( https://www.zhangzk.net )
- Files: the same as bobyqa's

## argmin 0.11.0

- Project: <https://github.com/argmin-rs/argmin>
- License: MIT OR Apache-2.0, used under MIT:
  [`LICENSES/argmin-MIT.txt`](LICENSES/argmin-MIT.txt)
- Copyright: the license file has no copyright line. The crate's author is
  Stefan Kroboth.
- Files: `src/Internal/NelderMead.php`, `src/Internal/NelderMeadError.php`,
  `src/Internal/NelderMeadResult.php`, and
  `src/Internal/NelderMeadTermination.php`

## cdflib 0.4.3

- Project: <https://github.com/vigna/CDFLIB-rs>
- License: MIT: [`LICENSES/cdflib-MIT.txt`](LICENSES/cdflib-MIT.txt)
- Copyright: the license file has no copyright line. The crate's author is
  Sebastiano Vigna.
- Files: `src/Distribution/Cdflib.php`, and the calls into it in
  `src/Distribution/Beta.php`, `FisherF.php`, `Gamma.php`, `Normal.php`,
  `Special.php`, `StudentT.php`, and `Support.php`

## statrs 0.19.1

- Project: <https://github.com/statrs-dev/statrs>
- License: MIT: [`LICENSES/statrs-MIT.txt`](LICENSES/statrs-MIT.txt)
- Copyright: Copyright (c) 2016 Michael Ma
- Files: `src/Internal/Statrs.php`, and its uses in `src/Ols.php`,
  `src/Internal/Glm.php`, `src/Internal/Family/`, `src/Internal/LinearMixed.php`,
  `src/Internal/RandomIntercept.php`, and
  `src/Internal/LogisticRandomIntercept.php`

## nalgebra 0.35.0

- Project: <https://github.com/dimforge/nalgebra>
- License: Apache-2.0:
  [`LICENSES/nalgebra-Apache-2.0.txt`](LICENSES/nalgebra-Apache-2.0.txt), the
  full license text as the crate ships it
- Copyright: Copyright 2020 Sébastien Crozet
- Files: the operation order of the matrix factorizations and products in
  `src/Internal/Dense.php`, `LeastSquares.php`, `Svd.php`,
  `SymmetricEigen.php`, and `RustSort.php`, and their uses in `src/Draws.php`,
  `src/Ols.php`, `src/Internal/Glm.php`, `src/Internal/LinearMixed.php`,
  `src/Internal/RandomIntercept.php`, and
  `src/Internal/LogisticRandomIntercept.php`

## matrixmultiply 0.3.11

- Project: <https://github.com/bluss/matrixmultiply>
- License: MIT OR Apache-2.0, used under MIT:
  [`LICENSES/matrixmultiply-MIT.txt`](LICENSES/matrixmultiply-MIT.txt)
- Copyright: Copyright (c) 2016 - 2023 Ulrik Sverdrup "bluss" and
  Copyright (c) 2021 DutchGhost [constparse.rs]
- Files: the fused multiply-add kernel order in `src/Internal/Dense.php`,
  `src/Internal/Fma.php`, `src/Internal/LeastSquares.php`, and `src/Draws.php`

## rand_chacha 0.9.0 and rand 0.9.5

- Projects: <https://github.com/rust-random/rand>
- License: MIT OR Apache-2.0, used under MIT:
  [`LICENSES/rand_chacha-MIT.txt`](LICENSES/rand_chacha-MIT.txt) and
  [`LICENSES/rand-MIT.txt`](LICENSES/rand-MIT.txt)
- Copyright: Copyright 2018 Developers of the Rand project and
  Copyright (c) 2014 The Rust Project Developers
- Files: `src/Internal/ChaCha20.php`, and the seeded sampling in
  `src/Sampling.php` and `src/Distribution/Support.php`

## The Rust standard library 1.99

- Project: <https://github.com/rust-lang/rust>
- License: MIT OR Apache-2.0, used under MIT:
  [`LICENSES/rust-MIT.txt`](LICENSES/rust-MIT.txt)
- Copyright: Copyright (c) The Rust Project Contributors
- Files: `src/Internal/RustSort.php`, which follows
  `slice::sort_unstable_by`
