# QuantFin

A C++17 project for single-threaded Monte Carlo pricing of path-dependent
options under geometric Brownian motion (GBM), with Matplot++ visualization.
The local project folder and executable are named `quant`; the CMake project
is named `QuantFin`.

This snapshot preserves the baseline before modularization and multithreading.
The implementation currently lives in `src/main.cpp`.

## Implemented features

- GBM path simulation using the exact exponential transition.
- Arithmetic-average Asian call and put payoffs. The average uses the N
  monitoring prices, excluding the initial price S0.
- Fixed-strike lookback call and put payoffs, including S0 in the path extrema.
- Single-threaded Monte Carlo pricing for all four options.
- Discounted standard errors calculated from sample payoff variances.
- Approximate 95% confidence intervals (`price +/- 1.96 * standard error`),
  displayed for Asian and lookback calls.
- Monte Carlo convergence demonstrations with 100, 1,000, 10,000 and 100,000 paths.
- Deterministic fixed-path payoff checks and a lookback/Asian call relation check.
- Matplot++ visualization of 100 GBM paths, saved as `gbm_paths.png`.

The existing baseline uses S0 = 100, strike = 100, risk-free rate = 0.05,
volatility = 0.20, maturity = 1 year, 252 steps and random seed 2023.
The primary pricing run uses 10,000 paths; the confidence interval run uses
100,000 paths. These settings are defined in the source.

## Requirements

- A C++17-compatible compiler.
- CMake 3.14 or newer and a supported build tool.
- Matplot++ installed with its CMake package configuration and
  `Matplot++::matplot` target available.
- A plotting backend supported by the installed Matplot++ configuration
  (typically gnuplot), with image output support and a graphical environment
  for the interactive display.

## Build

From the project root:

```sh
cmake -S . -B build
cmake --build build --parallel 2
```

If CMake cannot locate an existing Matplot++ installation, provide its prefix:

```sh
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/matplot-install-prefix
cmake --build build --parallel 2
```

The executable is `build/quant`. Build outputs and local compiled binaries are
excluded from version control.

## Run

From the project root:

```sh
./build/quant
```

The executable prints prices, standard errors, convergence results, confidence
intervals and embedded test results. It then saves `gbm_paths.png` in the current
working directory (overwriting that image if present) and opens an interactive
plot through `matplot::show()`. There is currently no headless/test-only command
line option. For build-only verification without plotting windows, run only
the build commands above.

## Validation and known limitation

The baseline was successfully configured, compiled and linked with CMake on
macOS. The executable was not run during backup verification, so the embedded
runtime tests were not re-executed. There is no standalone CTest suite configured.

The existing plotting function `polt_GBM_path` uses `i <= time.size()` when
writing the time vector. The final iteration writes past the vector's end,
which is undefined behavior and can cause a crash. This known issue is retained
in this snapshot to preserve the existing C++ implementation; a successful build
does not establish runtime correctness.

## Project layout

- `CMakeLists.txt`: builds `src/main.cpp` and links Matplot++.
- `src/main.cpp`: current simulation, payoffs, pricing, checks and visualization.
- Other files in `src/`, plus `include/`, `tests/test_pricing.cpp` and
  `benchmarks/benchmark.cpp`: existing empty placeholders, not built by CMake.
- `gbm_paths.png`: existing sample visualization preserved with the baseline.

## Reference

[QuantFin: option introduction and exotic option pricing](https://quantfin.net/posts/option-introduction/)
