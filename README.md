[QuantFin: option introduction and exotic option pricing](https://quantfin.net/posts/option-introduction/)

## Multithreaded Monte Carlo Performance

The C++11 pricing engine simulates GBM paths for Asian and fixed-strike lookback calls/puts. Each `std::thread` owns an `mt19937_64` seeded with `base_seed + worker_index`, a path quota and local payoff statistics. Remainder paths go to the first workers. All workers join before counts, sums and squared sums are merged in worker order; discounting and standard errors are computed from the pooled samples. Worker exceptions are rethrown after joining. No RNG or running accumulator is shared.

Build and test:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

To build without the optional Matplot++ visualization dependency, add `-DQUANT_BUILD_VISUALIZATION=OFF` to the configure command. Only the visualization executable uses C++17; the core, tests and benchmarks use strict C++11. The validated Release build passed **24/24 tests**, with no compiler warnings. AddressSanitizer/UndefinedBehaviorSanitizer passed 24/24; ThreadSanitizer passed all 8 parallel groups.

Run the benchmark (Release required):

```sh
./build/benchmark_parallel
# Same workload, exporting summary and every timed observation:
./build/benchmark_parallel --csv docs/benchmark_parallel_m3.csv --raw-csv docs/benchmark_parallel_m3_raw.csv
# Optional shorter smoke run (not comparable to the full workload):
./build/benchmark_parallel --paths 10003 --steps 32 --repeats 3
```

Default workload: **1,000,000 total paths per configuration, 252 time steps, S0 = K = 100, sigma = 0.20, r = 0.05, T = 1 year, base seed = 2023**. All four prices and errors are computed in each call. Each configuration receives one full warmup and five timed repetitions, resetting seeds each call. Configuration order rotates across rounds; with five repetitions each configuration occupies every position once. Median time includes path allocation, RNG initialization, pricing and, for the parallel API, worker launch/join/merge. Validation, logging, CSV output and visualization are outside timing. CSV files include workload metadata; raw output also includes run order and all prices/errors. Existing CSV destinations are overwritten.

**Measured on this Mac, 2026-11-25:** Apple M3, 8 physical/logical cores (4 performance + 4 efficiency), macOS 14.7 (23H124), Apple Clang 15.0.0 (`clang-1500.3.9.4`), arm64, CMake 4.5.0. Core/benchmark flags: `-O3 -DNDEBUG -std=c++11 -Wall -Wextra -Wpedantic`. Hardware was verified with `sysctl`; OS/compiler/build flags were inspected locally. These are measurements from one session, not guaranteed performance across machines or power/thermal conditions.

| Configuration | Median ms | Min–max ms (5 runs) | Speedup | Efficiency |
|---|---:|---:|---:|---:|
| Original single | 7593.66 | 7572.75–7762.14 | 1.0000 | 100.00% |
| Parallel 1 | 7603.30 | 7569.85–7956.04 | 0.9987 | 99.87% |
| Parallel 2 | 3867.73 | 3843.42–3996.32 | 1.9633 | 98.17% |
| Parallel 4 | 1994.74 | 1977.10–2024.02 | 3.8068 | 95.17% |
| Parallel 8 | 1307.94 | 1280.84–1376.38 | 5.8058 | 72.57% |

[Summary CSV](docs/benchmark_parallel_m3.csv) · [25 raw timed observations](docs/benchmark_parallel_m3_raw.csv). The one-worker median was about 0.13% slower than serial, smaller than the observed within-configuration timing ranges; this does not isolate thread overhead. Four workers achieved 3.81× and eight achieved 5.81×, so scaling improved while efficiency fell to 72.57% at eight. Thread creation/scheduling, RNG work, path allocation, cache/memory effects and unequal core speeds are possible contributors; their individual costs were not measured. Adjacent result slots may share cache lines, but each is assigned once per worker in the source, so false sharing is not an established bottleneck. No sanitizer or build job ran concurrently with this benchmark; other background activity and power/thermal state were not controlled.

`speedup(p) = median(original single) / median(parallel p)`; `efficiency(p) = speedup(p) / p`. Baseline speedup/efficiency are 1. Efficiency is a scaling ratio, not a measurement of CPU utilization. “Original single” calls the existing serial pricing API on the calling thread. “Parallel 1” adds a worker, result/error storage, join and merge; it is useful for assessing wrapper costs and is not the primary baseline.

For fixed parameters, seed, path count, active thread count and executable/toolchain, repeated results are exactly reproducible in the tests. Changing thread count changes random streams and floating-point grouping; cross-thread-count or cross-standard-library bitwise reproducibility is not promised. Independent-seed comparisons use six combined standard errors, not equality. Adjacent seeds retain the existing RNG strategy and do not formally guarantee independent streams. The raw sum-of-squares variance formula still has cancellation risk for nearly constant payoffs, and extreme finite inputs can overflow floating-point calculations.

See the [implementation review, tests, sanitizer commands and limitations](docs/parallel_validation.md), including the earlier user-supplied **preliminary** results. Next, profile RNG, allocation and payoff traversal; then measure an experiment that reuses a path buffer per worker while preserving random-number consumption. No bottleneck has yet been isolated by profiling.
