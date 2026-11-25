# Benchmark results — Apple M3

Measured on 25 November 2026; reused for this presentation without rerunning the full benchmark. Values below are rounded for display; the CSV files retain greater precision.

| Implementation | Threads | Median runtime (ms) | Min–max runtime (ms) | Speedup | Efficiency |
|---|---:|---:|---:|---:|---:|
| Original serial | 1 | 7593.66 | 7572.75–7762.14 | 1.0000 | 100.00% |
| Parallel engine | 1 | 7603.30 | 7569.85–7956.04 | 0.9987 | 99.87% |
| Parallel engine | 2 | 3867.73 | 3843.42–3996.32 | 1.9633 | 98.17% |
| Parallel engine | 4 | 1994.74 | 1977.10–2024.02 | 3.8068 | 95.17% |
| Parallel engine | 8 | 1307.94 | 1280.84–1376.38 | 5.8058 | 72.57% |

Speedup = serial median / parallel median. Efficiency = speedup / worker count, shown as a percentage. “Original serial” executes on the calling thread; “Parallel engine, 1” includes a worker launch, join and merge. The original serial API is the baseline for every ratio. Figures show runtime in **seconds**, while this table and the CSV use **milliseconds**.

## Recorded environment and workload

| Item | Value |
|---|---|
| CPU / architecture | Apple M3 / arm64 |
| Cores | 8 physical and logical: 4 performance + 4 efficiency |
| OS | macOS 14.7, build 23H124 |
| Compiler | Apple Clang 15.0.0, clang-1500.3.9.4 |
| CMake | 4.5.0 |
| Build / pricing language | Release / C++11 |
| Core and benchmark flags | `-O3 -DNDEBUG -std=c++11 -Wall -Wextra -Wpedantic` |
| Total paths per pricing call | 1,000,000, identical for all configurations |
| GBM time steps | 252; each path also includes its initial value |
| Initial price / strike | S0 = 100 / K = 100 |
| Volatility / risk-free rate | 0.20 / 0.05 |
| Maturity | 1 year |
| Base RNG seed | 2023; worker i uses 2023 + i |
| Payoffs | Asian call/put and fixed-strike lookback call/put |
| Observation convention | Asian average excludes S0; lookback extrema include S0 |
| Warmups / measured runs | 1 full warmup + 5 timed runs per configuration |
| Measurement order | Rotated across five rounds; each configuration occupies each position once |
| Timer | `std::chrono::steady_clock` |
| Timed work | RNG initialization, path generation, payoffs, statistics, finalization; parallel calls also include launch/join/merge |
| Excluded work | Validation, logging, CSV export and plotting |

All four prices and standard errors are computed per call. Seeds reset on each call; repetitions measure execution-time variation, not independent pricing experiments. Background activity and power/thermal state were not controlled. No build or sanitizer job ran concurrently with the recorded full benchmark.

## Data and figures

- [Summary CSV](benchmark_parallel_m3.csv): 5 configurations and workload metadata.
- [Raw observations CSV](benchmark_parallel_m3_raw.csv): 25 timings, order, prices and standard errors.
- [Execution time](figures/execution_time_vs_threads.png)
- [Speedup, including ideal linear scaling](figures/speedup_vs_threads.png)
- [Parallel efficiency](figures/parallel_efficiency_vs_threads.png)

The figures are generated from the summary CSV by the standalone C++ Matplot++ program in `benchmarks/plot_benchmark.cpp`. Blue markers represent measured parallel medians. Runtime includes a horizontal original-serial reference, not an extra parallel measurement. Dashed ideal curves are references, not observations. Connecting lines guide the eye; no interpolation experiments or confidence intervals are implied. The min–max ranges in the table are observed timing ranges, not confidence intervals.

[Experimental analysis](benchmark_analysis.md) · [Demo commands](DEMO.md) · [Validation details](parallel_validation.md)
