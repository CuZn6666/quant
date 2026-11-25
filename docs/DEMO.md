# Baumgarten demo — Monte Carlo pricing and parallel scaling

Run commands from the repository root. The short presentation uses the existing full-workload results and does not repeat the multi-minute benchmark.

## 1. Build once

Matplot++ and gnuplot must be installed (`brew install matplotplusplus gnuplot` if missing). PNG export uses the Matplot++ gnuplot backend; no interactive plot window is needed.

The plotter sets the `C` locale and `pngcairo` terminal only for its own process and child gnuplot processes, avoiding unavailable-locale initialization errors on macOS. Its single-axes backend omits the unnecessary multiplot wrapper to avoid gnuplot 6.0.6 inline-data warnings.

```sh
cd ~/Desktop/quant
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DQUANT_BUILD_VISUALIZATION=OFF -DQUANT_BUILD_BENCHMARK_PLOTS=ON
cmake --build build
```

The pricing library, tests, benchmark and `pricing_demo` compile as **C++11**. The installed Matplot++ package requires **C++17**, confined to the optional `plot_benchmark` executable. To build the pricing project without plotting dependencies, set both visualization and benchmark-plot options to `OFF`. The existing interactive `quant` visualization is not part of this demo.

## 2. Show pricing results

```sh
./build/pricing_demo
```

This quick demo prints all four prices and their Monte Carlo standard errors for the original serial API and parallel 1/2/4/8 workers. Every configuration uses 10,000 paths, 252 steps, S0 = K = 100, sigma = 0.20, r = 0.05, T = 1 and seed 2023. It is a pricing demonstration, not a runtime measurement. The one-worker result matches serial; other thread counts use different random streams and need not match exactly.

## 3. Show the 24 tests

```sh
ctest --test-dir build --output-on-failure
```

Coverage includes fixed-seed reproducibility, statistical consistency, path counts and moment merging, uneven workloads, deterministic cases, small samples, invalid inputs and overflow.

## 4. Present the recorded benchmark and regenerate figures

```sh
cat docs/benchmark_parallel_m3.csv
./build/plot_benchmark docs/benchmark_parallel_m3.csv docs/figures
```

Present the [complete results and environment table](benchmark_results.md), then the [short experimental analysis](benchmark_analysis.md). The plotter reads and validates CSV fields, workload consistency and baseline ratios; it runs no Monte Carlo simulations. It writes three 1600 × 1000 PNGs using one consistent style and overwrites the same image names on regeneration:

- `docs/figures/execution_time_vs_threads.png`
- `docs/figures/speedup_vs_threads.png`
- `docs/figures/parallel_efficiency_vs_threads.png`

Talking points: 7.59 s serial → 1.31 s with eight workers; 5.81× speedup; 72.57% efficiency. The single-worker difference is within observed timing variation. The causes of scaling losses have not been isolated by profiling.

## 5. Optional live benchmark

A quick execution demonstration uses a smaller workload; its timing must not be compared to the recorded full-workload table:

```sh
./build/benchmark_parallel --paths 10003 --steps 32 --repeats 3 \
  --csv build/demo_benchmark.csv --raw-csv build/demo_benchmark_raw.csv
```

Only if a new full experiment is wanted, run the unchanged default workload (1,000,000 paths, 252 steps, five timed runs plus one warmup per configuration):

```sh
./build/benchmark_parallel \
  --csv build/benchmark_new.csv --raw-csv build/benchmark_new_raw.csv
```

Keep new measurements separate from the archived `docs/` CSV files, record the new environment and avoid concurrent builds/tests while timing.
