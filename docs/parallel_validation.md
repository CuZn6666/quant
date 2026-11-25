# Parallel Monte Carlo validation

## Implementation review

- Workers construct separate `std::mt19937_64(base_seed + i)` generators. No RNG is shared. Consecutive seeds are retained; distinct seeds are not a proof of statistically independent, nonoverlapping streams.
- Active workers are `min(requested_threads, paths)`. For `P` paths and `W` active workers, `P / W` paths go to each worker and the first `P % W` workers get one extra. The assigned total is exactly `P`, including uneven and singleton batches. The quotient/remainder arithmetic does not overflow.
- Each worker returns local statistics into its own preallocated vector element; the thread index is captured by value. The same applies to exception slots. Vectors are not resized while workers use them.
- Every worker joins before errors are inspected or results read. A worker exception is saved with `exception_ptr` and rethrown on the caller after all joins. If thread construction throws, already-started workers are joined before propagation. Allocation failures before launch propagate normally. Thread-creation failure was reviewed but not fault-injected.
- Raw counts, four payoff sums and four squared sums are added in worker-index order. Count overflow is checked before addition. Finalization discounts pooled means and standard errors once. No averaging of per-worker standard errors occurs.
- The original serial API and the one-worker parallel API use the same generator seed, path generation and payoff calculations. Exact equality is tested within this build. Multiple workers change random streams and floating-point grouping.
- `accumulate_monte_carlo_parallel` exposes the existing merged raw statistics for direct count/moment verification; the pricing API wraps it with finalization. Its existing minimum of two total paths is retained. Singleton *worker* batches are supported.

## Issues addressed

1. CMake previously requested C++17, and tests used `std::string_view` plus aggregate initialization incompatible with the statistics structs' C++11 default member initializers. Core/tests/benchmarks now compile with `-std=c++11`; small constructors preserve zero initialization and brace-based initialization. Only the optional Matplot++ executable retains C++17. Test creation respects `BUILD_TESTING`.
2. `steps == SIZE_MAX` previously wrapped `steps + 1` to zero before `S[0]` was written. Shared GBM validation now rejects that overflow and lengths exceeding vector capacity before allocation/RNG use.
3. Invalid model values could previously produce NaN or nonsensical paths. GBM and batch accumulation now require finite `S0 > 0`, `sigma >= 0`, rate and `T >= 0`; accumulation also requires finite `K >= 0`. Negative rates, zero maturity, zero volatility and zero strike are valid. Finalization validates rate/maturity. Worker validation exceptions propagate through the existing join/rethrow path.
4. Previous parallel tests mostly checked a few finite results and approximate repeatability. They now verify all four prices/errors, exact fixed-configuration repeatability, raw counts/moments, statistical consistency and edge cases.
5. The benchmark previously ran all repetitions of each configuration consecutively, checked only one price, and only handled odd medians correctly. It now rotates configurations across rounds, consumes/checks all eight result fields outside timing, handles even medians, exports raw repetitions plus summary metadata, and rejects non-Release configurations. Its default workload and seed are unchanged.

Valid-input simulation, payoff conventions, variance formula and RNG strategy are unchanged. Existing uncommitted work was retained. No commits or pushes were made.

## Tests and checks

The initial Release build passed all 17 existing CTests. After changes:

| Check | Result |
|---|---|
| Release CTest | 24/24 passed |
| Debug AddressSanitizer + UndefinedBehaviorSanitizer | 24/24 passed; no sanitizer diagnostics |
| RelWithDebInfo ThreadSanitizer, parallel groups | 8/8 passed; no race diagnostics |
| Compiler warnings | None emitted with `-Wall -Wextra -Wpedantic` |
| Strict language mode | `-std=c++11` verified in core, test and benchmark compile flags |
| Benchmark CLI and CSV smoke checks | Passed: malformed/overflowing integers, invalid options, same destination rejection, non-Release rejection, 2/3 repetition medians, row counts, speedup and efficiency |

Eight separately registered parallel groups replace the former single test entry:

- `parallel_single_thread`: original serial vs one worker, two seeds, exact field equality.
- `parallel_reproducibility`: repeated identical configurations at 1/2/4/8 workers, two seeds, exact field equality.
- `parallel_statistical_consistency`: three seeds, 30,000 paths, 32 steps, all four prices at 1/2/4/8 workers. Baseline and parallel runs use different seed ranges. The tolerance is `6 * sqrt(SE_single^2 + SE_parallel^2) + 1e-12`; this is a conservative statistical regression check under an independent-sample assumption, not a proof of RNG independence or unbiasedness. Cross-thread results are not expected to match exactly.
- `parallel_partition_merge`: 10,003 paths, 1/2/3/4/8 workers, two seeds. Independently enumerates worker counts, regenerates each worker's paths, checks `n` and all raw moments against long-double sums, and checks final errors against a two-pass variance reference.
- `parallel_zero_volatility`: both call/put moneyness, constant and growing/declining deterministic paths, negative rates and zero maturity. Constant noninteger payoffs allow a small cancellation residual in SE.
- `parallel_small_samples`: 2/3/7/8/9 paths at 1/2/4/8/16 workers, plus `SIZE_MAX` requested workers capped to two paths.
- `parallel_invalid_parameters`: zero threads, 0/1 paths, zero steps, NaN/infinity, invalid signs, valid boundary inputs, worker exceptions, and unchanged caller RNG on rejected serial input.
- `parallel_overflow`: seed overflow and exact boundaries, `steps + 1` overflow, vector capacity, and both count-overflow merge orders.

The 16 existing serial/statistics tests are retained, including analytic one-step pricing, raw moments, unequal merges, empty/singleton statistics, discounting, and a two-pass serial reference.

## Reproduction

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
./build/benchmark_parallel
# Same default workload, with explicit exports:
./build/benchmark_parallel --csv docs/benchmark_parallel_m3.csv --raw-csv docs/benchmark_parallel_m3_raw.csv
```

For a fresh installation without Matplot++, add `-DQUANT_BUILD_VISUALIZATION=OFF` at configure time. The core benchmark does not require Matplot++.

Sanitizers (separate builds; do not use them for performance measurements):

```sh
cmake -S . -B build/sanitize -DCMAKE_BUILD_TYPE=Debug -DQUANT_BUILD_VISUALIZATION=OFF \
  '-DCMAKE_CXX_FLAGS=-fsanitize=address,undefined -fno-omit-frame-pointer'
cmake --build build/sanitize
ctest --test-dir build/sanitize --output-on-failure
cmake -S . -B build/tsan -DCMAKE_BUILD_TYPE=RelWithDebInfo -DQUANT_BUILD_VISUALIZATION=OFF \
  '-DCMAKE_CXX_FLAGS=-fsanitize=thread -fno-omit-frame-pointer'
cmake --build build/tsan --target test_parallel_pricing
ctest --test-dir build/tsan -R '^parallel_' --output-on-failure
```

## Remaining limits and optimization candidates

- Exact reproducibility requires the same model, seed, path count, active worker count, executable/toolchain and floating-point environment. Different thread counts change streams and grouping. `std::normal_distribution` need not produce identical values across standard-library implementations.
- Raw `sum_sq - n * mean^2` can lose precision for nearly constant payoffs; clamping negative variance to zero can conceal cancellation. Extreme finite inputs can still overflow GBM exponentials or floating-point moments. Manually supplied raw statistics are trusted; counts above `2^53` lose integer precision when converted to double. These are numerical limits, not repaired by this targeted concurrency validation.
- Worker count is capped by paths, not by hardware. Excessive requests can exhaust memory/thread resources. Resource exhaustion and OS join failures were not injected. Sanitizer passes cover exercised cases only.
- **Thread overhead and scheduling:** the parallel one-worker API adds vectors, thread creation, join and a merge. A single large-workload timing cannot isolate those costs from scheduling/frequency variation.
- **RNG and arithmetic:** every step draws a normal variate and evaluates an exponential. These are plausible compute costs; no profile has established their share of total runtime.
- **Memory/cache:** each simulated path allocates a vector and four payoff functions traverse it. Allocator contention, shared-cache pressure and locality may affect scaling; bandwidth/cache misses were not measured.
- **False sharing:** adjacent result/error slots may share cache lines, but source inspection shows the result vector is assigned only once per worker, after local accumulation. False sharing in those slots is not established as a bottleneck; inspect generated code/profile before adding padding.
- **Heterogeneous CPU:** verified hardware has four performance and four efficiency cores. Eight workers do not imply eight equal-speed cores; equal path quotas can create finish-time imbalance. Actual core placement and per-worker durations were not recorded, so this remains a possible explanation, not a measured cause.
- The next optimization experiment should profile RNG, allocation and payoff traversal, then test reusing one path buffer per worker while preserving RNG consumption and payoff conventions. Compare the full benchmark again. A stable-variance algorithm is a separate numerical improvement, not included here.
- Separate inspection found an existing out-of-bounds loop in `src/visualization.cpp` (`i <= time.size()`). That unrelated visualization routine was neither modified nor run; sanitizer results above apply to pricing tests, not the plotting executable.

## Earlier results supplied by the user (preliminary)

These were not independently reproduced under their original run conditions:

| Configuration | Median ms | Speedup | Efficiency |
|---|---:|---:|---:|
| Original single | 7556.11 | 1.000000 | 1.000000 |
| Parallel 1 | 7801.37 | 0.968562 | 0.968562 |
| Parallel 2 | 3944.13 | 1.915780 | 0.957892 |
| Parallel 4 | 2281.62 | 3.311730 | 0.827933 |
| Parallel 8 | 1471.39 | 5.135370 | 0.641921 |

The current measured run and raw samples are linked in the README. Differences between sessions do not establish a source-code speed improvement: run ordering, thermal/power state and background activity were not controlled across sessions.
