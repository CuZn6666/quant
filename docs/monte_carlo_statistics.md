# Mergeable Monte Carlo Statistics

## 数据结构

`include/pricing.hpp` 中的 `MonteCarloStatistics` 保存共享样本数 `n`，
以及 Asian Call、Asian Put、Lookback Call、Lookback Put 各自的
`PayoffStatistics { sum, sum_sq }`。所有字段默认初始化为零。

`sum` 和 `sum_sq` 分别是 **undiscounted payoff** 的和、平方和。
每条 path 同时计算四种 payoff，因此四种 option 使用相同的 `n`。

## API

- `accumulate_monte_carlo(initial_price, volatility, rate, maturity, strike,
  steps, paths, rng)`：生成 paths 并返回 raw statistics。`steps` 必须大于零；
  `paths` 可以为零或一。空批次不消耗 RNG；非空批次推进调用方提供的 RNG。
- `merge_statistics(left, right)`：返回两个批次的 `n/sum/sum_sq` 逐项相加结果，
  不修改输入。默认构造的空批次是 merge identity。样本数发生 `std::size_t`
  overflow 时抛出 `std::overflow_error`。
- `finalize_monte_carlo(statistics, rate, maturity)`：要求 `n >= 2`，否则抛出
  `std::invalid_argument`；返回现有 `OptionPrices`，不修改 raw statistics。
- 原 `monte_carlo_pricing(...)` 签名及 `M >= 2` 约束保持不变，内部依次调用
  accumulation 和 finalization。

## 统计计算

对每一种 option，使用原有公式：

```text
mean = sum / n
sample_variance = max(0, (sum_sq - n * mean * mean) / (n - 1))
discount = exp(-rate * maturity)
price = discount * mean
standard_error = discount * sqrt(sample_variance / n)
```

discounting 只在 finalization 发生。应先 merge raw statistics，再 finalize；
不要平均各批次的 standard errors，也不要对不同大小批次的 prices 简单平均。

## 使用约定

只合并来自相同模型参数、strike、time grid 和 payoff observation conventions
的批次。结构不存储这些参数，也不自动检查它们；调用方负责保持一致。
手动构造的 statistics 必须表示有效 payoff moments，空批次的所有 sums 应为零。

Phase 1 不包含 `std::thread`。后续并行调用时，每个 worker 必须独占 RNG 和
accumulator，避免重复使用相同 stream；归并时采用确定的 batch 顺序。

把同一 RNG stream 顺序拆成多个批次，不改变随机数消耗顺序；但 merge 改变了
浮点加法分组，因此结果应按 tolerance 比较，而非要求 bitwise equality。
不同 standard library 的 `std::normal_distribution` 结果也不要求完全相同。

此阶段保留原有 sum-of-squares variance 公式，包括其 cancellation 风险；
没有引入 Welford/Chan、额外依赖或新的 pricing models。

## 验证

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

新增七项 CTest，覆盖 raw statistics、split-and-merge、unequal batches、
empty/singleton batches、discounting、single-threaded regression 和 count overflow。
regression test 使用当前平台生成的 samples 和独立 two-pass variance reference，
同时检查 RNG state；没有写死跨平台 Monte Carlo 数值。
