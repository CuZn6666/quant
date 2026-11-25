#include "pricing.hpp"
#include "gbm.hpp"
#include "payoff.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}

template <typename Exception, typename Function>
void require_throws(Function function)
{
    try { function(); }
    catch (const Exception&) { return; }
    throw std::runtime_error("Expected exception was not thrown");
}

std::array<double, 8> values(const OptionPrices& p)
{
    return {{p.asian_call, p.asian_put, p.lookback_call, p.lookback_put,
        p.asian_call_std_error, p.asian_put_std_error,
        p.lookback_call_std_error, p.lookback_put_std_error}};
}
std::array<double, 8> moments(const MonteCarloStatistics& s)
{
    return {{s.asian_call.sum, s.asian_put.sum, s.lookback_call.sum, s.lookback_put.sum,
        s.asian_call.sum_sq, s.asian_put.sum_sq, s.lookback_call.sum_sq, s.lookback_put.sum_sq}};
}
void close(double actual, double expected, double tolerance)
{
    require(std::isfinite(actual) && std::isfinite(expected) &&
        std::abs(actual - expected) <= tolerance,
        "Mismatch: actual=" + std::to_string(actual) + " expected=" + std::to_string(expected));
}
void equal(const OptionPrices& a, const OptionPrices& b)
{
    const auto av = values(a), bv = values(b);
    for (std::size_t j = 0; j < av.size(); ++j) close(av[j], bv[j], 0.0);
}

void test_single_thread()
{
    for (std::uint64_t seed : {2023ULL, 91ULL}) {
        std::mt19937_64 rng(seed);
        equal(monte_carlo_pricing(100, .2, .05, 1, 100, 17, 10003, rng),
              monte_carlo_pricing_parallel(100, .2, .05, 1, 100, 17, 10003, 1, seed));
    }
}
void test_reproducibility()
{
    for (std::size_t threads : {1, 2, 4, 8}) {
        for (std::uint64_t seed : {2023ULL, 91ULL}) {
            const auto a = monte_carlo_pricing_parallel(100, .2, .05, 1, 100, 17, 10003, threads, seed);
            const auto b = monte_carlo_pricing_parallel(100, .2, .05, 1, 100, 17, 10003, threads, seed);
            equal(a, b); // Exact field equality, not struct padding / memcmp.
        }
    }
}
void test_statistical_consistency()
{
    for (std::uint64_t seed : {2023ULL, 91ULL, 76543ULL}) {
        std::mt19937_64 rng(seed);
        const auto single = values(monte_carlo_pricing(100, .2, .05, 1, 100, 32, 30000, rng));
        for (std::size_t threads : {1, 2, 4, 8}) {
            // Different seeds avoid sharing the baseline's first worker stream.
            const auto parallel = values(monte_carlo_pricing_parallel(
                100, .2, .05, 1, 100, 32, 30000, threads, seed + 1000000));
            for (std::size_t j = 0; j < 4; ++j) {
                require(single[j + 4] > 0 && parallel[j + 4] > 0, "Positive SE required");
                const double combined_se = std::hypot(single[j + 4], parallel[j + 4]);
                close(single[j], parallel[j], 6.0 * combined_se + 1e-12);
            }
        }
    }
}

void check_partition(std::size_t paths, std::size_t threads, std::uint64_t seed)
{
    const auto actual = accumulate_monte_carlo_parallel(100, .2, .05, 1, 100, 17, paths, threads, seed);
    require(actual.n == paths, "Parallel sample count must equal requested paths");
    const std::size_t active = std::min(paths, threads);
    // Independently enumerate the assignment instead of copying quotient/remainder code.
    std::vector<std::size_t> counts(active, 0);
    for (std::size_t i = 0; i < paths; ++i) ++counts[i % active];
    std::array<std::vector<double>, 4> samples;
    for (std::size_t worker = 0; worker < active; ++worker) {
        std::mt19937_64 rng(seed + worker);
        for (std::size_t i = 0; i < counts[worker]; ++i) {
            const auto path = simulation_GBM(100, .2, .05, 17, 1, rng);
            samples[0].push_back(Asian_call_payoff(path, 100));
            samples[1].push_back(Asian_put_payoff(path, 100));
            samples[2].push_back(Lookback_call_payoff(path, 100));
            samples[3].push_back(Lookback_put_payoff(path, 100));
        }
    }
    const auto raw = moments(actual);
    const auto prices = values(finalize_monte_carlo(actual, .05, 1));
    equal(finalize_monte_carlo(actual, .05, 1), monte_carlo_pricing_parallel(
        100, .2, .05, 1, 100, 17, paths, threads, seed));
    for (std::size_t j = 0; j < 4; ++j) {
        require(samples[j].size() == paths, "Reference count");
        long double sum = 0, sum_sq = 0;
        for (double x : samples[j]) { sum += x; sum_sq += static_cast<long double>(x) * x; }
        close(raw[j], static_cast<double>(sum), 1e-11 * std::max(1.0, static_cast<double>(sum)));
        close(raw[j + 4], static_cast<double>(sum_sq), 1e-11 * std::max(1.0, static_cast<double>(sum_sq)));
        const long double mean = sum / paths;
        long double deviations = 0;
        for (double x : samples[j]) deviations += (x - mean) * (x - mean);
        close(prices[j], std::exp(-.05) * static_cast<double>(mean), 1e-10);
        close(prices[j + 4], std::exp(-.05) * std::sqrt(static_cast<double>(deviations / (paths - 1) / paths)), 1e-7);
    }
}
void test_partition_merge()
{
    for (std::size_t threads : {1, 2, 3, 4, 8}) {
        for (std::uint64_t seed : {2023ULL, 91ULL}) check_partition(10003, threads, seed);
    }
}
void test_small_samples()
{
    for (std::size_t paths : {2, 3, 7, 8, 9}) {
        for (std::size_t threads : {1, 2, 4, 8, 16}) check_partition(paths, threads, 2023);
    }
    // A huge requested thread count is capped before allocating or calculating seeds.
    check_partition(2, std::numeric_limits<std::size_t>::max(), 2023);
}
void test_zero_volatility()
{
    for (std::size_t threads : {1, 2, 4, 8}) {
        for (double s0 : {80.0, 120.0}) {
            const auto p = values(monte_carlo_pricing_parallel(s0, 0, 0, 1, 100, 17, 19, threads, 2023));
            for (std::size_t j = 0; j < 4; ++j) {
                close(p[j], j % 2 == 0 ? std::max(s0 - 100, 0.0) : std::max(100 - s0, 0.0), 0);
                close(p[j + 4], 0, 0);
            }
        }
        for (double rate : {-.05, .05}) {
            double average = 0;
            for (int step = 1; step <= 17; ++step) average += 100 * std::exp(rate * step / 17.0) / 17;
            const auto p = values(monte_carlo_pricing_parallel(100, 0, rate, 1, 100, 17, 19, threads, 2023));
            const std::array<double, 4> expected = {{std::max(average - 100, 0.0), std::max(100 - average, 0.0),
                std::max(100 * std::exp(rate) - 100, 0.0), std::max(100 - 100 * std::exp(rate), 0.0)}};
            for (std::size_t j = 0; j < 4; ++j) {
                close(p[j], std::exp(-rate) * expected[j], 1e-10);
                // Raw sum-of-squares variance has cancellation at constant noninteger payoffs.
                close(p[j + 4], 0, 1e-6);
            }
        }
        const auto expiry = values(monte_carlo_pricing_parallel(120, .2, -.05, 0, 100, 17, 19, threads, 2023));
        close(expiry[0], 20, 0);
        close(expiry[4], 0, 0);
    }
}
void test_invalid_parameters()
{
    for (std::size_t threads : {1, 2, 4, 8}) {
        for (std::size_t paths : {0, 1}) require_throws<std::invalid_argument>([&] {
            monte_carlo_pricing_parallel(100, .2, .05, 1, 100, 17, paths, threads, 2023);
        });
        // These errors originate inside workers and must be rethrown after all joins.
        require_throws<std::invalid_argument>([&] {
            monte_carlo_pricing_parallel(100, .2, .05, 1, 100, 0, 19, threads, 2023);
        });
        for (std::size_t field = 0; field < 5; ++field) {
            for (double invalid : {std::numeric_limits<double>::quiet_NaN(),
                                   std::numeric_limits<double>::infinity(),
                                   -std::numeric_limits<double>::infinity(), -1.0}) {
                if (field == 2 && invalid == -1) continue; // Negative rates are valid.
                std::array<double, 5> p = {{100, .2, .05, 1, 100}};
                p[field] = invalid;
                std::mt19937_64 rng(2023);
                const auto before = rng;
                require_throws<std::invalid_argument>([&] {
                    monte_carlo_pricing(p[0], p[1], p[2], p[3], p[4], 17, 19, rng);
                });
                require(rng == before, "Invalid input must not consume caller RNG");
                require_throws<std::invalid_argument>([&] {
                    monte_carlo_pricing_parallel(p[0], p[1], p[2], p[3], p[4], 17, 19, threads, 2023);
                });
            }
        }
    }
    require_throws<std::invalid_argument>([] { monte_carlo_pricing_parallel(100, .2, .05, 1, 100, 17, 19, 0, 2023); });
    require_throws<std::invalid_argument>([] { monte_carlo_pricing_parallel(0, .2, .05, 1, 100, 17, 19, 4, 2023); });
    const auto valid = values(monte_carlo_pricing_parallel(100, 0, -.05, 0, 0, 1, 2, 8, 2023));
    close(valid[0], 100, 0); // Zero strike, zero maturity, negative rate are supported.
    MonteCarloStatistics stats;
    stats.n = 2;
    require_throws<std::invalid_argument>([&] { finalize_monte_carlo(stats, std::numeric_limits<double>::quiet_NaN(), 1); });
    require_throws<std::invalid_argument>([&] { finalize_monte_carlo(stats, .05, -1); });
}
void test_overflow()
{
    const auto max_seed = std::numeric_limits<std::uint64_t>::max();
    require_throws<std::invalid_argument>([&] { monte_carlo_pricing_parallel(100, .2, .05, 1, 100, 1, 2, 2, max_seed); });
    check_partition(2, 8, max_seed - 1); // Boundary is based on active, not requested workers.
    check_partition(2, 1, max_seed);
    const auto max_size = std::numeric_limits<std::size_t>::max();
    std::mt19937_64 rng(2023);
    const auto before = rng;
    require_throws<std::overflow_error>([&] { simulation_GBM(100, .2, .05, max_size, 1, rng); });
    require(rng == before, "Overflow rejected before RNG use");
    require_throws<std::overflow_error>([&] { monte_carlo_pricing_parallel(100, .2, .05, 1, 100, max_size, 2, 8, 2023); });
    const auto capacity = std::vector<double>().max_size();
    require_throws<std::length_error>([&] { monte_carlo_pricing_parallel(100, .2, .05, 1, 100, capacity, 2, 8, 2023); });
    MonteCarloStatistics largest, one;
    largest.n = max_size; one.n = 1;
    require_throws<std::overflow_error>([&] { merge_statistics(largest, one); });
    require_throws<std::overflow_error>([&] { merge_statistics(one, largest); });
}
}

int main(int argc, char** argv)
{
    struct Test { const char* name; void (*run)(); };
    const Test tests[] = {{"single_thread", test_single_thread}, {"reproducibility", test_reproducibility},
        {"statistical_consistency", test_statistical_consistency}, {"partition_merge", test_partition_merge},
        {"zero_volatility", test_zero_volatility}, {"small_samples", test_small_samples},
        {"invalid_parameters", test_invalid_parameters}, {"overflow", test_overflow}};
    bool matched = false;
    try {
        require(argc <= 2, "Usage: test_parallel_pricing [test_name]");
        for (const auto& test : tests) {
            if (argc == 2 && std::string(test.name) != argv[1]) continue;
            matched = true;
            test.run();
            std::cout << "PASS: " << test.name << '\n';
        }
        require(matched, "Unknown test name");
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
