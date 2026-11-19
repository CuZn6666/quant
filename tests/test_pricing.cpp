#include "gbm.hpp"
#include "payoff.hpp"
#include "pricing.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void require_close(double actual, double expected, double tolerance,
                   const std::string& message)
{
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::fabs(actual - expected) > tolerance) {
        std::ostringstream details;
        details << message << ": expected " << expected << ", got " << actual
                << ", tolerance " << tolerance;
        throw std::runtime_error(details.str());
    }
}

template <typename Callable>
void require_invalid_argument(Callable operation, const std::string& message)
{
    try {
        operation();
    } catch (const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error(message + ": expected std::invalid_argument");
}

void test_asian_call()
{
    require_close(Asian_call_payoff({100, 105, 110, 95, 120}, 100), 7.5, 1e-12,
                  "Asian call fixed path");
    require_close(Asian_call_payoff({1000, 90, 110}, 100), 0.0, 1e-12,
                  "Asian call excludes S0");
    require_close(Asian_call_payoff({100, 80, 90}, 100), 0.0, 1e-12,
                  "Asian call out of the money");
    require_close(Asian_call_payoff({100, 120}, 100), 20.0, 1e-12,
                  "Asian call single observation");
}

void test_asian_put()
{
    require_close(Asian_put_payoff({100, 95, 90, 105, 80}, 100), 7.5, 1e-12,
                  "Asian put fixed path");
    require_close(Asian_put_payoff({1, 90, 110}, 100), 0.0, 1e-12,
                  "Asian put excludes S0");
    require_close(Asian_put_payoff({100, 105, 110, 95, 120}, 100), 0.0, 1e-12,
                  "Asian put out of the money");
    require_close(Asian_put_payoff({100, 80}, 100), 20.0, 1e-12,
                  "Asian put single observation");
}

void test_lookback_call()
{
    require_close(Lookback_call_payoff({100, 105, 110, 95, 120}, 100), 20.0, 1e-12,
                  "Lookback call fixed path");
    require_close(Lookback_call_payoff({130, 90, 110}, 100), 30.0, 1e-12,
                  "Lookback call includes S0");
    require_close(Lookback_call_payoff({80, 90, 100}, 100), 0.0, 1e-12,
                  "Lookback call at the money");
    require_close(Lookback_call_payoff({80, 90}, 100), 0.0, 1e-12,
                  "Lookback call out of the money");
    require_close(Lookback_call_payoff({120}, 100), 20.0, 1e-12,
                  "Lookback call singleton path");
}

void test_lookback_put()
{
    require_close(Lookback_put_payoff({100, 105, 110, 95, 120}, 100), 5.0, 1e-12,
                  "Lookback put fixed path");
    require_close(Lookback_put_payoff({70, 90, 110}, 100), 30.0, 1e-12,
                  "Lookback put includes S0");
    require_close(Lookback_put_payoff({100, 110, 120}, 100), 0.0, 1e-12,
                  "Lookback put at the money");
    require_close(Lookback_put_payoff({110, 120}, 100), 0.0, 1e-12,
                  "Lookback put out of the money");
    require_close(Lookback_put_payoff({80}, 100), 20.0, 1e-12,
                  "Lookback put singleton path");
}

void test_gbm_properties()
{
    const std::size_t steps = 32;
    std::mt19937_64 rng(2023);
    std::mt19937_64 repeated_rng(2023);
    const auto initial_rng = rng;
    const auto path = simulation_GBM(100.0, 0.2, 0.05, steps, 1.0, rng);
    require(path.size() == steps + 1, "GBM path length includes S0");
    require_close(path.front(), 100.0, 1e-12, "GBM initial price");
    for (double price : path) {
        require(std::isfinite(price) && price > 0.0, "GBM prices are finite and positive");
    }
    require(rng != initial_rng, "GBM advances the caller's RNG");
    require(path == simulation_GBM(100.0, 0.2, 0.05, steps, 1.0, repeated_rng),
            "GBM fixed seed is reproducible within this implementation");

    const auto deterministic = simulation_GBM(100.0, 0.0, 0.05, steps, 2.0, rng);
    for (std::size_t index = 0; index < deterministic.size(); ++index) {
        const double time = 2.0 * static_cast<double>(index) / steps;
        require_close(deterministic[index], 100.0 * std::exp(0.05 * time), 1e-10,
                      "Zero-volatility GBM follows deterministic growth");
    }
    const auto single_step = simulation_GBM(100.0, 0.0, 0.05, 1, 2.0, rng);
    require(single_step.size() == 2, "GBM accepts one time step");
    require_close(single_step.back(), 100.0 * std::exp(0.1), 1e-10,
                  "Single-step GBM maturity");
}

void test_monte_carlo_properties()
{
    const std::size_t steps = 32;
    std::mt19937_64 rng(2023);
    const auto initial_rng = rng;
    const auto prices = monte_carlo_pricing(100.0, 0.2, 0.05, 1.0, 100.0,
                                           steps, 10000, rng);
    const std::array<double, 8> values = {
        prices.asian_call, prices.asian_put, prices.lookback_call, prices.lookback_put,
        prices.asian_call_std_error, prices.asian_put_std_error,
        prices.lookback_call_std_error, prices.lookback_put_std_error
    };
    for (double value : values) {
        require(std::isfinite(value) && value > 0.0,
                "Stochastic prices and standard errors are finite and positive");
    }
    require(rng != initial_rng, "Monte Carlo advances the caller's RNG");
    require(prices.lookback_call >= prices.asian_call,
            "Lookback call dominates Asian call on shared paths");
    require(prices.lookback_put >= prices.asian_put,
            "Lookback put dominates Asian put on shared paths");

    double expected_average = 0.0;
    for (std::size_t observation = 1; observation <= steps; ++observation) {
        expected_average += 100.0 * std::exp(0.05 * observation / steps) / steps;
    }
    const double expected_parity = std::exp(-0.05) * (expected_average - 100.0);
    const double tolerance = 6.0 * (prices.asian_call_std_error +
                                    prices.asian_put_std_error) + 0.01;
    require_close(prices.asian_call - prices.asian_put, expected_parity, tolerance,
                  "Asian put-call parity agrees with the analytic expected average");
}

void test_monte_carlo_zero_volatility()
{
    for (double initial_price : {80.0, 120.0}) {
        std::mt19937_64 rng(2023);
        const auto constant = monte_carlo_pricing(initial_price, 0.0, 0.0, 1.0,
                                                 100.0, 4, 8, rng);
        const double expected_call = initial_price > 100.0 ? 20.0 : 0.0;
        const double expected_put = initial_price < 100.0 ? 20.0 : 0.0;
        require_close(constant.asian_call, expected_call, 1e-12, "Constant Asian call");
        require_close(constant.asian_put, expected_put, 1e-12, "Constant Asian put");
        require_close(constant.lookback_call, expected_call, 1e-12, "Constant lookback call");
        require_close(constant.lookback_put, expected_put, 1e-12, "Constant lookback put");
        for (double error : {constant.asian_call_std_error, constant.asian_put_std_error,
                             constant.lookback_call_std_error, constant.lookback_put_std_error}) {
            require_close(error, 0.0, 1e-12, "Constant payoff has zero standard error");
        }
    }

    std::mt19937_64 rng(2023);
    const auto discounted = monte_carlo_pricing(100.0, 0.0, 0.05, 1.0, 90.0, 1, 2, rng);
    const double expected_call = 100.0 - 90.0 * std::exp(-0.05);
    require_close(discounted.asian_call, expected_call, 1e-10, "Asian call discounting");
    require_close(discounted.lookback_call, expected_call, 1e-10, "Lookback call discounting");
    require_close(discounted.asian_put, 0.0, 1e-12, "Deterministic Asian put");
    require_close(discounted.lookback_put, 0.0, 1e-12, "Deterministic lookback put");
}

void test_monte_carlo_one_step()
{
    std::mt19937_64 rng(8675309);
    const auto prices = monte_carlo_pricing(100.0, 0.2, 0.05, 1.0, 100.0, 1, 20000, rng);
    const auto normal_cdf = [](double value) {
        return 0.5 * std::erfc(-value / std::sqrt(2.0));
    };
    const double d_plus = (0.05 + 0.5 * 0.2 * 0.2) / 0.2;
    const double d_minus = d_plus - 0.2;
    const double expected_call = 100.0 * normal_cdf(d_plus) -
                                 100.0 * std::exp(-0.05) * normal_cdf(d_minus);
    const double expected_put = 100.0 * std::exp(-0.05) * normal_cdf(-d_minus) -
                                100.0 * normal_cdf(-d_plus);
    require(prices.asian_call_std_error > 0.0 && prices.asian_call_std_error < 0.5,
            "One-step call standard error is sensible");
    require(prices.asian_put_std_error > 0.0 && prices.asian_put_std_error < 0.5,
            "One-step put standard error is sensible");
    require_close(prices.asian_call, expected_call,
                  6.0 * prices.asian_call_std_error + 0.01, "One-step Black-Scholes call");
    require_close(prices.asian_put, expected_put,
                  6.0 * prices.asian_put_std_error + 0.01, "One-step Black-Scholes put");
    require_close(prices.lookback_call, prices.asian_call, 1e-12,
                  "One-step ATM lookback and Asian calls coincide");
    require_close(prices.lookback_put, prices.asian_put, 1e-12,
                  "One-step ATM lookback and Asian puts coincide");
    require_close(prices.lookback_call_std_error, prices.asian_call_std_error, 1e-12,
                  "One-step call standard errors coincide");
    require_close(prices.lookback_put_std_error, prices.asian_put_std_error, 1e-12,
                  "One-step put standard errors coincide");
}

void test_invalid_inputs()
{
    for (const auto& path : {std::vector<double>{}, std::vector<double>{100.0}}) {
        require_invalid_argument([&] { Asian_call_payoff(path, 100.0); },
                                 "Asian call rejects insufficient observations");
        require_invalid_argument([&] { Asian_put_payoff(path, 100.0); },
                                 "Asian put rejects insufficient observations");
    }
    require_invalid_argument([] { Lookback_call_payoff({}, 100.0); },
                             "Lookback call rejects empty paths");
    require_invalid_argument([] { Lookback_put_payoff({}, 100.0); },
                             "Lookback put rejects empty paths");
    std::mt19937_64 rng(2023);
    require_invalid_argument([&] { simulation_GBM(100.0, 0.2, 0.05, 0, 1.0, rng); },
                             "GBM rejects zero steps");
    for (std::size_t samples : {0, 1}) {
        require_invalid_argument([&] {
            monte_carlo_pricing(100.0, 0.2, 0.05, 1.0, 100.0, 32, samples, rng);
        }, "Monte Carlo rejects fewer than two samples");
    }
    require_invalid_argument([&] {
        monte_carlo_pricing(100.0, 0.2, 0.05, 1.0, 100.0, 0, 2, rng);
    }, "Monte Carlo propagates invalid GBM steps");
}

std::array<double, 8> statistics_values(const MonteCarloStatistics& statistics)
{
    return {statistics.asian_call.sum, statistics.asian_put.sum,
            statistics.lookback_call.sum, statistics.lookback_put.sum,
            statistics.asian_call.sum_sq, statistics.asian_put.sum_sq,
            statistics.lookback_call.sum_sq, statistics.lookback_put.sum_sq};
}

std::array<double, 8> price_values(const OptionPrices& prices)
{
    return {prices.asian_call, prices.asian_put, prices.lookback_call, prices.lookback_put,
            prices.asian_call_std_error, prices.asian_put_std_error,
            prices.lookback_call_std_error, prices.lookback_put_std_error};
}

void require_values_close(const std::array<double, 8>& actual,
                          const std::array<double, 8>& expected,
                          double tolerance, const std::string& message)
{
    for (std::size_t index = 0; index < actual.size(); ++index) {
        require_close(actual[index], expected[index],
                      tolerance * std::max(1.0, std::fabs(expected[index])),
                      message + " field " + std::to_string(index));
    }
}

void require_statistics_close(const MonteCarloStatistics& actual,
                              const MonteCarloStatistics& expected,
                              double tolerance, const std::string& message)
{
    require(actual.n == expected.n, message + " sample count");
    require_values_close(statistics_values(actual), statistics_values(expected),
                         tolerance, message);
}

void test_statistics_raw()
{
    std::mt19937_64 rng(2023);
    const auto rising = accumulate_monte_carlo(100.0, 0.0, std::log(4.0), 1.0,
                                              250.0, 2, 3, rng);
    const MonteCarloStatistics expected_rising{
        3, {150.0, 7500.0}, {0.0, 0.0}, {450.0, 67500.0}, {450.0, 67500.0}
    };
    require_statistics_close(rising, expected_rising, 1e-12,
                             "Raw statistics for paths {100, 200, 400}");

    const auto falling = accumulate_monte_carlo(100.0, 0.0, -std::log(4.0), 1.0,
                                               60.0, 2, 3, rng);
    const MonteCarloStatistics expected_falling{
        3, {0.0, 0.0}, {67.5, 1518.75}, {120.0, 4800.0}, {105.0, 3675.0}
    };
    require_statistics_close(falling, expected_falling, 1e-12,
                             "Raw statistics for paths {100, 50, 25}");
}

void test_statistics_split_merge()
{
    std::mt19937_64 whole_rng(2023);
    std::mt19937_64 split_rng(2023);
    const auto whole = accumulate_monte_carlo(100.0, 0.2, 0.05, 1.0, 100.0,
                                             17, 97, whole_rng);
    const auto first = accumulate_monte_carlo(100.0, 0.2, 0.05, 1.0, 100.0,
                                             17, 30, split_rng);
    const auto second = accumulate_monte_carlo(100.0, 0.2, 0.05, 1.0, 100.0,
                                              17, 67, split_rng);
    const auto merged = merge_statistics(first, second);
    require_statistics_close(merged, whole, 1e-12, "Split and merged raw statistics");
    require(whole_rng == split_rng, "Splitting preserves RNG consumption");
    require_values_close(price_values(finalize_monte_carlo(merged, 0.05, 1.0)),
                         price_values(finalize_monte_carlo(whole, 0.05, 1.0)),
                         1e-10, "Split and merged prices and standard errors");
}

void test_statistics_unequal_batches()
{
    const MonteCarloStatistics first{
        1, {10.0, 100.0}, {0.0, 0.0}, {40.0, 1600.0}, {20.0, 400.0}
    };
    const MonteCarloStatistics second{
        3, {0.0, 0.0}, {30.0, 300.0}, {60.0, 1200.0}, {120.0, 4800.0}
    };
    const auto first_before = first;
    const auto second_before = second;
    const auto merged = merge_statistics(first, second);
    const MonteCarloStatistics expected{
        4, {10.0, 100.0}, {30.0, 300.0}, {100.0, 2800.0}, {140.0, 5200.0}
    };
    require_statistics_close(merged, expected, 0.0, "Unequal batch counts and moments");
    require_statistics_close(merge_statistics(second, first), expected, 0.0,
                             "Reversed merge");
    require_statistics_close(first, first_before, 0.0, "Merge preserves first input");
    require_statistics_close(second, second_before, 0.0, "Merge preserves second input");
    require_values_close(price_values(finalize_monte_carlo(merged, 0.0, 1.0)),
                         {2.5, 7.5, 25.0, 35.0, 2.5, 2.5, 5.0, 5.0}, 1e-12,
                         "Weighted means and pooled sample variances");
}

void test_statistics_empty_singleton()
{
    std::mt19937_64 rng(2023);
    const auto initial_rng = rng;
    const auto empty = accumulate_monte_carlo(120.0, 0.0, 0.0, 1.0, 100.0, 8, 0, rng);
    const MonteCarloStatistics default_statistics;
    require_statistics_close(empty, default_statistics, 0.0, "Empty batch is all zero");
    require(rng == initial_rng, "Empty batch does not consume RNG");

    const auto singleton = accumulate_monte_carlo(120.0, 0.0, 0.0, 1.0, 100.0, 8, 1, rng);
    const MonteCarloStatistics expected_singleton{
        1, {20.0, 400.0}, {0.0, 0.0}, {20.0, 400.0}, {0.0, 0.0}
    };
    require_statistics_close(singleton, expected_singleton, 0.0, "Singleton raw statistics");
    require_statistics_close(merge_statistics(empty, singleton), singleton, 0.0,
                             "Left empty merge identity");
    require_statistics_close(merge_statistics(singleton, empty), singleton, 0.0,
                             "Right empty merge identity");
    require_statistics_close(merge_statistics(empty, empty), empty, 0.0,
                             "Merging empty batches");
    require_invalid_argument([&] { finalize_monte_carlo(empty, 0.0, 1.0); },
                             "Finalization rejects empty statistics");
    require_invalid_argument([&] { finalize_monte_carlo(singleton, 0.0, 1.0); },
                             "Finalization rejects singleton statistics");
    const auto another = accumulate_monte_carlo(120.0, 0.0, 0.0, 1.0, 100.0, 8, 1, rng);
    require_values_close(
        price_values(finalize_monte_carlo(merge_statistics(singleton, another), 0.0, 1.0)),
        {20.0, 0.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 1e-12,
        "Two singleton batches can be finalized");
    for (std::size_t paths : {0, 1}) {
        require_invalid_argument([&] {
            accumulate_monte_carlo(100.0, 0.2, 0.05, 1.0, 100.0, 0, paths, rng);
        }, "Accumulation rejects zero time steps even for small batches");
    }
}

void test_statistics_discounting()
{
    const MonteCarloStatistics statistics{
        3, {6.0, 14.0}, {12.0, 56.0}, {18.0, 126.0}, {24.0, 224.0}
    };
    const auto before = statistics;
    const double discount = std::exp(-0.05 * 2.0);
    const double base_error = 1.0 / std::sqrt(3.0);
    const std::array<double, 8> expected{
        2.0, 4.0, 6.0, 8.0, base_error, 2.0 * base_error,
        3.0 * base_error, 4.0 * base_error
    };
    auto discounted_expected = expected;
    for (double& value : discounted_expected) {
        value *= discount;
    }
    require_values_close(price_values(finalize_monte_carlo(statistics, 0.05, 2.0)),
                         discounted_expected, 1e-12, "Prices and SE discounted exactly once");
    require_values_close(price_values(finalize_monte_carlo(statistics, 0.0, 2.0)),
                         expected, 1e-12, "Zero-rate finalization");
    require_statistics_close(statistics, before, 0.0, "Finalization preserves raw statistics");
}

void test_single_threaded_regression()
{
    for (unsigned long long seed : {2023ULL, 8675309ULL}) {
        for (std::size_t paths : {2, 257}) {
            std::mt19937_64 reference_rng(seed);
            std::mt19937_64 pricing_rng(seed);
            std::array<std::vector<double>, 4> payoffs;
            for (std::size_t index = 0; index < paths; ++index) {
                const auto path = simulation_GBM(100.0, 0.2, 0.05, 17, 1.0, reference_rng);
                payoffs[0].push_back(Asian_call_payoff(path, 100.0));
                payoffs[1].push_back(Asian_put_payoff(path, 100.0));
                payoffs[2].push_back(Lookback_call_payoff(path, 100.0));
                payoffs[3].push_back(Lookback_put_payoff(path, 100.0));
            }
            std::array<double, 8> expected{};
            for (std::size_t option = 0; option < payoffs.size(); ++option) {
                const auto& samples = payoffs[option];
                const double mean = std::accumulate(samples.begin(), samples.end(), 0.0) / paths;
                double squared_deviations = 0.0;
                for (double payoff : samples) {
                    const double deviation = payoff - mean;
                    squared_deviations += deviation * deviation;
                }
                const double sample_variance = squared_deviations / (paths - 1);
                expected[option] = std::exp(-0.05) * mean;
                expected[option + 4] = std::exp(-0.05) * std::sqrt(sample_variance / paths);
            }
            const auto actual = monte_carlo_pricing(100.0, 0.2, 0.05, 1.0, 100.0,
                                                    17, paths, pricing_rng);
            require_values_close(price_values(actual), expected, 1e-10,
                                 "Legacy pricing agrees with a two-pass sample reference");
            require(pricing_rng == reference_rng, "Legacy pricing preserves RNG consumption");
        }
    }
    std::mt19937_64 rng(2023);
    const auto initial_rng = rng;
    for (std::size_t paths : {0, 1}) {
        require_invalid_argument([&] {
            monte_carlo_pricing(100.0, 0.2, 0.05, 1.0, 100.0, 17, paths, rng);
        }, "Legacy pricing still rejects fewer than two samples");
        require(rng == initial_rng, "Rejected legacy calls do not consume RNG");
    }
}

void test_statistics_count_overflow()
{
    MonteCarloStatistics largest{};
    largest.n = std::numeric_limits<std::size_t>::max();
    require_statistics_close(merge_statistics(largest, MonteCarloStatistics{}),
                             largest, 0.0, "Maximum count can merge with empty batch");
    MonteCarloStatistics singleton{};
    singleton.n = 1;
    try {
        merge_statistics(largest, singleton);
    } catch (const std::overflow_error&) {
        return;
    }
    throw std::runtime_error("Merge must reject sample-count overflow");
}

struct TestCase {
    std::string_view name;
    void (*run)();
};

}

int main(int argc, char* argv[])
{
    const std::array<TestCase, 16> tests = {{
        {"asian_call", test_asian_call},
        {"asian_put", test_asian_put},
        {"lookback_call", test_lookback_call},
        {"lookback_put", test_lookback_put},
        {"gbm_properties", test_gbm_properties},
        {"monte_carlo_properties", test_monte_carlo_properties},
        {"monte_carlo_zero_volatility", test_monte_carlo_zero_volatility},
        {"monte_carlo_one_step", test_monte_carlo_one_step},
        {"invalid_inputs", test_invalid_inputs},
        {"statistics_raw", test_statistics_raw},
        {"statistics_split_merge", test_statistics_split_merge},
        {"statistics_unequal_batches", test_statistics_unequal_batches},
        {"statistics_empty_singleton", test_statistics_empty_singleton},
        {"statistics_discounting", test_statistics_discounting},
        {"single_threaded_regression", test_single_threaded_regression},
        {"statistics_count_overflow", test_statistics_count_overflow}
    }};
    if (argc > 2) {
        std::cerr << "Usage: test_pricing [test_name]\n";
        return EXIT_FAILURE;
    }
    bool matched = false;
    bool failed = false;
    for (const auto& test : tests) {
        if (argc == 2 && test.name != argv[1]) {
            continue;
        }
        matched = true;
        try {
            test.run();
            std::cout << "PASS: " << test.name << '\n';
        } catch (const std::exception& error) {
            std::cerr << "FAIL: " << test.name << ": " << error.what() << '\n';
            failed = true;
        } catch (...) {
            std::cerr << "FAIL: " << test.name << ": unknown exception\n";
            failed = true;
        }
    }
    if (!matched) {
        std::cerr << "Unknown test name: " << argv[1] << '\n';
        return EXIT_FAILURE;
    }
    return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
