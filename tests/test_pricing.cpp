#include "gbm.hpp"
#include "payoff.hpp"
#include "pricing.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <initializer_list>
#include <iostream>
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

struct TestCase {
    std::string_view name;
    void (*run)();
};

}

int main(int argc, char* argv[])
{
    const std::array<TestCase, 9> tests = {{
        {"asian_call", test_asian_call},
        {"asian_put", test_asian_put},
        {"lookback_call", test_lookback_call},
        {"lookback_put", test_lookback_put},
        {"gbm_properties", test_gbm_properties},
        {"monte_carlo_properties", test_monte_carlo_properties},
        {"monte_carlo_zero_volatility", test_monte_carlo_zero_volatility},
        {"monte_carlo_one_step", test_monte_carlo_one_step},
        {"invalid_inputs", test_invalid_inputs}
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
