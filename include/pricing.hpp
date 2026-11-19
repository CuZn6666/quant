#pragma once

#include <cstddef>
#include <random>

struct OptionPrices
{
    //add put,before only has call
    double asian_call;//V Estimate
    double asian_put;
    double lookback_call;
    double lookback_put;
    double asian_call_std_error;//SE
    double asian_put_std_error;
    double lookback_call_std_error;//A measure of the statistical accuracy of simulation results.
    double lookback_put_std_error;
};

struct PayoffStatistics
{
    double sum = 0.0;
    double sum_sq = 0.0;
};

struct MonteCarloStatistics
{
    std::size_t n = 0;
    PayoffStatistics asian_call;
    PayoffStatistics asian_put;
    PayoffStatistics lookback_call;
    PayoffStatistics lookback_put;
};

MonteCarloStatistics accumulate_monte_carlo(double initial_price, double volatility,
    double rate, double maturity, double strike, std::size_t steps,
    std::size_t paths, std::mt19937_64& rng);

MonteCarloStatistics merge_statistics(const MonteCarloStatistics& left,
    const MonteCarloStatistics& right);

OptionPrices finalize_monte_carlo(const MonteCarloStatistics& statistics,
    double rate, double maturity);

OptionPrices monte_carlo_pricing(double s0, double sigma, double r, double T, double K, std::size_t steps, std::size_t M, std::mt19937_64& rng);
