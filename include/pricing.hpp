#pragma once

#include <cstddef>
#include <random>
#include<cstdint>//std::uint64_t Specify the random number seed explicitly

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
    double sum;
    double sum_sq;
    PayoffStatistics(double total = 0.0, double squares = 0.0)
        : sum(total), sum_sq(squares) {}
};

struct MonteCarloStatistics
{
    std::size_t n;
    PayoffStatistics asian_call;
    PayoffStatistics asian_put;
    PayoffStatistics lookback_call;
    PayoffStatistics lookback_put;
    MonteCarloStatistics(std::size_t count = 0,
        PayoffStatistics ac = {}, PayoffStatistics ap = {},
        PayoffStatistics lc = {}, PayoffStatistics lp = {})
        : n(count), asian_call(ac), asian_put(ap), lookback_call(lc), lookback_put(lp) {}
};

MonteCarloStatistics accumulate_monte_carlo(double initial_price, double volatility,
    double rate, double maturity, double strike, std::size_t steps,
    std::size_t paths, std::mt19937_64& rng);

MonteCarloStatistics merge_statistics(const MonteCarloStatistics& left,
    const MonteCarloStatistics& right);

OptionPrices finalize_monte_carlo(const MonteCarloStatistics& statistics,
    double rate, double maturity);

OptionPrices monte_carlo_pricing(double s0, double sigma, double r, double T, double K, std::size_t steps, std::size_t M, std::mt19937_64& rng);

OptionPrices monte_carlo_pricing_parallel(double s0,double sigma,double r,double T,double K,std::size_t steps,std::size_t M,std::size_t num_threads,std::uint64_t base_seed);
// Exposes merged raw moments for auditing; like parallel pricing, requires M >= 2.
MonteCarloStatistics accumulate_monte_carlo_parallel(double s0, double sigma, double r,
    double T, double K, std::size_t steps, std::size_t M,
    std::size_t num_threads, std::uint64_t base_seed);
//Each thread should have its own RNG to avoid race conditions.
//num_threads: The number of worker threads to use.
//base_seed: The starting point for assigning random number seeds to all workers.
