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

OptionPrices monte_carlo_pricing(double s0, double sigma, double r, double T, double K, std::size_t steps, std::size_t M, std::mt19937_64& rng);
