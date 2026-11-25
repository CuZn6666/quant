#include "pricing.hpp"

#include "gbm.hpp"
#include "gbm_validation.hpp"
#include "payoff.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

MonteCarloStatistics accumulate_monte_carlo(double initial_price, double volatility,double rate, double maturity, double strike, std::size_t steps,std::size_t paths, std::mt19937_64& rng)
{
    validate_gbm(initial_price, volatility, rate, steps, maturity);
    if (!std::isfinite(strike) || strike < 0.0) {
        throw std::invalid_argument("strike must be finite and nonnegative");
    }

    MonteCarloStatistics statistics{};
    for (std::size_t index = 0; index < paths; ++index) {
        const std::vector<double> path =
            simulation_GBM(initial_price, volatility, rate, steps, maturity, rng);
        const double asian_call = Asian_call_payoff(path, strike);
        const double asian_put = Asian_put_payoff(path, strike);
        const double lookback_call = Lookback_call_payoff(path, strike);
        const double lookback_put = Lookback_put_payoff(path, strike);

        statistics.asian_call.sum += asian_call;
        statistics.asian_put.sum += asian_put;
        statistics.lookback_call.sum += lookback_call;
        statistics.lookback_put.sum += lookback_put;

        statistics.asian_call.sum_sq += asian_call * asian_call;
        statistics.asian_put.sum_sq += asian_put * asian_put;
        statistics.lookback_call.sum_sq += lookback_call * lookback_call;
        statistics.lookback_put.sum_sq += lookback_put * lookback_put;
    }
    statistics.n = paths;
    return statistics;
}

MonteCarloStatistics merge_statistics(const MonteCarloStatistics& left,
    const MonteCarloStatistics& right)
{
    if (right.n > std::numeric_limits<std::size_t>::max() - left.n) {
        throw std::overflow_error("Merged sample count exceeds std::size_t");
    }

    MonteCarloStatistics result{};
    result.n = left.n + right.n;
    result.asian_call.sum = left.asian_call.sum + right.asian_call.sum;
    result.asian_put.sum = left.asian_put.sum + right.asian_put.sum;
    result.lookback_call.sum = left.lookback_call.sum + right.lookback_call.sum;
    result.lookback_put.sum = left.lookback_put.sum + right.lookback_put.sum;
    result.asian_call.sum_sq = left.asian_call.sum_sq + right.asian_call.sum_sq;
    result.asian_put.sum_sq = left.asian_put.sum_sq + right.asian_put.sum_sq;
    result.lookback_call.sum_sq = left.lookback_call.sum_sq + right.lookback_call.sum_sq;
    result.lookback_put.sum_sq = left.lookback_put.sum_sq + right.lookback_put.sum_sq;
    return result;
}

OptionPrices finalize_monte_carlo(const MonteCarloStatistics& statistics,
    double rate, double maturity)
{
    if (!std::isfinite(rate) || !std::isfinite(maturity) || maturity < 0.0) {
        throw std::invalid_argument("Finalization requires finite rate and T >= 0");
    }
    if (statistics.n < 2) {
        throw std::invalid_argument("At least two samples are required to calculate sample variance");
    }

    const std::size_t samples = statistics.n;
    const double asian_call_mean = statistics.asian_call.sum / samples;
    const double asian_put_mean = statistics.asian_put.sum / samples;
    const double lookback_call_mean = statistics.lookback_call.sum / samples;
    const double lookback_put_mean = statistics.lookback_put.sum / samples;

    const double asian_call_variance = std::max(0.0,
        (statistics.asian_call.sum_sq - samples * asian_call_mean * asian_call_mean) / (samples - 1));
    const double asian_put_variance = std::max(0.0,
        (statistics.asian_put.sum_sq - samples * asian_put_mean * asian_put_mean) / (samples - 1));
    const double lookback_call_variance = std::max(0.0,
        (statistics.lookback_call.sum_sq - samples * lookback_call_mean * lookback_call_mean) / (samples - 1));
    const double lookback_put_variance = std::max(0.0,
        (statistics.lookback_put.sum_sq - samples * lookback_put_mean * lookback_put_mean) / (samples - 1));

    const double discount = std::exp(-rate * maturity);
    OptionPrices result{};
    result.asian_call = discount * asian_call_mean;
    result.asian_put = discount * asian_put_mean;
    result.lookback_call = discount * lookback_call_mean;
    result.lookback_put = discount * lookback_put_mean;
    result.asian_call_std_error = discount * std::sqrt(asian_call_variance / samples);
    result.asian_put_std_error = discount * std::sqrt(asian_put_variance / samples);
    result.lookback_call_std_error = discount * std::sqrt(lookback_call_variance / samples);
    result.lookback_put_std_error = discount * std::sqrt(lookback_put_variance / samples);
    return result;
}

OptionPrices monte_carlo_pricing(double s0, double sigma, double r, double T, double K, std::size_t steps, std::size_t M, std::mt19937_64& rng)
{
    if (M < 2) 
        throw std::invalid_argument(" M must be at least 2 to calculate sample variance");
    return finalize_monte_carlo(accumulate_monte_carlo(s0, sigma, r, T, K, steps, M, rng), r, T);
}
