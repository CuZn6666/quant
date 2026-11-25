#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

// Shared by the path generator and batch accumulator (including empty batches).
inline void validate_gbm(double s0, double sigma, double rate,
                         std::size_t steps, double maturity)
{
    if (!std::isfinite(s0) || s0 <= 0.0 || !std::isfinite(sigma) || sigma < 0.0 ||
        !std::isfinite(rate) || !std::isfinite(maturity) || maturity < 0.0) {
        throw std::invalid_argument("GBM requires finite S0 > 0, sigma >= 0, rate, T >= 0");
    }
    if (steps == 0) {
        throw std::invalid_argument("steps must be positive");
    }
    if (steps == std::numeric_limits<std::size_t>::max()) {
        throw std::overflow_error("GBM steps + 1 overflows std::size_t");
    }
    if (steps + 1 > std::vector<double>().max_size()) {
        throw std::length_error("GBM path exceeds vector capacity");
    }
}
