#include "pricing.hpp"

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

namespace {
void print_result(const std::string& implementation, const OptionPrices& p)
{
    const char* names[] = {"Asian call", "Asian put", "Lookback call", "Lookback put"};
    const double prices[] = {p.asian_call, p.asian_put, p.lookback_call, p.lookback_put};
    const double errors[] = {p.asian_call_std_error, p.asian_put_std_error,
        p.lookback_call_std_error, p.lookback_put_std_error};
    for (std::size_t i = 0; i < 4; ++i) {
        if (!std::isfinite(prices[i]) || !std::isfinite(errors[i]) || prices[i] < 0 || errors[i] < 0)
            throw std::runtime_error("Invalid pricing result");
        std::cout << std::left << std::setw(17) << implementation << std::setw(16) << names[i]
                  << std::right << std::setw(12) << prices[i] << std::setw(14) << errors[i] << '\n';
    }
}
}

int main()
{
    try {
        // Deliberately small presentation workload; this is not a timing benchmark.
        const std::size_t paths = 10000, steps = 252;
        const unsigned long long seed = 2023;
        std::cout << "Monte Carlo pricing demo | 10000 paths | 252 steps | seed 2023\n"
                  << "S0=100, K=100, sigma=0.20, r=0.05, T=1 year\n"
                  << "Each estimate uses all 10000 paths. SE = discounted Monte Carlo standard error.\n\n"
                  << std::left << std::setw(17) << "Implementation" << std::setw(16) << "Option"
                  << std::right << std::setw(12) << "Price" << std::setw(14) << "SE" << '\n'
                  << std::fixed << std::setprecision(6);
        std::mt19937_64 rng(seed);
        print_result("Original single", monte_carlo_pricing(100, .2, .05, 1, 100, steps, paths, rng));
        for (std::size_t threads : {1, 2, 4, 8}) {
            print_result("Parallel " + std::to_string(threads), monte_carlo_pricing_parallel(
                100, .2, .05, 1, 100, steps, paths, threads, seed));
        }
        std::cout << "\nOne worker shares the serial stream. More workers change streams;\n"
                  << "cross-thread-count estimates need not be identical. See CTest for statistical checks.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Demo failed: " << error.what() << '\n';
        return 1;
    }
}
