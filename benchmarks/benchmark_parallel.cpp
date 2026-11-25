#include "pricing.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double S0 = 100.0, SIGMA = .20, R = .05, T = 1.0, K = 100.0;
constexpr std::uint64_t BASE_SEED = 2023;
const std::array<std::size_t, 5> THREADS = {{0, 1, 2, 4, 8}}; // 0 = original baseline
volatile double result_sink = 0;
struct Settings {
    std::size_t paths = 1000000, steps = 252, repeats = 5;
    std::string csv, raw_csv;
};
struct Measurement { std::size_t configuration, repeat, order; double ms; OptionPrices prices; };

std::size_t number(const std::string& text)
{
    if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
        throw std::invalid_argument("Expected a positive integer: " + text);
    const unsigned long long parsed = std::stoull(text);
    if (parsed == 0 || parsed > std::numeric_limits<std::size_t>::max())
        throw std::invalid_argument("Integer out of range: " + text);
    return static_cast<std::size_t>(parsed);
}
void consume(const OptionPrices& p)
{
    const double fields[] = {p.asian_call, p.asian_put, p.lookback_call, p.lookback_put,
        p.asian_call_std_error, p.asian_put_std_error, p.lookback_call_std_error, p.lookback_put_std_error};
    double sum = 0;
    for (double value : fields) {
        if (!std::isfinite(value) || value < 0) throw std::runtime_error("Invalid benchmark result");
        sum += value;
    }
    result_sink = sum;
}
Measurement measure(const Settings& settings, std::size_t config, std::size_t repeat, std::size_t order)
{
    // Include RNG initialization in both APIs. Parallel timing also includes launch,
    // scheduling, joins and merging: this measures end-to-end pricing latency.
    const auto start = std::chrono::steady_clock::now();
    OptionPrices p;
    if (config == 0) {
        std::mt19937_64 rng(BASE_SEED);
        p = monte_carlo_pricing(S0, SIGMA, R, T, K, settings.steps, settings.paths, rng);
    } else {
        p = monte_carlo_pricing_parallel(S0, SIGMA, R, T, K, settings.steps, settings.paths, THREADS[config], BASE_SEED);
    }
    const auto end = std::chrono::steady_clock::now();
    consume(p); // Validation / output / vector storage are outside the timed region.
    return {config, repeat, order, std::chrono::duration<double, std::milli>(end - start).count(), p};
}
double median(std::vector<double> times)
{
    std::sort(times.begin(), times.end());
    const std::size_t n = times.size();
    return n % 2 ? times[n / 2] : (times[n / 2 - 1] + times[n / 2]) / 2;
}
const char* name(std::size_t config) { return config == 0 ? "original_single" : "parallel"; }
void metadata(std::ostream& out, const Settings& settings)
{
    out << ',' << settings.paths << ',' << settings.steps << ',' << settings.repeats
        << ",1," << S0 << ',' << SIGMA << ',' << R << ',' << T << ',' << K << ',' << BASE_SEED
        << ',' << QUANT_BUILD_CONFIG << '\n';
}
const char* metadata_header = ",paths,steps,repeats,warmups,s0,sigma,rate,maturity,strike,base_seed,build\n";
void summary(std::ostream& out, const Settings& settings, const std::vector<Measurement>& measurements)
{
    std::array<std::vector<double>, 5> times;
    for (const auto& m : measurements) times[m.configuration].push_back(m.ms);
    const double baseline = median(times[0]);
    out << std::setprecision(12) << "implementation,threads,median_ms,min_ms,max_ms,speedup,efficiency" << metadata_header;
    for (std::size_t config = 0; config < THREADS.size(); ++config) {
        const std::size_t threads = config == 0 ? 1 : THREADS[config];
        const double elapsed = median(times[config]);
        const auto bounds = std::minmax_element(times[config].begin(), times[config].end());
        out << name(config) << ',' << threads << ',' << elapsed << ',' << *bounds.first << ',' << *bounds.second
            << ',' << baseline / elapsed << ',' << baseline / elapsed / threads;
        metadata(out, settings);
    }
}
void raw(std::ostream& out, const Settings& settings, const std::vector<Measurement>& measurements)
{
    out << std::setprecision(17) << "implementation,threads,repeat,order,elapsed_ms,asian_call,asian_put,lookback_call,lookback_put,asian_call_se,asian_put_se,lookback_call_se,lookback_put_se" << metadata_header;
    for (const auto& m : measurements) {
        const auto& p = m.prices;
        out << name(m.configuration) << ',' << (m.configuration == 0 ? 1 : THREADS[m.configuration])
            << ',' << m.repeat << ',' << m.order << ',' << m.ms << ',' << p.asian_call << ',' << p.asian_put
            << ',' << p.lookback_call << ',' << p.lookback_put << ',' << p.asian_call_std_error
            << ',' << p.asian_put_std_error << ',' << p.lookback_call_std_error << ',' << p.lookback_put_std_error;
        metadata(out, settings);
    }
}
}

int main(int argc, char** argv)
{
    try {
        Settings settings;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--help") {
                std::cout << "Usage: benchmark_parallel [--paths N] [--steps N] [--repeats N] [--csv FILE] [--raw-csv FILE]\n"
                          << "Defaults: 1000000 paths, 252 steps, 5 timed repetitions, 1 warmup/configuration. Release only.\n";
                return 0;
            }
            if (i + 1 == argc) throw std::invalid_argument("Missing value for " + arg);
            const std::string value = argv[++i];
            if (arg == "--paths") settings.paths = number(value);
            else if (arg == "--steps") settings.steps = number(value);
            else if (arg == "--repeats") settings.repeats = number(value);
            else if (arg == "--csv") settings.csv = value;
            else if (arg == "--raw-csv") settings.raw_csv = value;
            else throw std::invalid_argument("Unknown option: " + arg);
        }
        if (std::string(QUANT_BUILD_CONFIG) != "Release") throw std::runtime_error("Configure with -DCMAKE_BUILD_TYPE=Release (or build --config Release)");
        if (settings.paths < 8) throw std::invalid_argument("Use at least 8 paths to exercise all requested workers");
        if (settings.steps >= std::vector<double>().max_size()) throw std::invalid_argument("steps exceeds path capacity");
        if (settings.repeats > std::vector<Measurement>().max_size() / THREADS.size()) throw std::invalid_argument("Too many repetitions");
        if (!settings.csv.empty() && settings.csv == settings.raw_csv) throw std::invalid_argument("CSV destinations must differ");
        std::ofstream csv, raw_csv;
        if (!settings.csv.empty()) { csv.open(settings.csv.c_str()); if (!csv) throw std::runtime_error("Cannot open summary CSV"); }
        if (!settings.raw_csv.empty()) { raw_csv.open(settings.raw_csv.c_str()); if (!raw_csv) throw std::runtime_error("Cannot open raw CSV"); }
        std::vector<Measurement> measurements;
        measurements.reserve(settings.repeats * THREADS.size());
        std::cerr << "Release benchmark: " << settings.paths << " paths, " << settings.steps << " steps, "
                  << settings.repeats << " timed repetitions, 1 warmup per configuration.\n";
        for (std::size_t config = 0; config < THREADS.size(); ++config) {
            measure(settings, config, 0, config);
            std::cerr << "Warmup complete: " << name(config) << " threads=" << (config ? THREADS[config] : 1) << '\n';
        }
        // Rotate configurations each round: with five repeats each occupies each position once.
        for (std::size_t repeat = 0; repeat < settings.repeats; ++repeat) {
            for (std::size_t order = 0; order < THREADS.size(); ++order) {
                const std::size_t config = (repeat % THREADS.size() + order) % THREADS.size();
                measurements.push_back(measure(settings, config, repeat + 1, order + 1));
            }
            std::cerr << "Completed repetition " << repeat + 1 << '/' << settings.repeats << '\n';
        }
        summary(std::cout, settings, measurements);
        if (csv.is_open()) { summary(csv, settings, measurements); csv.flush(); if (!csv) throw std::runtime_error("Summary CSV write failed"); }
        if (raw_csv.is_open()) { raw(raw_csv, settings, measurements); raw_csv.flush(); if (!raw_csv) throw std::runtime_error("Raw CSV write failed"); }
        if (!std::cout) throw std::runtime_error("Standard output write failed");
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Benchmark failed: " << e.what() << '\n';
        return 1;
    }
}
