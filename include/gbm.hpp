#pragma once

#include <cstddef>
#include <random>
#include <vector>

std::vector<double> simulation_GBM(double s0, double sigma, double r, std::size_t steps, double T, std::mt19937_64& rng);
