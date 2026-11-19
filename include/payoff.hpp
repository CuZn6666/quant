#pragma once

#include <vector>

double Asian_call_payoff(const std::vector<double>& path, double k);
double Asian_put_payoff(const std::vector<double>& path, double k);
double Lookback_call_payoff(const std::vector<double>& path, double k);
double Lookback_put_payoff(const std::vector<double>& path, double k);
