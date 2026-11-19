#include "pricing.hpp"

#include "gbm.hpp"
#include "payoff.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include <stdexcept>
#include <vector>

OptionPrices monte_carlo_pricing(double s0,double sigma,double r,double T,double K,std::size_t steps,std::size_t M,std::mt19937_64& rng)
{
    if(M < 2 )// M !=1 because M-1 != 0 for variance_caculate
        throw std::invalid_argument(" M must be at least 2 to calculate sample variance");
    /*
    double asian_sum = 0.0;
    double lookback_sum = 0.0;
    double asian_sq_sum = 0.0;//Add a sum-of-squares variable
    double lookback_sq_sum = 0.0;
    for(std::size_t i=0;i<M;i++)//get "SUM" of monte carlo
    {
        // Generate one GBM path
        std::vector<double> path = simulation_GBM(s0,sigma,r,steps,T,rng);
        // Calculate both payoffs
        //asian_sum += Asian_call_payoff(path,K);
        //lookback_sum += Lookback_call_payoff(path,K);
       const double asian_payoff = Asian_call_payoff(path,K);
       const double lookback_payoff = Lookback_call_payoff(path,K);
       // To get Sum(pi)
       asian_sum += asian_payoff;
       lookback_sum += lookback_payoff;
       // To get Sum(pi^2)
       asian_sq_sum += asian_payoff * asian_payoff;
       lookback_sq_sum += lookback_payoff * lookback_payoff;
    }
    //------4th--version------add standard error

    const double Exponent = std::exp(-r*T);// // Discount factor
    //mean
    const double asian_mean = asian_sum / M;
    const double lookback_mean = lookback_sum / M;
    //variance
    //Important Metho(std::max)d!!!! Handling Minor Negative Variances Caused by Floating-Point Errors
    const double asian_variance = std::max(0.0,(asian_sq_sum - M * asian_mean * asian_mean) / (M -1));
    const double lookback_variance = std::max(0.0,(lookback_sq_sum - M * lookback_mean * lookback_mean) / (M -1));

    OptionPrices result = {};// initialize all of its members to zero.
    //---5th--version--update and add formular
    //result.asian_call = Exponent * asian_sum / M;
    //result.lookback_call = Exponent * lookback_sum / M;
    result.asian_call = Exponent * asian_mean;//update
    result.lookback_call = Exponent * lookback_mean;
    result.asian_call_std_error = std::sqrt(asian_variance / M) * Exponent;
    result.lookback_call_std_error = std::sqrt(lookback_variance / M) * Exponent;
    return result;
    */
// update new version , add put before only has call
// -------- Sum of Payoffs --------
    double asian_call_sum = 0.0;
    double asian_put_sum = 0.0;
    double lookback_call_sum = 0.0;
    double lookback_put_sum = 0.0;
// -------- Sum of Squared Payoffs --------
    double asian_call_sq_sum = 0.0;
    double asian_put_sq_sum = 0.0;
    double lookback_call_sq_sum = 0.0;
    double lookback_put_sq_sum = 0.0;
    for (std::size_t i = 0; i < M; ++i)
    {
    // Generate one GBM path
        std::vector<double> path =
        simulation_GBM(s0, sigma, r, steps, T, rng);
    // Calculate four payoffs
        const double AC = Asian_call_payoff(path, K);//Asian_Call
        const double AP = Asian_put_payoff(path, K);//Put
        const double LC = Lookback_call_payoff(path, K);
        const double LP = Lookback_put_payoff(path, K);

    // -------- Sum of Payoffs --------
        asian_call_sum += AC;
        asian_put_sum += AP;
        lookback_call_sum += LC;
        lookback_put_sum += LP;
    // -------- Sum of Squared Payoffs --------
        asian_call_sq_sum += AC * AC;
        asian_put_sq_sum += AP * AP;
        lookback_call_sq_sum += LC * LC;
        lookback_put_sq_sum += LP * LP;
    }
// -------- Mean Payoffs --------
    const double asian_call_mean = asian_call_sum / M;
    const double asian_put_mean = asian_put_sum / M;
    const double lookback_call_mean = lookback_call_sum / M;
    const double lookback_put_mean = lookback_put_sum / M;
// -------- Sample Variances --------
    const double asian_call_variance = std::max(0.0,(asian_call_sq_sum - M * asian_call_mean * asian_call_mean)/ (M - 1));
    const double asian_put_variance = std::max(0.0,(asian_put_sq_sum - M * asian_put_mean * asian_put_mean)/ (M - 1));
    const double lookback_call_variance = std::max(0.0,(lookback_call_sq_sum - M * lookback_call_mean * lookback_call_mean)/ (M - 1));
    const double lookback_put_variance = std::max(0.0,(lookback_put_sq_sum - M * lookback_put_mean * lookback_put_mean)/ (M - 1));
// -------- Discount Factor --------
    const double discount = std::exp(-r * T);
    OptionPrices result{};  
// -------- Option Prices --------
    result.asian_call = discount * asian_call_mean;
    result.asian_put = discount * asian_put_mean;
    result.lookback_call = discount * lookback_call_mean;
    result.lookback_put = discount * lookback_put_mean;
// -------- Standard Errors --------
    result.asian_call_std_error = discount * std::sqrt(asian_call_variance / M);
    result.asian_put_std_error = discount * std::sqrt(asian_put_variance / M);
    result.lookback_call_std_error = discount * std::sqrt(lookback_call_variance / M);
    result.lookback_put_std_error = discount * std::sqrt(lookback_put_variance / M);
    return result;
}
