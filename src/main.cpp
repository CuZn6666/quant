#include "gbm.hpp"
#include "payoff.hpp"
#include "pricing.hpp"
#include "visualization.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <random>
#include <vector>

int main()
{
    const double S0 = 100;
    const double r = 0.05;//risk_free_rate
    const double Sigma = 0.20;//volatility
    const double T = 1.0;
    const size_t steps = 252;//0---252 in total = 253
    std::mt19937_64 rng(2023);
    std::vector<double> path = simulation_GBM(S0,Sigma,r,steps,T,rng);
    std::cout<<"start price = "<<path.front()<<std::endl;
    std::cout<<"end price = "<<path.back()<<std::endl;
    // Calculate option payoffs
    double K = 100.0;//K is Strike Price（执行价格）

    double asian = Asian_call_payoff(path, K);
    double lookback = Lookback_call_payoff(path, K);

    std::cout << "Asian Call Payoff = " << asian << std::endl;
    std::cout << "Lookback Call Payoff = " << lookback << std::endl;

    // -------- Monte Carlo Pricing Test 1--------
    const std::size_t M = 10000;

// Separate RNG for Monte Carlo pricing
    std::mt19937_64 mc_rng(2023);
    //Call a Function
    OptionPrices prices = monte_carlo_pricing(S0, Sigma, r, T, K, steps, M, mc_rng);
    std::cout << "\n===== Monte Carlo Pricing =====\n";
    std::cout << "Number of paths = " << M << '\n';
    std::cout << "Asian Call Price = "<< prices.asian_call << '\n';
    std::cout << "Lookback Call Price = "<< prices.lookback_call << '\n';//or you can use std::endl
    std::cout << "Asian_Call Standard Error = "<<prices.asian_call_std_error<<std::endl;
    std::cout << "Lookback_CALL Standard Error ="<<prices.lookback_call_std_error<<std::endl;

// -------- Monte Carlo Convergence Test 2--------
//Test different numbers of simulation paths and observe the convergence of Monte Carlo pricing

    std::vector<std::size_t> test_M = {100,1000,10000,100000};
    std::cout << "\n===== Convergence Test =====\n";
    for (std::size_t i = 0; i < test_M.size(); i++)
    {
        std::mt19937_64 test_rng(2023);//The same random seed, 2023, is used here every time.
        OptionPrices result = monte_carlo_pricing(S0, Sigma, r, T, K, steps, test_M[i], test_rng);
        std::cout << "M = " << test_M[i]<< " | Asian_Call = " << result.asian_call<< " | Lookback_Call = " << result.lookback_call<< '\n';
    }
// -------- Monte Carlo Convergence Test 3--------<<one of interesting example>>
// 95% confidence interval [if Z ~ N(0,1)--->p(-1.96<=z<=1.96) almost = 0.95]
    std::cout << "\n===== 95% Confidence Interval Test =====\n";
    const std::size_t ci_M = 100000;
    std::mt19937_64  ci_rng(2023);
    OptionPrices ci_result = monte_carlo_pricing(S0,Sigma,r,T,K,steps,ci_M,ci_rng);
// Asian call
    std::cout <<"Asian Call Price"<<" "<<ci_result.asian_call<<std::endl;
    std::cout <<"Asian_Call Standard Error"<<" "<<ci_result.asian_call_std_error<<std::endl;
    std::cout <<"Asian_Call 95% Confidence interval"<<" "<<ci_result.asian_call - 1.96*ci_result.asian_call_std_error<<","<<ci_result.asian_call + 1.96 * ci_result.asian_call_std_error<<std::endl;
// Lookback call
    std::cout <<"Lookback Call Price"<<" "<<ci_result.lookback_call<<std::endl;
    std::cout <<"Lookback_Call Standard Error"<<" "<<ci_result.lookback_call_std_error<<std::endl;
    std::cout <<"Lookback_Call 95% Confidence interval"<<" "<<ci_result.lookback_call - 1.96*ci_result.lookback_call_std_error<<","<<ci_result.lookback_call + 1.96*ci_result.lookback_call_std_error<<std::endl;
// -------- Overall Payoff Test --------

    std::vector<double> test_path = {100, 105, 110, 95, 120};

    double test_K = 100.0;

    std::cout << "\n===== Four Payoff Test =====\n";
    std::cout << "Asian Call = "<< Asian_call_payoff(test_path, test_K) << '\n';
    std::cout << "Asian Put = "<< Asian_put_payoff(test_path, test_K) << '\n';
    std::cout << "Lookback Call = "<< Lookback_call_payoff(test_path, test_K) << '\n';
    std::cout << "Lookback Put = "<< Lookback_put_payoff(test_path, test_K) << '\n';
// -------- Monte Carlo Convergence Test 4---new test case after I add put 
    std::cout << "\n===== Monte Carlo Pricing =====\n";
    std::cout << "Asian Call Price = "<< prices.asian_call << '\n';
    std::cout << "Asian Put Price = "<< prices.asian_put << '\n';
    std::cout << "Lookback Call Price = "<< prices.lookback_call << '\n';
    std::cout << "Lookback Put Price = "<< prices.lookback_put << '\n';
    std::cout << "\n===== Standard Errors =====\n";
    std::cout << "Asian Call SE = "<< prices.asian_call_std_error << '\n';
    std::cout << "Asian Put SE = "<< prices.asian_put_std_error << '\n';
    std::cout << "Lookback Call SE = "<< prices.lookback_call_std_error << '\n';
    std::cout << "Lookback Put SE = "<< prices.lookback_put_std_error << '\n';
// -------- Fixed Path Test new --------
// Test four payoff functions using a deterministic stock price path
    std::cout << "\n===== Fixed Path Test =====\n";
    const std::vector<double> fixed_path = {100, 105, 110, 95, 120};
    const double fixed_K = 100.0;
// Calculate four payoffs
    const double AC = Asian_call_payoff(fixed_path, fixed_K);
    const double AP = Asian_put_payoff(fixed_path, fixed_K);
    const double LC = Lookback_call_payoff(fixed_path, fixed_K);
    const double LP = Lookback_put_payoff(fixed_path, fixed_K);
// Print results
    std::cout << "Asian Call Payoff = " << AC << '\n';
    std::cout << "Asian Put Payoff = " << AP << '\n';
    std::cout << "Lookback Call Payoff = " << LC << '\n';
    std::cout << "Lookback Put Payoff = " << LP << '\n';
// Verify results
    const double tolerance = 1e-10;
        bool passed = std::fabs(AC - 7.5) < tolerance && std::fabs(AP - 0.0) < tolerance && std::fabs(LC - 20.0) < tolerance && std::fabs(LP - 5.0) < tolerance;
    if (passed)
        std::cout << "Fixed Path Test: PASSED\n";
    else
        std::cout << "Fixed Path Test: FAILED\n";
//##################intersting examples
/*
This is because the maximum price of a stock's price path 
cannot be less than the average price of that path.
*/
    if (prices.lookback_call >= prices.asian_call)
        std::cout << "Payoff relation test: PASSED\n";
    else
        std::cout << "Payoff relation test: FAILED\n";
// -------- figures Test : GBM Visualization --------
    std::cout << "\n===== GBM Visualization =====\n";
    polt_GBM_path(S0, Sigma, r, steps, T);
    std::cout << "GBM Paths saved to gbm_paths.png\n";
    return 0;
}
