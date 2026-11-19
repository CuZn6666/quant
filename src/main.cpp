//-------------first---version-----------------
//Single-Threaded Monte Carlo Simulation
#include<iostream>
#include<vector>
#include<cmath>
#include<random>//Random Number Engines & Probability Distributions
#include<cstddef>// size_t--> unsigned
#include<algorithm>
#include<numeric>//Sigma
#include <stdexcept>//Handling Exceptions
#include <fstream>// to get visualization path
#include <matplot/matplot.h>//to get plot
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

std::vector<double> simulation_GBM(double s0,double sigma,double r,std::size_t steps,double T,std::mt19937_64 & rng)//rng =: Random Number Generator
//Exact Discretization not Euler–Maruyama!
{
    if (steps == 0)
    throw std::invalid_argument("steps must be positive");
    const double dt = T / steps;//deta t
    std::normal_distribution<double> normal(0.0,1.0);
    std::vector<double> S(steps+1);//t0----tn
    S[0] = s0;
    for(size_t i=0;i<steps;i++)//debug i<step rather than r<S.size()
    {
        const double Z = normal(rng);//Generate a random number Z that follows a normal distribution.
        const double dirft = (r - 0.5*sigma*sigma)*dt;
        const double diffusion = sigma * std::sqrt(dt) * Z;
        S[i+1] = S[i] * std::exp(dirft + diffusion);
    }
    return S;// return an array
}

//--------------second---version-----------------
double Asian_call_payoff(const std::vector<double>& path,double k)
{
    //const double sum = std::accumulate(path.begin()+1,path.end(),0.0);
    //const double average = sum / (path.size() - 1);
    if(path.size() < 2)
        throw std::invalid_argument("Path must contain at least two prices");
    double sum = 0.0;
    for(std::size_t i=1;i<path.size();i++)//debug here i should start with 1 rather than 0
        sum += path[i];
    const double average = sum / (path.size() - 1);
    return std::max(average - k,0.0);
}

double Asian_put_payoff(const std::vector<double>& path,double k)
{
    if(path.size() < 2)
        throw std::invalid_argument("Path must contain at least two prices");
    double sum = 0.0;
    for(std::size_t i=1;i<path.size();i++)
        sum += path[i];
    const double average = sum / (path.size() - 1);
    return std::max(k - average,0.0);
}
//Fixed-Strike Lookback Option
double Lookback_call_payoff(const std::vector<double>& path,double k)
{
    //auto it = std::max_element(path.begin(),path.end());
    //double max_price = *it;//Here, * is a dereference, not a pointer.
    //Returns the position (iterator) of the largest element rather than the element itself
    if (path.empty())
        throw std::invalid_argument("Path must not be empty");
    double max_price = path[0];
    for(std::size_t i=0;i<path.size();i++)
    {
        if(path[i] > max_price)
            max_price = path[i];
    }
    return std::max(max_price - k,0.0);
}
double Lookback_put_payoff(const std::vector<double>& path,double k)
{
    if (path.empty())
        throw std::invalid_argument("Path must not be empty");
    double min_price = path[0];
    for(std::size_t i=0;i<path.size();i++)
    {
        if(path[i] < min_price)
            min_price = path[i];
    }
    return std::max(k - min_price,0.0);
}
//-----3rd----version-----
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
// -------- GBM Visualization --------
void polt_GBM_path(double s0,double sigma,double r,std::size_t steps,double T)
{
    if (steps == 0)
        throw std::invalid_argument("steps must be positive");
    const std::size_t N = 100;//Plot the top 100 GBM paths
    // Separate random number generator
    std::mt19937_64 polt_rng(2023);
    std::vector<double> time(steps+1);
    for(std::size_t i=0;i<=time.size();i++)
         time[i] = static_cast<double>(i);
    // Create figure
    matplot::figure(true);
    // Draw multiple paths on the same figure
    matplot::hold(matplot::on);
    for(std::size_t j=0;j<N;j++)
    {
        std::vector<double> path = simulation_GBM(s0,sigma,r,steps,T,polt_rng);
        matplot::plot(time,path);
    }
    //Image Coordinate Axes
    matplot::title("GBM Simulation - 100 Paths");
    matplot::xlabel("Time steps");
    matplot::ylabel("Stock Price");
    // Set X-axis range: 0 to 252
    matplot::xlim({0, static_cast<double>(steps)});
    //save image
    matplot::save("gbm_paths.png");
    //Display figure
    matplot::show();
}
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