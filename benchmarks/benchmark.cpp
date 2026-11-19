#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <cstddef>//A time-tracking tool.

#include "pricing.hpp"

int main()
{
    const double S0 = 100.0;
    const double sigma = 0.20;
    const double r = 0.05;
    const double T = 1.0;
    const double K = 100.0;
    std::size_t steps = 252;
    // Different numbers of Monte Carlo paths to test
    std::vector<std::size_t> M_values = {1000,10000,100000};
    std::size_t repeats = 3;
    std::cout<< "M,Average Time (seconds),"<< "Asian Call,Asian Put,"<< "Lookback Call,Lookback Put\n";
    // ------ Benchmark Loop --------
    for(std::size_t i : M_values)//i is the number of Monte Carlo paths
    {
        double total_time  = 0.0;
        // Accumulate option prices across repeated runs
        OptionPrices price_sum = {};
        for(std::size_t j=0;j<repeats;j++)//j is the number of the replicate experiment.
        {
            std::mt19937_64 rng(2023 + j);//// Independent RNG for each run
            auto start = std::chrono::steady_clock::now();//start time
            OptionPrices price = monte_carlo_pricing(S0,sigma,r,T,K,steps,i,rng);
            auto end = std::chrono::steady_clock::now();//end time
            std::chrono::duration<double> elapsed = end - start;
            total_time += elapsed.count();
            price_sum.asian_call += price.asian_call;
            price_sum.asian_put += price.asian_put;

            price_sum.lookback_call += price.lookback_call;
            price_sum.lookback_put += price.lookback_put;
        }//debug --// Inner loop ends here!
        // ------ Calculate Averages --------
            const double average_time = total_time / repeats;
            const double asian_call_avg = price_sum.asian_call / repeats;
            const double asian_put_avg = price_sum.asian_put / repeats;
            const double lookback_call_avg = price_sum.lookback_call / repeats;
            const double lookback_put_avg = price_sum.lookback_put / repeats;           
        
        // ------- Output Results --------
            std::cout<< i << ","<< average_time << ","<< asian_call_avg << ","<< asian_put_avg << ","<< lookback_call_avg << ","<< lookback_put_avg << std::endl;
    }
}
