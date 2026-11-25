#include "pricing.hpp"

#include <algorithm>
#include <cstdint>//define std::uint64_t base_seed
#include <exception>//Saving and Passing Exceptions in Threads
#include <limits>//Check the value range and infinity
#include <random>//Generate a Monte Carlo random sample
#include <stdexcept>//Report Invalid Input
#include <thread>//Parallel Execution of Monte Carlo Simulations
#include <vector>

// --------------------------------------------------
// Parallel Monte Carlo Pricing
// --------------------------------------------------
MonteCarloStatistics accumulate_monte_carlo_parallel(double s0,double sigma,double r,double T,double K,std::size_t steps,std::size_t M,std::size_t num_threads,std::uint64_t base_seed)
{
    //input validation
    if(M < 2)
        throw std::invalid_argument("Monte Carlo requires at least 2 paths");
    if(num_threads == 0)
        throw std::invalid_argument("Number of threads must be positive");
    //Path partitioning
    const std::size_t actual_threads = std::min(num_threads,M);
    //The largest seed is: base_seed + (actual_threads - 1) || actual_threads - 1 meansMaximum thread number (numbered starting from 0)
    //--->base_seed + (actual_threads - 1) <= UNIT64_MAX = 2^64 - 1
    //--->exchange term --->actual_threads - 1 <= UNIT64_MAX - base_seed
    if(actual_threads - 1 > std::numeric_limits<std::uint64_t>::max() - base_seed)
    //std::numeric_limits<std::uint64_t>::max()
    //Get the largest integer that can be represented by the std::uint64_t type
        throw std::invalid_argument("Seed range overflow");
    const std::size_t base = M / actual_threads;
    const std::size_t remainder = M % actual_threads;
    std::vector<std::size_t> path_per_threads(actual_threads);
    for(std::size_t i=0;i<actual_threads;i++)
    {
        if(i < remainder)
            path_per_threads[i] = base + 1;
        else
            path_per_threads[i] = base;
    }
    //Prepare thread-local result storage
    std::vector<MonteCarloStatistics> partial_result(actual_threads);
    //Exceptions from individual workers
    std::vector<std::exception_ptr> errors(actual_threads);

    //-------------------------------------------
    //creat threads
    //--------------------------------------------
    std::vector<std::thread> workers;
    workers.reserve(actual_threads);//预留,避免重复线程开销
    try
    {
        for(std::size_t i=0;i<actual_threads;i++)
        {
            //Create a new element at the end of the array
            //&: Other external local variables used are captured by reference by default.
            //i: Captured by value; each thread maintains its own thread ID.
            workers.emplace_back([&,i]()
            {
                try// Each worker owns a separate RNG
                {
                    //Create a 64-bit pseudorandom number generator (rng) and use `base_seed + i` as the random seed.
                    std::mt19937_64 rng(base_seed + static_cast<std::uint64_t>(i));//static_cast --->Explicit Type Conversion
                    // Each worker calculates its own statistics
                    partial_result[i] = accumulate_monte_carlo(s0,sigma,r,T,K,steps,path_per_threads[i],rng);
                }
            catch(...)
            {
                // Save the exception for the main thread
                errors[i] = std::current_exception();
            }
            });
        }
    }
    catch(...)
    {
        // Thread creation may fail.
        // Join already-created threads before propagating.
        for(auto & worker : workers)//For each worker where 0 ≤ worker ≤ workers
        {
            if(worker.joinable())
                worker.join();
        }
        throw;
    }
    //----------------------------
    //join all threads
    //----------------------------   
    for(auto & worker : workers)
    {
        worker.join();
    }
    //Propagate worker exceptions after joining.
    for(const auto & error : errors)
    {
        if(error)
            std::rethrow_exception(error);
    }
    //-----------------------------
    //merge reults
    //-----------------------------
    MonteCarloStatistics total_statistics;
    for(std::size_t i=0;i<actual_threads;i++)
    {
        total_statistics  = merge_statistics(total_statistics,partial_result[i]);
    }
    //------------------------
    //finally
    //------------------------
    return total_statistics;
}

OptionPrices monte_carlo_pricing_parallel(double s0, double sigma, double r, double T,
    double K, std::size_t steps, std::size_t M, std::size_t num_threads,
    std::uint64_t base_seed)
{
    return finalize_monte_carlo(accumulate_monte_carlo_parallel(
        s0, sigma, r, T, K, steps, M, num_threads, base_seed), r, T);
}
