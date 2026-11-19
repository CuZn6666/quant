#include "visualization.hpp"

#include "gbm.hpp"

#include <cstddef>
#include <random>
#include <stdexcept>
#include <vector>

#include <matplot/matplot.h>

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
